#!/bin/bash

# =====================================================
# Script: mapping_aa_to_cg.sh
# Description:
# - Runs cg-martini3 to generate CG mapping
# - Maps atomistic trajectory to coarse-grained (CG)
# - Generates topology, runs GROMACS preprocessing
# - Analyzes bonds/angles/dihedrals from CG trajectories
# =====================================================

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(dirname "$SCRIPT_DIR")"

# -----------------------
# Usage
# -----------------------
usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "ESSENTIAL ARGUMENTS:"
    echo "  --aa_tpr FILE             AA .tpr file"
    echo "  --aa_xtc FILE             AA trajectory"
    echo "  --aa_gro FILE             AA .gro file"
    echo "  --aa_itp FILE             AA .itp file"
    echo "  --beads_json FILE         Beads definition JSON file"
    echo "  --force_application STR   Force application (e.g., 'random=[1250,30;30,4]')"
    echo "  --beads_position STR      Beads position ('com' or 'geom')"
    echo "  --input_mdp FILE          MDP file for grompp"
    echo "  --path_ff DIR             Force field directory"
    echo ""
    echo "OPTIONAL ARGUMENTS:"
    echo "  --output_dir DIR          Output directory (default: results)"
    echo "  --name_molecule NAME      Molecule name (default: molecule)"
    echo "  --number_molecule NUM     Number of molecules (default: 1)"
    echo "  --ff FILE                 Force field ITP (default: martini_v3.0.0.itp)"
    echo "  --ions FILE               Ions ITP (default: martini_v3.0.0_ions_v1.itp)"
    echo "  --solvent FILE            Solvent ITP (default: martini_v3.0.0_solvents_v1.itp)"
    echo "  --title_comments TEXT     Topology comments"
    echo "  --title_system TEXT       System title"
    echo "  --output_topol FILE       Output topology (default: topol_cg.top)"
    echo "  --default_martini         Use default Martini 3 masses (72) and zero charges"
    echo "  --maxwarn N               Max warnings for grompp (default: 1)"
    echo "  --remove_pbc              Remove PBC (default: true)"
    echo "  --no_pbc                  Skip PBC removal"
    echo "  --skip_grompp             Skip grompp step"
    echo "  --skip_analysis           Skip bonds/angles/dihedrals analysis"
    echo "  --analyze_remove_pbc      Remove PBC before analysis"
    echo "  --keep_temp               Keep temporary files"
    echo "  --verbose                 Verbose output"
    echo ""
    echo "Example:"
    echo "  $0 --aa_tpr setup/md.tpr --aa_xtc setup/md.xtc --aa_gro setup/md.gro \\"
    echo "     --aa_itp setup/carb.itp --beads_json json/beads_config.json \\"
    echo "     --force_application 'random=[1250,30;30,4]' --beads_position com \\"
    echo "     --input_mdp mdp/minimization.mdp --path_ff ff_files/"
    exit 1
}

# -----------------------
# Default values
# -----------------------
# Essential args (no defaults)
AA_TPR=""
AA_XTC=""
AA_GRO=""
AA_ITP=""
BEADS_JSON=""
FORCE_APPLICATION=""
BEADS_POSITION=""
INPUT_MDP=""
PATH_FF=""

# Optional args with defaults
OUTPUT_DIR="results"
NAME_MOLECULE="molecule"
NUMBER_MOLECULE=1
FF="martini_v3.0.0.itp"
IONS="martini_v3.0.0_ions_v1.itp"
SOLVENT="martini_v3.0.0_solvents_v1.itp"
TITLE_COMMENTS="Topology system in Martini 3"
TITLE_SYSTEM="molecule in aqueous solution"
OUTPUT_TOPOL="topol_cg.top"
DEFAULT_MARTINI="false"
MAXWARN=1
REMOVE_PBC="true"
SKIP_GROMPP="false"
SKIP_ANALYSIS="false"
ANALYZE_REMOVE_PBC="false"
ANALYZE_GROUP_1="System"
ANALYZE_GROUP_2="System"
KEEP_INTERMEDIATE="false"
KEEP_TEMP="false"
VERBOSE="false"

# Internal variables
PYTHON_SCRIPTS_DIR=""
C_FILE=""
O_FILE=""
P_FILE=""

# -----------------------
# Parse arguments
# -----------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --aa_tpr) AA_TPR="$2"; shift 2 ;;
        --aa_xtc) AA_XTC="$2"; shift 2 ;;
        --aa_gro) AA_GRO="$2"; shift 2 ;;
        --aa_itp) AA_ITP="$2"; shift 2 ;;
        --beads_json) BEADS_JSON="$2"; shift 2 ;;
        --force_application) FORCE_APPLICATION="$2"; shift 2 ;;
        --beads_position) BEADS_POSITION="$2"; shift 2 ;;
        --input_mdp) INPUT_MDP="$2"; shift 2 ;;
        --path_ff) PATH_FF="$2"; shift 2 ;;
        
        --output_dir) OUTPUT_DIR="$2"; shift 2 ;;
        --name_molecule) NAME_MOLECULE="$2"; shift 2 ;;
        --number_molecule) NUMBER_MOLECULE="$2"; shift 2 ;;
        --ff) FF="$2"; shift 2 ;;
        --ions) IONS="$2"; shift 2 ;;
        --solvent) SOLVENT="$2"; shift 2 ;;
        --title_comments) TITLE_COMMENTS="$2"; shift 2 ;;
        --title_system) TITLE_SYSTEM="$2"; shift 2 ;;
        --output_topol) OUTPUT_TOPOL="$2"; shift 2 ;;
        --default_martini) DEFAULT_MARTINI="true"; shift ;;
        --maxwarn) MAXWARN="$2"; shift 2 ;;
        --remove_pbc) REMOVE_PBC="true"; shift ;;
        --no_pbc) REMOVE_PBC="false"; shift ;;
        --skip_grompp) SKIP_GROMPP="true"; shift ;;
        --skip_analysis) SKIP_ANALYSIS="true"; shift ;;
        --analyze_remove_pbc) ANALYZE_REMOVE_PBC="true"; shift ;;
        --analyze_group_1) ANALYZE_GROUP_1="$2"; shift 2 ;;
        --analyze_group_2) ANALYZE_GROUP_2="$2"; shift 2 ;;
        --keep_intermediate) KEEP_INTERMEDIATE="true"; shift ;;
        --keep_temp) KEEP_TEMP="true"; shift ;;
        --verbose) VERBOSE="true"; shift ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

# -----------------------
# Validate essential arguments
# -----------------------
if [ -z "$AA_TPR" ] || [ -z "$AA_XTC" ] || [ -z "$AA_GRO" ] || \
   [ -z "$AA_ITP" ] || [ -z "$BEADS_JSON" ] || [ -z "$FORCE_APPLICATION" ] || \
   [ -z "$BEADS_POSITION" ] || [ -z "$INPUT_MDP" ] || [ -z "$PATH_FF" ]; then
    echo "Error: Missing essential arguments"
    echo ""
    usage
fi

# -----------------------
# Auto-detect paths
# -----------------------
ORIGINAL_DIR="$(pwd)"

get_abs_path() {
    if [[ "$1" = /* ]]; then
        echo "$1"
    else
        echo "$ORIGINAL_DIR/$1"
    fi
}

AA_TPR_ABS=$(get_abs_path "$AA_TPR")
AA_XTC_ABS=$(get_abs_path "$AA_XTC")
AA_GRO_ABS=$(get_abs_path "$AA_GRO")
AA_ITP_ABS=$(get_abs_path "$AA_ITP")
BEADS_JSON_ABS=$(get_abs_path "$BEADS_JSON")
INPUT_MDP_ABS=$(get_abs_path "$INPUT_MDP")
PATH_FF_ABS=$(get_abs_path "$PATH_FF")

# Auto-detect python scripts directory
if [ -d "$PACKAGE_DIR/scripts" ]; then
    PYTHON_SCRIPTS_DIR="$PACKAGE_DIR/scripts"
elif [ -d "$SCRIPT_DIR/../scripts" ]; then
    PYTHON_SCRIPTS_DIR="$(cd "$SCRIPT_DIR/../scripts" && pwd)"
else
    echo "Error: Could not auto-detect python scripts directory"
    exit 1
fi

# -----------------------
# Create output directories
# -----------------------
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/GMX"
mkdir -p "$OUTPUT_DIR/ANALYSIS"
mkdir -p "$OUTPUT_DIR/NDX"

cd "$OUTPUT_DIR" || exit 1

# Function to check errors
check_error() {
    if [ $? -ne 0 ]; then
        echo "Error: $1 failed"
        cd "$ORIGINAL_DIR"
        exit 1
    fi
}

# Verbose output function
log_verbose() {
    if [ "$VERBOSE" = "true" ]; then
        echo "$1"
    fi
}

# -----------------------
# Display configuration
# -----------------------
echo "=========================================="
echo "CG Mapping Pipeline Configuration"
echo "=========================================="
echo "Essential arguments:"
echo "  AA TPR:            $AA_TPR_ABS"
echo "  AA XTC:            $AA_XTC_ABS"
echo "  AA GRO:            $AA_GRO_ABS"
echo "  AA ITP:            $AA_ITP_ABS"
echo "  Beads JSON:        $BEADS_JSON_ABS"
echo "  Force application: $FORCE_APPLICATION"
echo "  Beads position:    $BEADS_POSITION"
echo "  Input MDP:         $INPUT_MDP_ABS"
echo "  FF directory:      $PATH_FF_ABS"
echo ""
echo "Output directories:"
echo "  Main:              $OUTPUT_DIR"
echo "  GMX files:         $OUTPUT_DIR/GMX/"
echo "  Analysis:          $OUTPUT_DIR/ANALYSIS/"
echo "  NDX files:         $OUTPUT_DIR/NDX/"
echo ""
echo "Optional arguments:"
echo "  Molecule name:     $NAME_MOLECULE"
echo "  Number molecules:  $NUMBER_MOLECULE"
echo "  Remove PBC:        $REMOVE_PBC"
echo "  Skip grompp:       $SKIP_GROMPP"
echo "  Skip analysis:     $SKIP_ANALYSIS"
echo "  Analyze remove PBC:$ANALYZE_REMOVE_PBC"
echo "  Verbose:           $VERBOSE"
echo "=========================================="

# -----------------------
# Step 0: Run cg-martini3 to generate CG mapping
# -----------------------
echo ""
echo "Step 0/4: Running cg-martini3 to generate CG mapping"

# Check if cg-martini3 is available
if ! command -v cg-martini3 &> /dev/null; then
    echo "Error: cg-martini3 not found. Please install it first."
    echo "  pip install carb_param_automated_martini3"
    exit 1
fi

CG_MARTINI_CMD="cg-martini3 \
    --input_gro $AA_GRO_ABS \
    --input_beads_definitions $BEADS_JSON_ABS \
    --output_dir OUTPUT_CG \
    --beads_position $BEADS_POSITION \
    --aa_itp $AA_ITP_ABS \
    --force_application \"$FORCE_APPLICATION\""

if [ "$VERBOSE" = "true" ]; then
    CG_MARTINI_CMD="$CG_MARTINI_CMD --verbose"
fi

log_verbose "Running: $CG_MARTINI_CMD"
eval $CG_MARTINI_CMD
check_error "cg-martini3"

# Copy generated files to output directories
if [ -d "OUTPUT_CG/GRO" ]; then
    cp OUTPUT_CG/GRO/*.gro . 2>/dev/null
    cp OUTPUT_CG/GRO/*.gro GMX/ 2>/dev/null
fi

if [ -d "OUTPUT_CG/ITP" ]; then
    cp OUTPUT_CG/ITP/*.itp . 2>/dev/null
    cp OUTPUT_CG/ITP/*.itp GMX/ 2>/dev/null
fi

if [ -d "OUTPUT_CG/NDX" ]; then
    cp OUTPUT_CG/NDX/*.ndx NDX/ 2>/dev/null
    cp OUTPUT_CG/NDX/cg.ndx . 2>/dev/null
fi

if [ -d "OUTPUT_CG/JSON" ]; then
    cp OUTPUT_CG/JSON/*.json . 2>/dev/null
fi

# Set paths for generated files
CG_NDX="cg.ndx"
BONDS_NDX="NDX/bonds.ndx"
ANGLES_NDX="NDX/angles.ndx"
DIHEDRALS_NDX="NDX/dihedrals.ndx"
CG_GRO="cg.gro"
CG_ITP="cg.itp"

log_verbose "✓ cg-martini3 completed successfully"

# -----------------------
# Step 1: Mapping
# -----------------------
echo ""
echo "Step 1/4: Mapping AA → CG trajectory"

MAP_CMD="python $PYTHON_SCRIPTS_DIR/1-map_aa_traj_to_cg.py"
MAP_CMD_ARGS="--index_cg $CG_NDX --aa_tpr $AA_TPR_ABS --aa_xtc $AA_XTC_ABS --output_mapped mapped.xtc --output_cg_gro $CG_GRO"

if [ "$REMOVE_PBC" = "true" ]; then
    MAP_CMD_ARGS="$MAP_CMD_ARGS --remove_pbc"
else
    MAP_CMD_ARGS="$MAP_CMD_ARGS --no_pbc"
fi

if [ "$VERBOSE" = "true" ]; then
    MAP_CMD_ARGS="$MAP_CMD_ARGS --verbose"
fi

log_verbose "Running: $MAP_CMD $MAP_CMD_ARGS"
$MAP_CMD $MAP_CMD_ARGS
check_error "AA→CG mapping"

# Move mapping outputs to GMX directory
mv mapped.xtc GMX/ 2>/dev/null
mv $CG_GRO GMX/ 2>/dev/null
cp GMX/mapped.xtc . 2>/dev/null
cp GMX/$CG_GRO . 2>/dev/null

echo "  ✓ Generated: GMX/mapped.xtc"
echo "  ✓ Generated: GMX/$CG_GRO"

# -----------------------
# Step 2: Generate ITP
# -----------------------
echo ""
echo "Step 2/4: Generating CG ITP"

GENERATE_ITP_CMD="python $PYTHON_SCRIPTS_DIR/4-defining_atoms_type.py"
GENERATE_ITP_ARGS="--cg_gro GMX/$CG_GRO --cg_ndx $CG_NDX --aa_itp $AA_ITP_ABS --output $CG_ITP --name_molecule $NAME_MOLECULE"

# Find the mapping JSON file
DEF_JSON=""
if [ -f "$PACKAGE_DIR/../data/definitions_atoms_ff_martini3.json" ]; then
    DEF_JSON="$PACKAGE_DIR/../data/definitions_atoms_ff_martini3.json"
elif [ -f "$SCRIPT_DIR/../../data/definitions_atoms_ff_martini3.json" ]; then
    DEF_JSON="$SCRIPT_DIR/../../data/definitions_atoms_ff_martini3.json"
fi

if [ -n "$DEF_JSON" ]; then
    GENERATE_ITP_ARGS="$GENERATE_ITP_ARGS --def_json $DEF_JSON"
fi

if [ "$DEFAULT_MARTINI" = "true" ]; then
    GENERATE_ITP_ARGS="$GENERATE_ITP_ARGS --default"
fi

log_verbose "Running: $GENERATE_ITP_CMD $GENERATE_ITP_ARGS"
$GENERATE_ITP_CMD $GENERATE_ITP_ARGS
check_error "ITP generation"

# Move ITP to GMX directory
mv $CG_ITP GMX/ 2>/dev/null
cp GMX/$CG_ITP . 2>/dev/null

echo "  ✓ Generated: GMX/$CG_ITP"

# -----------------------
# Step 3: Generate topology and run grompp
# -----------------------
if [ "$SKIP_GROMPP" != "true" ]; then
    
    echo ""
    echo "Step 3/4: Generating topology and running grompp"
    
    # Build command for topology generation
    TOP_CMD="python $PYTHON_SCRIPTS_DIR/2-obtaining_cg_top.py"
    TOP_CMD_ARGS="--path_ff $PATH_FF_ABS"
    TOP_CMD_ARGS="$TOP_CMD_ARGS --ff $FF"
    TOP_CMD_ARGS="$TOP_CMD_ARGS --ions $IONS"
    TOP_CMD_ARGS="$TOP_CMD_ARGS --solvent $SOLVENT"
    TOP_CMD_ARGS="$TOP_CMD_ARGS --itp_ligand GMX/$CG_ITP"
    TOP_CMD_ARGS="$TOP_CMD_ARGS --name_molecule $NAME_MOLECULE"
    TOP_CMD_ARGS="$TOP_CMD_ARGS --number_molecule $NUMBER_MOLECULE"
    TOP_CMD_ARGS="$TOP_CMD_ARGS --title_comments \"$TITLE_COMMENTS\""
    TOP_CMD_ARGS="$TOP_CMD_ARGS --title_system \"$TITLE_SYSTEM\""
    TOP_CMD_ARGS="$TOP_CMD_ARGS --output_topol $OUTPUT_TOPOL"
    
    log_verbose "Running: $TOP_CMD $TOP_CMD_ARGS"
    eval $TOP_CMD $TOP_CMD_ARGS
    check_error "Topology generation"
    
    # Move topology to GMX directory
    mv $OUTPUT_TOPOL GMX/ 2>/dev/null
    cp GMX/$OUTPUT_TOPOL . 2>/dev/null
    echo "  ✓ Generated: GMX/$OUTPUT_TOPOL"
    
    # Set default files for grompp
    [ -z "$C_FILE" ] && C_FILE="GMX/$CG_GRO"
    [ -z "$P_FILE" ] && P_FILE="GMX/$OUTPUT_TOPOL"
    [ -z "$O_FILE" ] && O_FILE="GMX/CG.tpr"
    
    # Run grompp
    echo ""
    echo "Running grompp..."
    
    GROMPP_CMD="gmx grompp -f $INPUT_MDP_ABS -c $C_FILE -p $P_FILE -o $O_FILE -maxwarn $MAXWARN"
    
    if [ "$VERBOSE" = "true" ]; then
        GROMPP_CMD="$GROMPP_CMD -v"
    fi
    
    log_verbose "Running: $GROMPP_CMD"
    $GROMPP_CMD
    check_error "grompp"
    echo "  ✓ Generated: $O_FILE"
    
else
    echo ""
    echo "Step 3/4: Skipping grompp (--skip_grompp)"
fi

# -----------------------
# Step 4: Analyze bonds, angles, dihedrals
# -----------------------
if [ "$SKIP_ANALYSIS" != "true" ]; then
    
    echo ""
    echo "Step 4/4: Analyzing bonds, angles, and dihedrals"
    
    # Check if index files exist
    if [ ! -f "$BONDS_NDX" ]; then
        echo "Warning: Bonds index file not found: $BONDS_NDX"
        echo "Skipping bonds analysis"
    fi
    
    if [ ! -f "$ANGLES_NDX" ]; then
        echo "Warning: Angles index file not found: $ANGLES_NDX"
        echo "Skipping angles analysis"
    fi
    
    if [ ! -f "$DIHEDRALS_NDX" ]; then
        echo "Warning: Dihedrals index file not found: $DIHEDRALS_NDX"
        echo "Skipping dihedrals analysis"
    fi
    
    # Determine which trajectory to use for analysis
    if [ "$ANALYZE_REMOVE_PBC" = "true" ]; then
        ANALYZE_XTC="GMX/mapped.xtc"
        ANALYZE_TPR="GMX/CG.tpr"
        ANALYZE_ARGS="--remove_pbc --group_1 \"$ANALYZE_GROUP_1\" --group_2 \"$ANALYZE_GROUP_2\""
        if [ "$KEEP_INTERMEDIATE" = "true" ]; then
            ANALYZE_ARGS="$ANALYZE_ARGS --keep_intermediate"
        fi
    else
        ANALYZE_XTC="GMX/mapped.xtc"
        ANALYZE_TPR="GMX/CG.tpr"
        ANALYZE_ARGS=""
    fi
    
    # Run analysis if TPR exists
    if [ -f "$ANALYZE_TPR" ] && [ -f "$BONDS_NDX" ] && [ -f "$ANGLES_NDX" ] && [ -f "$DIHEDRALS_NDX" ]; then
        ANALYZE_CMD="python $PYTHON_SCRIPTS_DIR/generate_bonds_angles_dihedrals.py"
        ANALYZE_CMD_ARGS="--bonds_ndx $BONDS_NDX --angles_ndx $ANGLES_NDX --dihedrals_ndx $DIHEDRALS_NDX --xtc_file $ANALYZE_XTC --tpr_file $ANALYZE_TPR"
        
        if [ -n "$ANALYZE_ARGS" ]; then
            ANALYZE_CMD_ARGS="$ANALYZE_CMD_ARGS $ANALYZE_ARGS"
        fi
        
        # Run analysis from ANALYSIS directory
        cd ANALYSIS
        log_verbose "Running: $ANALYZE_CMD $ANALYZE_CMD_ARGS"
        eval $ANALYZE_CMD $ANALYZE_CMD_ARGS
        check_error "Bonds/angles/dihedrals analysis"
        
        # Move results to ANALYSIS directory (they're already there)
        cd ..
        echo "  ✓ Analysis completed"
        echo "    • Results in: ANALYSIS/bonds/, ANALYSIS/angles/, ANALYSIS/dihedrals/"
    else
        echo "Warning: Cannot run analysis. Missing required files:"
        [ ! -f "$ANALYZE_TPR" ] && echo "  - $ANALYZE_TPR"
        [ ! -f "$BONDS_NDX" ] && echo "  - $BONDS_NDX"
        [ ! -f "$ANGLES_NDX" ] && echo "  - $ANGLES_NDX"
        [ ! -f "$DIHEDRALS_NDX" ] && echo "  - $DIHEDRALS_NDX"
    fi
    
else
    echo ""
    echo "Step 4/4: Skipping bonds/angles/dihedrals analysis (--skip_analysis)"
fi

# -----------------------
# Cleanup temporary files
# -----------------------
if [ "$KEEP_TEMP" != "true" ]; then
    log_verbose "Cleaning up temporary files..."
    rm -rf OUTPUT_CG 2>/dev/null
    find . -name "*.pyc" -type f -delete 2>/dev/null
    find . -name "__pycache__" -type d -delete 2>/dev/null
    find . -name "#*" -type f -delete 2>/dev/null
    find . -name "*.1#" -type f -delete 2>/dev/null
    find . -name "*.2#" -type f -delete 2>/dev/null
fi

# -----------------------
# Create summary file
# -----------------------
cat > SUMMARY.txt << EOF
==========================================
CG Mapping Pipeline Summary
==========================================

Input files:
  AA TPR:           $AA_TPR_ABS
  AA XTC:           $AA_XTC_ABS
  AA GRO:           $AA_GRO_ABS
  AA ITP:           $AA_ITP_ABS
  Beads JSON:       $BEADS_JSON_ABS
  Force application:$FORCE_APPLICATION
  Beads position:   $BEADS_POSITION
  Input MDP:        $INPUT_MDP_ABS
  FF directory:     $PATH_FF_ABS

Output files:
  Main directory:   $(pwd)
  
GMX files (results/GMX/):
  - mapped.xtc      Mapped CG trajectory
  - cg.gro          CG coordinates
  - cg.itp          CG topology
  - topol_cg.top    GROMACS topology
  - CG.tpr          GROMACS TPR file (if grompp ran)

Analysis files (results/ANALYSIS/):
  - bonds/          Bond distance analysis
  - angles/         Angle analysis
  - dihedrals/      Dihedral analysis

NDX files (results/NDX/):
  - cg.ndx          CG bead mapping
  - bonds.ndx       Bond definitions
  - angles.ndx      Angle definitions
  - dihedrals.ndx   Dihedral definitions

==========================================
EOF

# -----------------------
# Done
# -----------------------
echo ""
echo "=========================================="
echo "✓ Pipeline completed successfully!"
echo "=========================================="
echo "Output directory: $(pwd)"
echo ""
echo "Directory structure:"
echo "  ├── GMX/        - GROMACS files (trajectory, topology, TPR)"
echo "  ├── ANALYSIS/   - Bonds, angles, dihedrals analysis"
echo "  └── NDX/        - Index files for CG mapping"
echo ""
echo "Summary saved in: SUMMARY.txt"
echo "=========================================="

cd "$ORIGINAL_DIR"
