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

# -----------------------
# Clean previous output directory
# -----------------------
if [ -d "$OUTPUT_DIR" ]; then
    echo "Cleaning previous output directory: $OUTPUT_DIR"
    rm -rf "$OUTPUT_DIR"
fi

# Create clean output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/GMX"
mkdir -p "$OUTPUT_DIR/CG_MARTINI3"
mkdir -p "$OUTPUT_DIR/XVG"
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
echo "  CG-MARTINI3 files: $OUTPUT_DIR/CG_MARTINI3/"
echo "  XVG files:         $OUTPUT_DIR/XVG/"
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
    --output_dir CG_MARTINI3_TEMP \
    --beads_position $BEADS_POSITION \
    --aa_itp $AA_ITP_ABS \
    --force_application \"$FORCE_APPLICATION\""

if [ "$VERBOSE" = "true" ]; then
    CG_MARTINI_CMD="$CG_MARTINI_CMD --verbose"
fi

log_verbose "Running: $CG_MARTINI_CMD"
eval $CG_MARTINI_CMD
check_error "cg-martini3"

# Move cg-martini3 outputs to CG_MARTINI3 directory
if [ -d "CG_MARTINI3_TEMP" ]; then
    mv CG_MARTINI3_TEMP/* CG_MARTINI3/ 2>/dev/null
    rm -rf CG_MARTINI3_TEMP
fi

# Copy NDX files to NDX directory
if [ -d "CG_MARTINI3/NDX" ]; then
    cp CG_MARTINI3/NDX/*.ndx NDX/ 2>/dev/null
    cp CG_MARTINI3/NDX/*.map NDX/ 2>/dev/null
fi

# Set paths for generated files
CG_NDX="NDX/cg.ndx"
BONDS_NDX="NDX/bonds.ndx"
ANGLES_NDX="NDX/angles.ndx"
DIHEDRALS_NDX="NDX/dihedrals.ndx"
CG_GRO="CG_MARTINI3/GRO/cg.gro"
CG_ITP_SRC="CG_MARTINI3/ITP/final_cg.itp"

log_verbose "✓ cg-martini3 completed successfully"

# -----------------------
# Step 1: Mapping (using map_aa_to_cg.py)
# -----------------------
echo ""
echo "Step 1/4: Mapping AA → CG trajectory"

# Try to use the module first, fallback to script
if python -c "import bayesian_potentials.scripts.map_aa_to_cg" 2>/dev/null; then
    MAP_CMD="python -m bayesian_potentials.scripts.map_aa_to_cg"
    log_verbose "Using Python module for mapping"
else
    MAP_SCRIPT=""
    if [ -f "$PACKAGE_DIR/scripts/map_aa_to_cg.py" ]; then
        MAP_SCRIPT="$PACKAGE_DIR/scripts/map_aa_to_cg.py"
    elif [ -f "$SCRIPT_DIR/../scripts/map_aa_to_cg.py" ]; then
        MAP_SCRIPT="$SCRIPT_DIR/../scripts/map_aa_to_cg.py"
    else
        echo "Error: Could not find map_aa_to_cg.py script"
        exit 1
    fi
    MAP_CMD="python $MAP_SCRIPT"
    log_verbose "Using script: $MAP_SCRIPT"
fi

MAP_CMD_ARGS="--index_cg $CG_NDX --aa_tpr $AA_TPR_ABS --aa_xtc $AA_XTC_ABS --output_mapped GMX/mapped.xtc --output_cg_gro GMX/cg.gro"

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

echo "  ✓ Generated: GMX/mapped.xtc"
echo "  ✓ Generated: GMX/cg.gro"

# -----------------------
# Step 2: Copy CG ITP
# -----------------------
echo ""
echo "Step 2/4: Copying CG ITP"

# Copy ITP from CG_MARTINI3 to GMX
if [ -f "$CG_ITP_SRC" ]; then
    cp "$CG_ITP_SRC" GMX/cg.itp
    echo "  ✓ Copied: GMX/cg.itp"
else
    echo "Warning: ITP file not found at $CG_ITP_SRC"
fi

# -----------------------
# Step 3: Generate topology and run grompp
# -----------------------
if [ "$SKIP_GROMPP" != "true" ]; then
    
    echo ""
    echo "Step 3/4: Generating topology and running grompp"
    
    # Create a temporary directory for topology generation
    TOP_TEMP_DIR="topology_temp"
    mkdir -p "$TOP_TEMP_DIR"
    
    # Copy necessary files to temp directory
    cp GMX/cg.itp "$TOP_TEMP_DIR/"
    cp GMX/cg.gro "$TOP_TEMP_DIR/"
    
    # Copy force field files to temp directory
    if [ -d "$PATH_FF_ABS" ]; then
        cp "$PATH_FF_ABS"/*.itp "$TOP_TEMP_DIR/" 2>/dev/null
        echo "  ✓ Copied force field files to temp directory"
    fi
    
    cd "$TOP_TEMP_DIR"
    
    # Try to use the module for topology generation
    if python -c "import bayesian_potentials.scripts.generate_cg_top" 2>/dev/null; then
        TOP_CMD="python -m bayesian_potentials.scripts.generate_cg_top"
        log_verbose "Using Python module for topology generation"
    else
        TOP_SCRIPT=""
        if [ -f "$PACKAGE_DIR/scripts/generate_cg_top.py" ]; then
            TOP_SCRIPT="$PACKAGE_DIR/scripts/generate_cg_top.py"
        elif [ -f "$SCRIPT_DIR/../scripts/generate_cg_top.py" ]; then
            TOP_SCRIPT="$SCRIPT_DIR/../scripts/generate_cg_top.py"
        else
            echo "Error: Could not find generate_cg_top.py script"
            exit 1
        fi
        TOP_CMD="python $TOP_SCRIPT"
        log_verbose "Using script: $TOP_SCRIPT"
    fi
    
    TOP_CMD_ARGS="--path_ff ."
    TOP_CMD_ARGS="$TOP_CMD_ARGS --ff $FF"
    TOP_CMD_ARGS="$TOP_CMD_ARGS --ions $IONS"
    TOP_CMD_ARGS="$TOP_CMD_ARGS --solvent $SOLVENT"
    TOP_CMD_ARGS="$TOP_CMD_ARGS --itp_ligand cg.itp"
    TOP_CMD_ARGS="$TOP_CMD_ARGS --name_molecule $NAME_MOLECULE"
    TOP_CMD_ARGS="$TOP_CMD_ARGS --number_molecule $NUMBER_MOLECULE"
    TOP_CMD_ARGS="$TOP_CMD_ARGS --title_comments \"$TITLE_COMMENTS\""
    TOP_CMD_ARGS="$TOP_CMD_ARGS --title_system \"$TITLE_SYSTEM\""
    TOP_CMD_ARGS="$TOP_CMD_ARGS --output_topol $OUTPUT_TOPOL"
    
    log_verbose "Running: $TOP_CMD $TOP_CMD_ARGS"
    eval $TOP_CMD $TOP_CMD_ARGS
    check_error "Topology generation"
    
    # Fix the topology file to remove any path prefixes from include statements
    if [ -f "$OUTPUT_TOPOL" ]; then
        # Remove any path prefixes from include statements
        sed -i 's|#include ".*/|#include "|g' "$OUTPUT_TOPOL"
        echo "  ✓ Generated and fixed: $OUTPUT_TOPOL"
    fi
    
    # Set default files for grompp (now in current directory)
    [ -z "$C_FILE" ] && C_FILE="cg.gro"
    [ -z "$P_FILE" ] && P_FILE="$OUTPUT_TOPOL"
    [ -z "$O_FILE" ] && O_FILE="CG.tpr"
    
    # Run grompp in the temp directory
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
    
    # Copy TPR back to GMX directory
    cp "$O_FILE" ../GMX/
    cp "$OUTPUT_TOPOL" ../GMX/
    
    cd ..
    rm -rf "$TOP_TEMP_DIR"
    
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
        # Try to use the module for analysis
        if python -c "import bayesian_potentials.scripts.generate_bonds_angles_dihedrals" 2>/dev/null; then
            ANALYZE_CMD="python -m bayesian_potentials.scripts.generate_bonds_angles_dihedrals"
            log_verbose "Using Python module for analysis"
        else
            ANALYZE_SCRIPT=""
            if [ -f "$PACKAGE_DIR/scripts/generate_bonds_angles_dihedrals.py" ]; then
                ANALYZE_SCRIPT="$PACKAGE_DIR/scripts/generate_bonds_angles_dihedrals.py"
            elif [ -f "$SCRIPT_DIR/../scripts/generate_bonds_angles_dihedrals.py" ]; then
                ANALYZE_SCRIPT="$SCRIPT_DIR/../scripts/generate_bonds_angles_dihedrals.py"
            else
                echo "Warning: Could not find analysis script"
                ANALYZE_SCRIPT=""
            fi
            ANALYZE_CMD="python $ANALYZE_SCRIPT"
        fi
        
        if [ -n "$ANALYZE_CMD" ]; then
            ANALYZE_CMD_ARGS="--bonds_ndx $BONDS_NDX --angles_ndx $ANGLES_NDX --dihedrals_ndx $DIHEDRALS_NDX --xtc_file $ANALYZE_XTC --tpr_file $ANALYZE_TPR"
            
            if [ -n "$ANALYZE_ARGS" ]; then
                ANALYZE_CMD_ARGS="$ANALYZE_CMD_ARGS $ANALYZE_ARGS"
            fi
            
            # Create temporary directory for analysis to avoid clutter
            ANALYSIS_TEMP_DIR="analysis_temp"
            mkdir -p "$ANALYSIS_TEMP_DIR"
            cd "$ANALYSIS_TEMP_DIR"
            
            eval $ANALYZE_CMD $ANALYZE_CMD_ARGS
            
            # Move XVG files to XVG directory
            if [ -d "bonds" ]; then
                mv bonds/*.xvg ../XVG/ 2>/dev/null
                mv bonds/distr_*.xvg ../XVG/ 2>/dev/null
                cp bonds/*.txt ../XVG/ 2>/dev/null
            fi
            
            if [ -d "angles" ]; then
                mv angles/*.xvg ../XVG/ 2>/dev/null
                mv angles/distr_*.xvg ../XVG/ 2>/dev/null
                cp angles/*.txt ../XVG/ 2>/dev/null
            fi
            
            if [ -d "dihedrals" ]; then
                mv dihedrals/*.xvg ../XVG/ 2>/dev/null
                mv dihedrals/distr_*.xvg ../XVG/ 2>/dev/null
                cp dihedrals/*.txt ../XVG/ 2>/dev/null
            fi
            
            cd ..
            rm -rf "$ANALYSIS_TEMP_DIR"
            
            echo "  ✓ Analysis completed"
            echo "    • XVG files in: XVG/"
        fi
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
  
GMX files (GMX/):
  - mapped.xtc      Mapped CG trajectory
  - cg.gro          CG coordinates
  - cg.itp          CG topology
  - topol_cg.top    GROMACS topology
  - CG.tpr          GROMACS TPR file

CG-MARTINI3 files (CG_MARTINI3/):
  - GRO/            CG structure files
  - ITP/            CG topology files
  - JSON/           Beads definition files
  - NDX/            Index files

XVG files (XVG/):
  - bond_*.xvg      Bond distance time series
  - distr_bond_*.xvg Bond distance distributions
  - ang_*.xvg       Angle time series
  - distr_ang_*.xvg Angle distributions
  - dih_*.xvg       Dihedral time series
  - distr_dih_*.xvg Dihedral distributions

NDX files (NDX/):
  - cg.ndx          CG bead mapping
  - bonds.ndx       Bond definitions
  - angles.ndx      Angle definitions
  - dihedrals.ndx   Dihedral definitions
  - cg.map          Atom to bead mapping

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
echo "  ├── GMX/           - GROMACS files (trajectory, topology, TPR)"
echo "  ├── CG_MARTINI3/   - cg-martini3 outputs (GRO, ITP, JSON, NDX)"
echo "  ├── XVG/           - Analysis XVG files (bonds, angles, dihedrals)"
echo "  └── NDX/           - Index files for CG mapping"
echo ""
echo "Summary saved in: SUMMARY.txt"
echo "=========================================="

cd "$ORIGINAL_DIR"
