#!/bin/bash

# =====================================================
# Script: mapping_aa_to_cg.sh
# Description:
# - Runs cg-martini3 to generate CG mapping
# - Maps atomistic trajectory to coarse-grained (CG)
# - Generates topology, runs GROMACS preprocessing
# - Analyzes bonds/angles/dihedrals from CG trajectories
# - Calculates distribution statistics (optional)
# - PREPARES CG SYSTEM FOR SIMULATION (bp_prep.py)
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
    echo "  --force_application STR   Force application (default: 'random=[1250,30;30,4]')"
    echo "  --beads_position STR      Beads position (default: 'com')"
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
    echo "  --cycle_restr STR         Cycle constraints: 'none' (default), 'fix=3,mode=cycle', or 'mode=linear'"
    echo "  --default_martini         Use default Martini 3 masses (72) and zero charges"
    echo "  --maxwarn N               Max warnings for grompp (default: 1)"
    echo "  --remove_pbc              Remove PBC (default: true)"
    echo "  --no_pbc                  Skip PBC removal"
    echo "  --skip_adapt_itp          Skip ITP adaptation step"
    echo "  --skip_grompp             Skip grompp step"
    echo "  --skip_analysis           Skip bonds/angles/dihedrals analysis"
    echo "  --skip_distributions      Skip distribution statistics calculation"
    echo "  --skip_cg_prep            Skip CG system preparation (bp_prep.py)"
    echo "  --analyze_remove_pbc      Remove PBC before analysis"
    echo "  --keep_temp               Keep temporary files"
    echo "  --verbose                 Verbose output"
    echo ""
    echo "CG SYSTEM PREPARATION ARGUMENTS (bp_prep.py):"
    echo "  --cg_prep_use_distance     Use distance from atom for box size"
    echo "  --cg_prep_distance NM      Distance from atom in nm (default: 1.2)"
    echo "  --cg_prep_box_size NM      Fixed box size in nm"
    echo "  --cg_prep_max_solvent N    Maximum number of solvent molecules"
    echo "  --cg_prep_salt CONC        Salt concentration in M (default: 0.15)"
    echo "  --cg_prep_water_gro FILE   Water GRO filename (default: water.gro)"
    echo "  --cg_prep_water_dir DIR    Water GRO directory"
    echo "  --cg_prep_ions_mdp DIR     Ions MDP directory"
    echo "  --cg_prep_ions_file FILE   Ions MDP filename (default: ions.mdp)"
    echo "  --cg_prep_min_mdp FILE     Minimization MDP filename (default: minimization.mdp)"
    echo ""
    echo "DISTRIBUTION STATISTICS ARGUMENTS (optional):"
    echo "  --run_distributions       Run distribution statistics calculation"
    echo "  --dist_bonds_dir DIR      Directory for bond XVG files (default: bonds)"
    echo "  --dist_angles_dir DIR     Directory for angle XVG files (default: angles)"
    echo "  --dist_dihedrals_dir DIR  Directory for dihedral XVG files (default: dihedrals)"
    echo "  --dist_output_dir DIR     Output directory for statistics (default: STATISTICS)"
    echo "  --dist_bond_out FILE      Bond statistics output (default: bond_statistics.tsv)"
    echo "  --dist_angle_out FILE     Angle statistics output (default: angle_statistics.tsv)"
    echo "  --dist_dihedral_out FILE  Dihedral statistics output (default: dihedral_statistics.tsv)"
    echo ""
    echo "Examples:"
    echo "  # Default run (with random potentials and cycle constraints)"
    echo "  $0 --aa_tpr setup/md.tpr --aa_xtc setup/md.xtc --aa_gro setup/md.gro \\"
    echo "     --aa_itp setup/carb.itp --beads_json json/beads_config.json \\"
    echo "     --input_mdp mdp/minimization.mdp --path_ff ff_files/"
    echo ""
    echo "  # Run with CG system preparation"
    echo "  $0 --aa_tpr setup/md.tpr --aa_xtc setup/md.xtc --aa_gro setup/md.gro \\"
    echo "     --aa_itp setup/carb.itp --beads_json json/beads_config.json \\"
    echo "     --input_mdp mdp/minimization.mdp --path_ff ff_files/ \\"
    echo "     --cg_prep_use_distance --cg_prep_distance 1.2 --cg_prep_salt 0.15"
    echo ""
    echo "  # Custom water and ions configuration"
    echo "  $0 --aa_tpr setup/md.tpr --aa_xtc setup/md.xtc --aa_gro setup/md.gro \\"
    echo "     --aa_itp setup/carb.itp --beads_json json/beads_config.json \\"
    echo "     --input_mdp mdp/minimization.mdp --path_ff ff_files/ \\"
    echo "     --cg_prep_use_distance --cg_prep_water_dir /path/to/water \\"
    echo "     --cg_prep_water_file my_water.gro --cg_prep_ions_mdp /path/to/mdp"
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
# Default values for cg-martini3
FORCE_APPLICATION="random=[1250,30;30,4]"  # Default: random potentials
BEADS_POSITION="com"                       # Default: center of mass
CYCLE_RESTR="fix=3,mode=cycle"            # Default: 3-member cycles
DEFAULT_MARTINI="false"
MAXWARN=1
REMOVE_PBC="true"
SKIP_ADAPT_ITP="false"
SKIP_GROMPP="false"
SKIP_ANALYSIS="false"
SKIP_DISTRIBUTIONS="false"
SKIP_CG_PREP="false"                       # NEW: Skip CG preparation
RUN_DISTRIBUTIONS="false"
ANALYZE_REMOVE_PBC="false"
ANALYZE_GROUP_1="System"
ANALYZE_GROUP_2="System"
KEEP_INTERMEDIATE="false"
KEEP_TEMP="false"
VERBOSE="false"

# CG Preparation defaults (for bp_prep.py)
CG_PREP_USE_DISTANCE="false"
CG_PREP_DISTANCE="1.2"
CG_PREP_BOX_SIZE=""
CG_PREP_MAX_SOLVENT=""
CG_PREP_SALT="0.15"
CG_PREP_WATER_GRO="water.gro"
CG_PREP_WATER_DIR=""
CG_PREP_IONS_MDP_DIR=""
CG_PREP_IONS_MDP_FILE="ions.mdp"
CG_PREP_MIN_MDP_FILE="minimization.mdp"

# Distribution statistics defaults
DIST_BONDS_DIR="bonds"
DIST_ANGLES_DIR="angles"
DIST_DIHEDRALS_DIR="dihedrals"
DIST_OUTPUT_DIR="STATISTICS"
DIST_BOND_OUT="bond_statistics.tsv"
DIST_ANGLE_OUT="angle_statistics.tsv"
DIST_DIHEDRAL_OUT="dihedral_statistics.tsv"

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
        --cycle_restr) CYCLE_RESTR="$2"; shift 2 ;;
        --default_martini) DEFAULT_MARTINI="true"; shift ;;
        --maxwarn) MAXWARN="$2"; shift 2 ;;
        --remove_pbc) REMOVE_PBC="true"; shift ;;
        --no_pbc) REMOVE_PBC="false"; shift ;;
        --skip_adapt_itp) SKIP_ADAPT_ITP="true"; shift ;;
        --skip_grompp) SKIP_GROMPP="true"; shift ;;
        --skip_analysis) SKIP_ANALYSIS="true"; shift ;;
        --skip_distributions) SKIP_DISTRIBUTIONS="true"; shift ;;
        --skip_cg_prep) SKIP_CG_PREP="true"; shift ;;
        --run_distributions) RUN_DISTRIBUTIONS="true"; shift ;;
        --analyze_remove_pbc) ANALYZE_REMOVE_PBC="true"; shift ;;
        --analyze_group_1) ANALYZE_GROUP_1="$2"; shift 2 ;;
        --analyze_group_2) ANALYZE_GROUP_2="$2"; shift 2 ;;
        --keep_intermediate) KEEP_INTERMEDIATE="true"; shift ;;
        --keep_temp) KEEP_TEMP="true"; shift ;;
        --verbose) VERBOSE="true"; shift ;;
        
        # CG Preparation arguments
        --cg_prep_use_distance) CG_PREP_USE_DISTANCE="true"; shift ;;
        --cg_prep_distance) CG_PREP_DISTANCE="$2"; shift 2 ;;
        --cg_prep_box_size) CG_PREP_BOX_SIZE="$2"; shift 2 ;;
        --cg_prep_max_solvent) CG_PREP_MAX_SOLVENT="$2"; shift 2 ;;
        --cg_prep_salt) CG_PREP_SALT="$2"; shift 2 ;;
        --cg_prep_water_gro) CG_PREP_WATER_GRO="$2"; shift 2 ;;
        --cg_prep_water_dir) CG_PREP_WATER_DIR="$2"; shift 2 ;;
        --cg_prep_ions_mdp) CG_PREP_IONS_MDP_DIR="$2"; shift 2 ;;
        --cg_prep_ions_file) CG_PREP_IONS_MDP_FILE="$2"; shift 2 ;;
        --cg_prep_min_mdp) CG_PREP_MIN_MDP_FILE="$2"; shift 2 ;;
        
        # Distribution statistics arguments
        --dist_bonds_dir) DIST_BONDS_DIR="$2"; shift 2 ;;
        --dist_angles_dir) DIST_ANGLES_DIR="$2"; shift 2 ;;
        --dist_dihedrals_dir) DIST_DIHEDRALS_DIR="$2"; shift 2 ;;
        --dist_output_dir) DIST_OUTPUT_DIR="$2"; shift 2 ;;
        --dist_bond_out) DIST_BOND_OUT="$2"; shift 2 ;;
        --dist_angle_out) DIST_ANGLE_OUT="$2"; shift 2 ;;
        --dist_dihedral_out) DIST_DIHEDRAL_OUT="$2"; shift 2 ;;
        
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

# -----------------------
# Validate essential arguments
# -----------------------
if [ -z "$AA_TPR" ] || [ -z "$AA_XTC" ] || [ -z "$AA_GRO" ] || \
   [ -z "$AA_ITP" ] || [ -z "$BEADS_JSON" ] || \
   [ -z "$INPUT_MDP" ] || [ -z "$PATH_FF" ]; then
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
mkdir -p "$OUTPUT_DIR/STATISTICS"
mkdir -p "$OUTPUT_DIR/MDRUN_CG"  # NEW: Directory for bp_prep.py output

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
echo "  Cycle restraints:  $CYCLE_RESTR"
echo "  Input MDP:         $INPUT_MDP_ABS"
echo "  FF directory:      $PATH_FF_ABS"
echo ""
echo "Output directories:"
echo "  Main:              $OUTPUT_DIR"
echo "  GMX files:         $OUTPUT_DIR/GMX/"
echo "  CG-MARTINI3 files: $OUTPUT_DIR/CG_MARTINI3/"
echo "  XVG files:         $OUTPUT_DIR/XVG/"
echo "  NDX files:         $OUTPUT_DIR/NDX/"
echo "  Statistics:        $OUTPUT_DIR/STATISTICS/"
echo "  CG Run:            $OUTPUT_DIR/MDRUN_CG/"
echo ""
echo "Optional arguments:"
echo "  Molecule name:     $NAME_MOLECULE"
echo "  Number molecules:  $NUMBER_MOLECULE"
echo "  Remove PBC:        $REMOVE_PBC"
echo "  Skip adapt ITP:    $SKIP_ADAPT_ITP"
echo "  Skip grompp:       $SKIP_GROMPP"
echo "  Skip analysis:     $SKIP_ANALYSIS"
echo "  Skip CG prep:      $SKIP_CG_PREP"
echo "  Run distributions: $RUN_DISTRIBUTIONS"
echo "  Analyze remove PBC:$ANALYZE_REMOVE_PBC"
echo "  Verbose:           $VERBOSE"
echo ""
if [ "$RUN_DISTRIBUTIONS" = "true" ]; then
    echo "Distribution statistics:"
    echo "  Bonds dir:         $DIST_BONDS_DIR"
    echo "  Angles dir:        $DIST_ANGLES_DIR"
    echo "  Dihedrals dir:     $DIST_DIHEDRALS_DIR"
    echo "  Output dir:        $DIST_OUTPUT_DIR"
    echo "  Bond output:       $DIST_BOND_OUT"
    echo "  Angle output:      $DIST_ANGLE_OUT"
    echo "  Dihedral output:   $DIST_DIHEDRAL_OUT"
fi
if [ "$SKIP_CG_PREP" != "true" ]; then
    echo ""
    echo "CG System Preparation (bp_prep.py):"
    if [ "$CG_PREP_USE_DISTANCE" = "true" ]; then
        echo "  Box method:        distance from atom ($CG_PREP_DISTANCE nm)"
    elif [ -n "$CG_PREP_BOX_SIZE" ]; then
        echo "  Box method:        fixed size ($CG_PREP_BOX_SIZE nm)"
    else
        echo "  Box method:        default padding (1.5 nm)"
    fi
    echo "  Salt concentration: $CG_PREP_SALT M"
    echo "  Water file:         $CG_PREP_WATER_GRO"
    [ -n "$CG_PREP_MAX_SOLVENT" ] && echo "  Max solvent:        $CG_PREP_MAX_SOLVENT"
fi
echo "=========================================="

# -----------------------
# Step 0: Run cg-martini3 to generate CG mapping
# -----------------------
echo ""
echo "Step 0/7: Running cg-martini3 to generate CG mapping"

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
    --force_application \"$FORCE_APPLICATION\" \
    --cycle_restr \"$CYCLE_RESTR\" \
    --name_molecule $NAME_MOLECULE"

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
echo "Step 1/7: Mapping AA → CG trajectory"

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
echo "Step 2/7: Copying CG ITP"

# Copy ITP from CG_MARTINI3 to GMX
if [ -f "$CG_ITP_SRC" ]; then
    cp "$CG_ITP_SRC" GMX/cg.itp
    echo "  ✓ Copied: GMX/cg.itp (molecule name: $NAME_MOLECULE)"
else
    echo "Warning: ITP file not found at $CG_ITP_SRC"
fi

# -----------------------
# Step 2.5: Adapt ITP to match GRO atom names
# -----------------------
if [ "$SKIP_ADAPT_ITP" != "true" ] && [ -f "GMX/cg.itp" ] && [ -f "GMX/cg.gro" ]; then
    
    echo ""
    echo "Step 2.5/7: Adapting ITP atom names to match GRO reference"
    
    # Create a backup of original ITP
    cp GMX/cg.itp GMX/cg.itp.original
    
    # Try to use the module for ITP adaptation
    if python -c "import bayesian_potentials.scripts.adaptation_gro_itp" 2>/dev/null; then
        ADAPT_CMD="python -m bayesian_potentials.scripts.adaptation_gro_itp"
        log_verbose "Using Python module for ITP adaptation"
    else
        ADAPT_SCRIPT=""
        if [ -f "$PACKAGE_DIR/scripts/adaptation_gro_itp.py" ]; then
            ADAPT_SCRIPT="$PACKAGE_DIR/scripts/adaptation_gro_itp.py"
        elif [ -f "$SCRIPT_DIR/../scripts/adaptation_gro_itp.py" ]; then
            ADAPT_SCRIPT="$SCRIPT_DIR/../scripts/adaptation_gro_itp.py"
        else
            echo "Warning: Could not find adaptation_gro_itp.py script"
            ADAPT_SCRIPT=""
        fi
        ADAPT_CMD="python $ADAPT_SCRIPT"
    fi
    
    if [ -n "$ADAPT_CMD" ]; then
        ADAPT_CMD_ARGS="--input_itp GMX/cg.itp"
        ADAPT_CMD_ARGS="$ADAPT_CMD_ARGS --input_gro_ref GMX/cg.gro"
        ADAPT_CMD_ARGS="$ADAPT_CMD_ARGS --output_itp_adapted GMX/cg.itp.adapted"
        
        if [ "$VERBOSE" = "true" ]; then
            ADAPT_CMD_ARGS="$ADAPT_CMD_ARGS --verbose"
        fi
        
        log_verbose "Running: $ADAPT_CMD $ADAPT_CMD_ARGS"
        eval $ADAPT_CMD $ADAPT_CMD_ARGS
        
        if [ $? -eq 0 ]; then
            # Replace original ITP with adapted version
            mv GMX/cg.itp.adapted GMX/cg.itp
            echo "  ✓ ITP successfully adapted"
            if [ "$VERBOSE" = "true" ]; then
                echo "    Original ITP backed up to: GMX/cg.itp.original"
            fi
        else
            echo "  ✗ ITP adaptation failed, using original ITP"
            rm -f GMX/cg.itp.adapted
        fi
    else
        echo "  ⚠ ITP adaptation script not found, skipping"
    fi
    
else
    if [ "$SKIP_ADAPT_ITP" = "true" ]; then
        echo ""
        echo "Step 2.5/7: Skipping ITP adaptation (--skip_adapt_itp)"
    else
        echo ""
        echo "Step 2.5/7: Cannot adapt ITP - missing files"
        [ ! -f "GMX/cg.itp" ] && echo "  ✗ Missing: GMX/cg.itp"
        [ ! -f "GMX/cg.gro" ] && echo "  ✗ Missing: GMX/cg.gro"
    fi
fi

# -----------------------
# Step 3: Generate topology and run grompp
# -----------------------
if [ "$SKIP_GROMPP" != "true" ]; then
    
    echo ""
    echo "Step 3/7: Generating topology and running grompp"
    
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
    
    # Copy TPR and topology back to GMX directory
    cp "$O_FILE" ../GMX/
    cp "$OUTPUT_TOPOL" ../GMX/
    
    cd ..
    rm -rf "$TOP_TEMP_DIR"
    
else
    echo ""
    echo "Step 3/7: Skipping grompp (--skip_grompp)"
fi

# -----------------------
# Step 4: Analyze bonds, angles, dihedrals
# -----------------------
if [ "$SKIP_ANALYSIS" != "true" ]; then
    
    echo ""
    echo "Step 4/7: Analyzing bonds, angles, and dihedrals"
    
    # Ensure NDX files are in the correct location
    if [ -d "CG_MARTINI3/NDX" ]; then
        echo "  Copying NDX files from CG_MARTINI3/NDX to NDX/..."
        mkdir -p NDX
        cp CG_MARTINI3/NDX/*.ndx NDX/ 2>/dev/null
        cp CG_MARTINI3/NDX/*.map NDX/ 2>/dev/null
    fi
    
    # Check if index files exist
    BONDS_NDX="NDX/bonds.ndx"
    ANGLES_NDX="NDX/angles.ndx"
    DIHEDRALS_NDX="NDX/dihedrals.ndx"
    
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
    
    # Run analysis if all required files exist
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
            # Change to results directory to run analysis
            cd "$ORIGINAL_DIR/$OUTPUT_DIR"
            
            ANALYZE_CMD_ARGS="--bonds_ndx $BONDS_NDX --angles_ndx $ANGLES_NDX --dihedrals_ndx $DIHEDRALS_NDX --xtc_file $ANALYZE_XTC --tpr_file $ANALYZE_TPR"
            
            if [ -n "$ANALYZE_ARGS" ]; then
                ANALYZE_CMD_ARGS="$ANALYZE_CMD_ARGS $ANALYZE_ARGS"
            fi
            
            echo "  Running: $ANALYZE_CMD $ANALYZE_CMD_ARGS"
            eval $ANALYZE_CMD $ANALYZE_CMD_ARGS
            
            if [ $? -eq 0 ]; then
                echo "  ✓ Analysis completed successfully"
                
                # Create subdirectories in XVG/
                echo "  Organizing XVG files into subdirectories..."
                mkdir -p XVG/bonds XVG/angles XVG/dihedrals
                
                # Move bond XVG files
                if [ -d "bonds" ]; then
                    echo "    Moving bond files to XVG/bonds/"
                    for file in bonds/*.xvg; do
                        if [ -f "$file" ]; then
                            mv "$file" XVG/bonds/ 2>/dev/null
                        fi
                    done
                    for file in bonds/distr_*.xvg; do
                        if [ -f "$file" ]; then
                            mv "$file" XVG/bonds/ 2>/dev/null
                        fi
                    done
                    # Copy text files
                    cp bonds/*.txt XVG/bonds/ 2>/dev/null
                    rm -rf bonds
                fi
                
                # Move angle XVG files
                if [ -d "angles" ]; then
                    echo "    Moving angle files to XVG/angles/"
                    for file in angles/*.xvg; do
                        if [ -f "$file" ]; then
                            mv "$file" XVG/angles/ 2>/dev/null
                        fi
                    done
                    for file in angles/distr_*.xvg; do
                        if [ -f "$file" ]; then
                            mv "$file" XVG/angles/ 2>/dev/null
                        fi
                    done
                    # Copy text files
                    cp angles/*.txt XVG/angles/ 2>/dev/null
                    rm -rf angles
                fi
                
                # Move dihedral XVG files
                if [ -d "dihedrals" ]; then
                    echo "    Moving dihedral files to XVG/dihedrals/"
                    for file in dihedrals/*.xvg; do
                        if [ -f "$file" ]; then
                            mv "$file" XVG/dihedrals/ 2>/dev/null
                        fi
                    done
                    for file in dihedrals/distr_*.xvg; do
                        if [ -f "$file" ]; then
                            mv "$file" XVG/dihedrals/ 2>/dev/null
                        fi
                    done
                    # Copy text files
                    cp dihedrals/*.txt XVG/dihedrals/ 2>/dev/null
                    rm -rf dihedrals
                fi
                
                # Count files in each subdirectory
                BOND_COUNT=$(ls -1 XVG/bonds/*.xvg 2>/dev/null | wc -l)
                ANGLE_COUNT=$(ls -1 XVG/angles/*.xvg 2>/dev/null | wc -l)
                DIH_COUNT=$(ls -1 XVG/dihedrals/*.xvg 2>/dev/null | wc -l)
                
                echo "    ✓ Bond files: $BOND_COUNT"
                echo "    ✓ Angle files: $ANGLE_COUNT"
                echo "    ✓ Dihedral files: $DIH_COUNT"
            else
                echo "  ✗ Analysis failed with error code $?"
            fi
            
            cd "$ORIGINAL_DIR"
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
    echo "Step 4/7: Skipping bonds/angles/dihedrals analysis (--skip_analysis)"
fi

# -----------------------
# Step 5: Calculate distribution statistics (optional)
# -----------------------
if [ "$RUN_DISTRIBUTIONS" = "true" ] && [ "$SKIP_DISTRIBUTIONS" != "true" ] && [ "$SKIP_ANALYSIS" != "true" ]; then
    
    echo ""
    echo "Step 5/7: Calculating distribution statistics"
    
    # Change to results directory
    cd "$ORIGINAL_DIR/$OUTPUT_DIR"
    
    # Check if XVG files exist
    if [ -d "XVG" ]; then
        
        # Create symbolic links or copy XVG files to expected directories
        # The bp_distributions.py script expects bonds/, angles/, dihedrals/ directories
        # with XVG files inside
        
        echo "  Preparing XVG files for statistics calculation..."
        
        # Create directories if they don't exist
        mkdir -p "$DIST_BONDS_DIR"
        mkdir -p "$DIST_ANGLES_DIR"
        mkdir -p "$DIST_DIHEDRALS_DIR"
        
        # Copy/link bond XVG files
        if [ -d "XVG/bonds" ]; then
            echo "    Copying bond files to $DIST_BONDS_DIR/"
            cp XVG/bonds/*.xvg "$DIST_BONDS_DIR/" 2>/dev/null
        fi
        
        # Copy/link angle XVG files
        if [ -d "XVG/angles" ]; then
            echo "    Copying angle files to $DIST_ANGLES_DIR/"
            cp XVG/angles/*.xvg "$DIST_ANGLES_DIR/" 2>/dev/null
        fi
        
        # Copy/link dihedral XVG files
        if [ -d "XVG/dihedrals" ]; then
            echo "    Copying dihedral files to $DIST_DIHEDRALS_DIR/"
            cp XVG/dihedrals/*.xvg "$DIST_DIHEDRALS_DIR/" 2>/dev/null
        fi
        
        # Run distribution statistics script
        if python -c "import bayesian_potentials.scripts.bp_distributions" 2>/dev/null; then
            DIST_CMD="python -m bayesian_potentials.scripts.bp_distributions"
            log_verbose "Using Python module for distribution statistics"
        else
            DIST_SCRIPT=""
            if [ -f "$PACKAGE_DIR/scripts/bp_distributions.py" ]; then
                DIST_SCRIPT="$PACKAGE_DIR/scripts/bp_distributions.py"
            elif [ -f "$SCRIPT_DIR/../scripts/bp_distributions.py" ]; then
                DIST_SCRIPT="$SCRIPT_DIR/../scripts/bp_distributions.py"
            else
                echo "Warning: Could not find bp_distributions.py script"
                DIST_SCRIPT=""
            fi
            DIST_CMD="python $DIST_SCRIPT"
        fi
        
        if [ -n "$DIST_CMD" ]; then
            # Create statistics output directory
            mkdir -p "$DIST_OUTPUT_DIR"
            
            DIST_CMD_ARGS="--bonds_dir $DIST_BONDS_DIR"
            DIST_CMD_ARGS="$DIST_CMD_ARGS --angles_dir $DIST_ANGLES_DIR"
            DIST_CMD_ARGS="$DIST_CMD_ARGS --dihedrals_dir $DIST_DIHEDRALS_DIR"
            DIST_CMD_ARGS="$DIST_CMD_ARGS --dir_to_output $DIST_OUTPUT_DIR"
            DIST_CMD_ARGS="$DIST_CMD_ARGS --bond_out $DIST_BOND_OUT"
            DIST_CMD_ARGS="$DIST_CMD_ARGS --angle_out $DIST_ANGLE_OUT"
            DIST_CMD_ARGS="$DIST_CMD_ARGS --dihedral_out $DIST_DIHEDRAL_OUT"
            
            echo "  Running: $DIST_CMD $DIST_CMD_ARGS"
            eval $DIST_CMD $DIST_CMD_ARGS
            
            if [ $? -eq 0 ]; then
                echo "  ✓ Distribution statistics completed successfully"
                
                # Count statistics files
                STAT_COUNT=$(ls -1 "$DIST_OUTPUT_DIR"/*.tsv 2>/dev/null | wc -l)
                echo "    Generated $STAT_COUNT statistics files in $DIST_OUTPUT_DIR/"
                
                # List generated files
                if [ "$VERBOSE" = "true" ]; then
                    echo "    Statistics files:"
                    ls -la "$DIST_OUTPUT_DIR"/*.tsv 2>/dev/null | awk '{print "      - " $9}'
                fi
            else
                echo "  ✗ Distribution statistics failed with error code $?"
            fi
        fi
        
        # Clean up temporary directories if not keeping temp
        if [ "$KEEP_TEMP" != "true" ]; then
            echo "  Cleaning up temporary directories..."
            rm -rf "$DIST_BONDS_DIR" "$DIST_ANGLES_DIR" "$DIST_DIHEDRALS_DIR"
        fi
        
    else
        echo "Warning: XVG directory not found. Cannot calculate distribution statistics."
    fi
    
    cd "$ORIGINAL_DIR"
    
else
    if [ "$RUN_DISTRIBUTIONS" = "true" ] && [ "$SKIP_ANALYSIS" = "true" ]; then
        echo ""
        echo "Step 5/7: Skipping distribution statistics (analysis was skipped)"
    elif [ "$RUN_DISTRIBUTIONS" != "true" ]; then
        echo ""
        echo "Step 5/7: Skipping distribution statistics (not requested)"
    fi
fi

# -----------------------
# Step 6: Prepare CG system for simulation (bp_prep.py)
# -----------------------
if [ "$SKIP_CG_PREP" != "true" ]; then
    
    echo ""
    echo "Step 6/7: Preparing CG system for simulation (bp_prep.py)"
    
    # Set paths for bp_prep.py input files
    MDRUN_CG_DIR="$OUTPUT_DIR/MDRUN_CG"
    
    # Determine input_ref_dir (where GMX files are)
    INPUT_REF_DIR="$ORIGINAL_DIR/$OUTPUT_DIR/GMX"
    
    # Check if required files exist in GMX directory
    if [ ! -f "$INPUT_REF_DIR/cg.gro" ]; then
        echo "  ✗ Error: cg.gro not found in $INPUT_REF_DIR"
        echo "  Skipping CG system preparation"
    elif [ ! -f "$INPUT_REF_DIR/cg.itp" ]; then
        echo "  ✗ Error: cg.itp not found in $INPUT_REF_DIR"
        echo "  Skipping CG system preparation"
    elif [ ! -f "$INPUT_REF_DIR/topol_cg.top" ] && [ ! -f "$INPUT_REF_DIR/$OUTPUT_TOPOL" ]; then
        echo "  ✗ Error: Topology file not found in $INPUT_REF_DIR"
        echo "  Skipping CG system preparation"
    else
        # Find topology file
        if [ -f "$INPUT_REF_DIR/topol_cg.top" ]; then
            TOPOL_FILE="topol_cg.top"
        else
            TOPOL_FILE="$OUTPUT_TOPOL"
        fi
        
        # Build bp_prep.py command
        BP_PREP_CMD="python $SCRIPT_DIR/../bp_prep.py"
        
        # Essential arguments
        BP_PREP_CMD="$BP_PREP_CMD --input_ref_dir $INPUT_REF_DIR"
        BP_PREP_CMD="$BP_PREP_CMD --input_gro cg.gro"
        BP_PREP_CMD="$BP_PREP_CMD --input_itp cg.itp"
        BP_PREP_CMD="$BP_PREP_CMD --input_topol $TOPOL_FILE"
        BP_PREP_CMD="$BP_PREP_CMD --output_dir $MDRUN_CG_DIR"
        
        # Force field files
        BP_PREP_CMD="$BP_PREP_CMD --ff $FF"
        BP_PREP_CMD="$BP_PREP_CMD --ions $IONS"
        BP_PREP_CMD="$BP_PREP_CMD --solvent $SOLVENT"
        
        # Path to ff_files (use the provided path_ff)
        if [ -n "$PATH_FF_ABS" ]; then
            BP_PREP_CMD="$BP_PREP_CMD --input_ff_dir $PATH_FF_ABS"
        fi
        
        # Water file
        if [ -n "$CG_PREP_WATER_DIR" ]; then
            BP_PREP_CMD="$BP_PREP_CMD --water_dir $CG_PREP_WATER_DIR"
        fi
        BP_PREP_CMD="$BP_PREP_CMD --water_file_gro $CG_PREP_WATER_GRO"
        
        # Ions MDP
        if [ -n "$CG_PREP_IONS_MDP_DIR" ]; then
            BP_PREP_CMD="$BP_PREP_CMD --ions_mdp_dir $CG_PREP_IONS_MDP_DIR"
        fi
        BP_PREP_CMD="$BP_PREP_CMD --ions_file_mdp $CG_PREP_IONS_MDP_FILE"
        
        # Minimization MDP
        BP_PREP_CMD="$BP_PREP_CMD --input_name_file_mdp $CG_PREP_MIN_MDP_FILE"
        
        # Box definition
        if [ "$CG_PREP_USE_DISTANCE" = "true" ]; then
            BP_PREP_CMD="$BP_PREP_CMD --use_distance_from_atom"
            BP_PREP_CMD="$BP_PREP_CMD --distance_from_atom $CG_PREP_DISTANCE"
        elif [ -n "$CG_PREP_BOX_SIZE" ]; then
            BP_PREP_CMD="$BP_PREP_CMD --box_size $CG_PREP_BOX_SIZE"
        fi
        
        # Solvation parameters
        if [ -n "$CG_PREP_MAX_SOLVENT" ]; then
            BP_PREP_CMD="$BP_PREP_CMD --max_solvent $CG_PREP_MAX_SOLVENT"
        fi
        
        # Salt concentration
        BP_PREP_CMD="$BP_PREP_CMD --salt $CG_PREP_SALT"
        
        # Verbose mode
        if [ "$VERBOSE" = "true" ]; then
            BP_PREP_CMD="$BP_PREP_CMD --verbose"
        fi
        
        echo "  Running: $BP_PREP_CMD"
        eval $BP_PREP_CMD
        
        if [ $? -eq 0 ]; then
            echo "  ✓ CG system preparation completed successfully"
            echo "    Output in: $MDRUN_CG_DIR/"
        else
            echo "  ✗ CG system preparation failed"
        fi
    fi
    
else
    echo ""
    echo "Step 6/7: Skipping CG system preparation (--skip_cg_prep)"
fi

# -----------------------
# Step 7: Cleanup temporary files
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
cat > "$OUTPUT_DIR/SUMMARY.txt" << EOF
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
  Cycle restraints: $CYCLE_RESTR
  Input MDP:        $INPUT_MDP_ABS
  FF directory:     $PATH_FF_ABS

Output files:
  Main directory:   $(pwd)
  
GMX files (GMX/):
  - mapped.xtc      Mapped CG trajectory
  - cg.gro          CG coordinates
  - cg.itp          CG topology (adapted)
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

EOF

# Add statistics section if distributions were run
if [ "$RUN_DISTRIBUTIONS" = "true" ] && [ "$SKIP_ANALYSIS" != "true" ]; then
    cat >> "$OUTPUT_DIR/SUMMARY.txt" << EOF

Statistics files (STATISTICS/):
  - $DIST_BOND_OUT     Bond distribution statistics
  - $DIST_ANGLE_OUT    Angle distribution statistics
  - $DIST_DIHEDRAL_OUT Dihedral distribution statistics

EOF
fi

# Add CG preparation section
if [ "$SKIP_CG_PREP" != "true" ]; then
    cat >> "$OUTPUT_DIR/SUMMARY.txt" << EOF

CG System for Simulation (MDRUN_CG/):
  - solv_ions_CG.gro  Solvated and neutralized system
  - topol_cg.top      Updated topology with ions
  - em.tpr            Energy minimization input
  - ff_files/         Martini force field parameters
  - mdp/              MDP files for simulation

EOF
fi

echo "✓ Summary saved in: $OUTPUT_DIR/SUMMARY.txt"
cd "$ORIGINAL_DIR"

# -----------------------
# Final message
# -----------------------
echo ""
echo "=========================================="
echo "Pipeline completed successfully!"
echo "=========================================="
echo "Output directory: $OUTPUT_DIR"
echo ""
echo "To view results:"
echo "  cd $OUTPUT_DIR"
echo "  cat SUMMARY.txt"
echo ""
if [ "$SKIP_CG_PREP" != "true" ]; then
    echo "CG system ready for simulation in: $OUTPUT_DIR/MDRUN_CG/"
    echo "To run energy minimization:"
    echo "  cd $OUTPUT_DIR/MDRUN_CG"
    echo "  gmx mdrun -deffnm em -v"
fi
if [ "$RUN_DISTRIBUTIONS" = "true" ] && [ "$SKIP_ANALYSIS" != "true" ]; then
    echo ""
    echo "Statistics available in: $OUTPUT_DIR/$DIST_OUTPUT_DIR/"
    echo "  - bond_statistics.tsv"
    echo "  - angle_statistics.tsv"
    echo "  - dihedral_statistics.tsv"
fi
echo "=========================================="
