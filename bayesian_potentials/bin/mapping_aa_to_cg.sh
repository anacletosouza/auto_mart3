#!/bin/bash

# =====================================================
# Script: mapping_aa_to_cg.sh
# Description:
# Maps atomistic trajectory to coarse-grained (CG),
# generates topology and runs GROMACS preprocessing.
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
    echo "Required:"
    echo "  --cg_ndx FILE             CG index file (used for both mapping and ITP generation)"
    echo "  --aa_tpr FILE             AA .tpr file"
    echo "  --aa_xtc FILE             AA trajectory"
    echo "  --output_mapped FILE      Output mapped trajectory"
    echo "  --output_cg_gro FILE      Output CG .gro"
    echo "  --aa_itp FILE             AA .itp file"
    echo "  --output_cg_itp FILE      Output CG .itp"
    echo ""
    echo "Optional:"
    echo "  --python_scripts_dir DIR  Python scripts directory (auto-detected)"
    echo "  --input_mdp FILE          MDP file for grompp (skips grompp if not provided)"
    echo "  --path_ff DIR             Force field directory (auto-detects Martini3)"
    echo "  --ff FILE                 Force field ITP filename (default: martini_v3.0.0.itp)"
    echo "  --ions FILE               Ions ITP filename (default: martini_v3.0.0_ions_v1.itp)"
    echo "  --solvent FILE            Solvent ITP filename (default: martini_v3.0.0_solvents_v1.itp)"
    echo "  --name_molecule NAME      Molecule name (default: CG)"
    echo "  --number_molecule NUM     Number of molecules (default: 1)"
    echo "  --title_comments TEXT     Topology comments"
    echo "  --title_system TEXT       System title"
    echo "  --output_topol FILE       Output topology filename (default: topol_cg.top)"
    echo "  --default_martini         Use default Martini 3 masses (72) and zero charges"
    echo "  --maxwarn N               Max warnings for grompp (default: 1)"
    echo "  --output_dir DIR          Output directory (default: results)"
    echo "  --remove_pbc              Remove PBC (default)"
    echo "  --no_pbc                  Skip PBC removal"
    echo "  --skip_grompp             Skip grompp step to generate CG.tpr (if True)"
    echo "  --keep_temp               Keep temporary files"
    echo "  --verbose                 Verbose output"
    exit 1
}

# -----------------------
# Default values
# -----------------------
CG_NDX=""
PYTHON_SCRIPTS_DIR=""
AA_TPR=""
AA_XTC=""
OUTPUT_MAPPED=""
OUTPUT_CG_GRO=""
AA_ITP=""
OUTPUT_CG_ITP=""
INPUT_MDP=""
OUTPUT_TOPOL="topol_cg.top"
OUTPUT_TPR="CG.tpr"

# Topology options
PATH_FF=""
FF="martini_v3.0.0.itp"
IONS="martini_v3.0.0_ions_v1.itp"
SOLVENT="martini_v3.0.0_solvents_v1.itp"
NAME_MOLECULE="CG"
NUMBER_MOLECULE=1
TITLE_COMMENTS="Topology system in Martini 3"
TITLE_SYSTEM="molecule in aqueous solution"
DEFAULT_MARTINI="false"

# Other options
REMOVE_PBC="true"
SKIP_GROMPP="false"
KEEP_TEMP="false"
VERBOSE="false"
MAXWARN=1
OUTPUT_DIR="results"

C_FILE=""
O_FILE=""
P_FILE=""

# -----------------------
# Parse arguments
# -----------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --cg_ndx) CG_NDX="$2"; shift 2 ;;
        --python_scripts_dir) PYTHON_SCRIPTS_DIR="$2"; shift 2 ;;
        --aa_tpr) AA_TPR="$2"; shift 2 ;;
        --aa_xtc) AA_XTC="$2"; shift 2 ;;
        --output_mapped) OUTPUT_MAPPED="$2"; shift 2 ;;
        --remove_pbc) REMOVE_PBC="true"; shift ;;
        --no_pbc) REMOVE_PBC="false"; shift ;;
        --output_cg_gro) OUTPUT_CG_GRO="$2"; shift 2 ;;
        --aa_itp) AA_ITP="$2"; shift 2 ;;
        --output_cg_itp) OUTPUT_CG_ITP="$2"; shift 2 ;;
        --input_mdp) INPUT_MDP="$2"; shift 2 ;;
        --output_topol) OUTPUT_TOPOL="$2"; shift 2 ;;
        --output_tpr) OUTPUT_TPR="$2"; shift 2 ;;
        --skip_grompp) SKIP_GROMPP="true"; shift ;;
        --keep_temp) KEEP_TEMP="true"; shift ;;
        --verbose) VERBOSE="true"; shift ;;
        --default_martini) DEFAULT_MARTINI="true"; shift ;;
        
        # Topology options
        --path_ff) PATH_FF="$2"; shift 2 ;;
        --ff) FF="$2"; shift 2 ;;
        --ions) IONS="$2"; shift 2 ;;
        --solvent) SOLVENT="$2"; shift 2 ;;
        --name_molecule) NAME_MOLECULE="$2"; shift 2 ;;
        --number_molecule) NUMBER_MOLECULE="$2"; shift 2 ;;
        --title_comments) TITLE_COMMENTS="$2"; shift 2 ;;
        --title_system) TITLE_SYSTEM="$2"; shift 2 ;;
        
        --maxwarn) MAXWARN="$2"; shift 2 ;;
        --output_dir) OUTPUT_DIR="$2"; shift 2 ;;
        -c) C_FILE="$2"; shift 2 ;;
        -o) O_FILE="$2"; shift 2 ;;
        -p) P_FILE="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

# -----------------------
# Validation (only essential required args)
# -----------------------
if [ -z "$CG_NDX" ] || \
   [ -z "$AA_TPR" ] || [ -z "$AA_XTC" ] || \
   [ -z "$OUTPUT_MAPPED" ] || [ -z "$OUTPUT_CG_GRO" ] || \
   [ -z "$AA_ITP" ] || [ -z "$OUTPUT_CG_ITP" ]; then
    echo "Error: Missing required arguments"
    echo ""
    usage
fi

# Auto-detect python scripts directory if not provided
if [ -z "$PYTHON_SCRIPTS_DIR" ]; then
    if [ -d "$PACKAGE_DIR/scripts" ]; then
        PYTHON_SCRIPTS_DIR="$PACKAGE_DIR/scripts"
    elif [ -d "$SCRIPT_DIR/../scripts" ]; then
        PYTHON_SCRIPTS_DIR="$(cd "$SCRIPT_DIR/../scripts" && pwd)"
    else
        echo "Error: Could not auto-detect python scripts directory"
        echo "Please provide --python_scripts_dir"
        exit 1
    fi
fi

# Auto-detect force field directory if not provided
if [ -z "$PATH_FF" ]; then
    if [ -d "$PACKAGE_DIR/../data/ff_files" ]; then
        PATH_FF="$(cd "$PACKAGE_DIR/../data/ff_files" && pwd)"
    elif [ -d "$SCRIPT_DIR/../../data/ff_files" ]; then
        PATH_FF="$(cd "$SCRIPT_DIR/../../data/ff_files" && pwd)"
    else
        echo "Warning: Could not auto-detect force field directory"
    fi
fi

# -----------------------
# Convert to absolute paths
# -----------------------
ORIGINAL_DIR="$(pwd)"

get_abs_path() {
    if [[ "$1" = /* ]]; then
        echo "$1"
    else
        echo "$ORIGINAL_DIR/$1"
    fi
}

CG_NDX_ABS=$(get_abs_path "$CG_NDX")
AA_TPR_ABS=$(get_abs_path "$AA_TPR")
AA_XTC_ABS=$(get_abs_path "$AA_XTC")
AA_ITP_ABS=$(get_abs_path "$AA_ITP")

if [ -n "$INPUT_MDP" ]; then
    INPUT_MDP_ABS=$(get_abs_path "$INPUT_MDP")
fi

if [ -n "$PATH_FF" ]; then
    PATH_FF_ABS=$(get_abs_path "$PATH_FF")
fi

# -----------------------
# Prepare output directory
# -----------------------
mkdir -p "$OUTPUT_DIR"
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
echo "Input files:"
echo "  CG index:          $CG_NDX_ABS"
echo "  AA TPR:            $AA_TPR_ABS"
echo "  AA XTC:            $AA_XTC_ABS"
echo "  AA ITP:            $AA_ITP_ABS"
echo ""
echo "Output files:"
echo "  Mapped trajectory: $OUTPUT_MAPPED"
echo "  CG GRO:            $OUTPUT_CG_GRO"
echo "  CG ITP:            $OUTPUT_CG_ITP"
echo "  Output directory:  $OUTPUT_DIR"
echo ""
echo "Options:"
echo "  Remove PBC:        $REMOVE_PBC"
echo "  Skip grompp:       $SKIP_GROMPP"
echo "  Default Martini:   $DEFAULT_MARTINI"
echo "  Verbose:           $VERBOSE"
echo "  Max warnings:      $MAXWARN"
echo ""
if [ -n "$PATH_FF_ABS" ]; then
    echo "Force field:"
    echo "  Directory:        $PATH_FF_ABS"
    echo "  FF file:          $FF"
    echo "  Ions file:        $IONS"
    echo "  Solvent file:     $SOLVENT"
    echo "  Molecule:         $NAME_MOLECULE ($NUMBER_MOLECULE)"
fi
echo "=========================================="

# -----------------------
# Set up Python commands
# -----------------------
# Check if we can use the module directly
if python -c "import bayesian_potentials.scripts.map_aa_to_cg" 2>/dev/null; then
    MAP_CMD="python -m bayesian_potentials.scripts.map_aa_to_cg"
    TOP_CMD="python -m bayesian_potentials.scripts.generate_cg_top"
    log_verbose "Using Python module: bayesian_potentials.scripts"
else
    # Fallback to direct script execution
    MAP_SCRIPT="$PYTHON_SCRIPTS_DIR/1-map_aa_traj_to_cg.py"
    TOP_SCRIPT="$PYTHON_SCRIPTS_DIR/2-obtaining_cg_top.py"
    
    if [ ! -f "$MAP_SCRIPT" ]; then
        echo "Error: Mapping script not found: $MAP_SCRIPT"
        exit 1
    fi
    if [ ! -f "$TOP_SCRIPT" ]; then
        echo "Error: Topology script not found: $TOP_SCRIPT"
        exit 1
    fi
    
    MAP_CMD="python $MAP_SCRIPT"
    TOP_CMD="python $TOP_SCRIPT"
    log_verbose "Using direct scripts: $PYTHON_SCRIPTS_DIR"
fi

# Check for cg-generate-itp command
if command -v cg-generate-itp &> /dev/null; then
    GENERATE_ITP_CMD="cg-generate-itp"
elif python -c "import bayesian_potentials.wrappers; bayesian_potentials.wrappers.generate_itp" 2>/dev/null; then
    GENERATE_ITP_CMD="python -c 'import sys; from bayesian_potentials.wrappers import generate_itp; sys.argv[0]=\"cg-generate-itp\"; generate_itp()'"
else
    GENERATE_ITP_CMD="python -m bayesian_potentials.scripts.4-defining_atoms_type"
fi

# -----------------------
# Step 1: Mapping
# -----------------------
echo ""
echo "Step 1/3: Mapping AA → CG trajectory"

MAP_CMD_ARGS="--index_cg $CG_NDX_ABS --aa_tpr $AA_TPR_ABS --aa_xtc $AA_XTC_ABS --output_mapped $OUTPUT_MAPPED --output_cg_gro $OUTPUT_CG_GRO"

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
echo "  ✓ Generated: $OUTPUT_MAPPED"
echo "  ✓ Generated: $OUTPUT_CG_GRO"

# -----------------------
# Step 2: Generate ITP
# -----------------------
echo ""
echo "Step 2/3: Generating CG ITP"

GENERATE_ITP_ARGS="--cg_gro $OUTPUT_CG_GRO --cg_ndx $CG_NDX_ABS --aa_itp $AA_ITP_ABS --output $OUTPUT_CG_ITP --name_molecule $NAME_MOLECULE"

# Find the mapping JSON file
DEF_JSON=""
if [ -f "$PACKAGE_DIR/../data/definitions_atoms_ff_martini3.json" ]; then
    DEF_JSON="$PACKAGE_DIR/../data/definitions_atoms_ff_martini3.json"
elif [ -f "$SCRIPT_DIR/../../data/definitions_atoms_ff_martini3.json" ]; then
    DEF_JSON="$SCRIPT_DIR/../../data/definitions_atoms_ff_martini3.json"
fi

if [ -n "$DEF_JSON" ]; then
    GENERATE_ITP_ARGS="$GENERATE_ITP_ARGS --def_json $DEF_JSON"
else
    echo "Warning: Could not find definitions_atoms_ff_martini3.json"
fi

if [ "$DEFAULT_MARTINI" = "true" ]; then
    GENERATE_ITP_ARGS="$GENERATE_ITP_ARGS --default"
fi

log_verbose "Running: $GENERATE_ITP_CMD $GENERATE_ITP_ARGS"
eval $GENERATE_ITP_CMD $GENERATE_ITP_ARGS
check_error "ITP generation"
echo "  ✓ Generated: $OUTPUT_CG_ITP"

# -----------------------
# Step 3: Generate topology and run grompp
# -----------------------
if [ "$SKIP_GROMPP" != "true" ] && [ -n "$INPUT_MDP" ] && [ -n "$PATH_FF_ABS" ]; then
    
    echo ""
    echo "Step 3/3: Generating topology and running grompp"
    
    # Build command for topology generation
    TOP_CMD_ARGS="--path_ff $PATH_FF_ABS"
    TOP_CMD_ARGS="$TOP_CMD_ARGS --ff $FF"
    TOP_CMD_ARGS="$TOP_CMD_ARGS --ions $IONS"
    TOP_CMD_ARGS="$TOP_CMD_ARGS --solvent $SOLVENT"
    TOP_CMD_ARGS="$TOP_CMD_ARGS --itp_ligand $OUTPUT_CG_ITP"
    TOP_CMD_ARGS="$TOP_CMD_ARGS --name_molecule $NAME_MOLECULE"
    TOP_CMD_ARGS="$TOP_CMD_ARGS --number_molecule $NUMBER_MOLECULE"
    TOP_CMD_ARGS="$TOP_CMD_ARGS --title_comments \"$TITLE_COMMENTS\""
    TOP_CMD_ARGS="$TOP_CMD_ARGS --title_system \"$TITLE_SYSTEM\""
    TOP_CMD_ARGS="$TOP_CMD_ARGS --output_topol $OUTPUT_TOPOL"
    
    log_verbose "Running: $TOP_CMD $TOP_CMD_ARGS"
    eval $TOP_CMD $TOP_CMD_ARGS
    check_error "Topology generation"
    echo "  ✓ Generated: $OUTPUT_TOPOL"
    
    # Set default files for grompp
    [ -z "$C_FILE" ] && C_FILE="$OUTPUT_CG_GRO"
    [ -z "$P_FILE" ] && P_FILE="$OUTPUT_TOPOL"
    [ -z "$O_FILE" ] && O_FILE="$OUTPUT_TPR"
    
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
    echo "Step 3/3: Skipping grompp (--skip_grompp or missing MDP/FF)"
    echo "  Topology files generated:"
    echo "    • $OUTPUT_CG_ITP - CG topology"
    echo "    • $OUTPUT_CG_GRO - CG coordinates"
    if [ -f "$OUTPUT_TOPOL" ]; then
        echo "    • $OUTPUT_TOPOL - GROMACS topology"
    fi
fi

# -----------------------
# Cleanup temporary files
# -----------------------
if [ "$KEEP_TEMP" != "true" ]; then
    log_verbose "Cleaning up temporary files..."
    find . -name "*.pyc" -type f -delete 2>/dev/null
    find . -name "__pycache__" -type d -delete 2>/dev/null
fi

# -----------------------
# Done
# -----------------------
echo ""
echo "=========================================="
echo "✓ Pipeline completed successfully!"
echo "=========================================="
echo "Output directory: $(pwd)"
echo ""
echo "Generated files:"
echo "  • $OUTPUT_MAPPED - Mapped trajectory"
echo "  • $OUTPUT_CG_GRO - CG coordinates"
echo "  • $OUTPUT_CG_ITP - CG topology"
if [ -f "$OUTPUT_TOPOL" ]; then
    echo "  • $OUTPUT_TOPOL - GROMACS topology"
fi
if [ -n "$O_FILE" ] && [ -f "$O_FILE" ]; then
    echo "  • $O_FILE - GROMACS TPR file"
fi
echo "=========================================="

cd "$ORIGINAL_DIR"
