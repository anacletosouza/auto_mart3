#!/bin/bash

set -e
set -u
set -o pipefail

# =====================================================
# Script: autoparam_CG.sh
# Description:
# - Automated iterative optimization of CG bonded parameters
# - Runs MD simulations, analyzes distributions, updates forces
# - Uses Bayesian update with simulated annealing
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
    echo "  --cg_md_dir DIR           Directory containing CG simulation files"
    echo "  --bonds_ref FILE          Reference bonds TSV file (AA reference)"
    echo "  --angles_ref FILE         Reference angles TSV file (AA reference)"
    echo "  --dihedrals_ref FILE      Reference dihedrals TSV file (AA reference)"
    echo "  --aa_xvg_bonds DIR        Directory with AA bond XVG files (for plotting)"
    echo "  --aa_xvg_angles DIR       Directory with AA angle XVG files (for plotting)"
    echo "  --aa_xvg_dihedrals DIR    Directory with AA dihedral XVG files (for plotting)"
    echo "  --ndx_bonds FILE          NDX file with bond indices"
    echo "  --ndx_angles FILE         NDX file with angle indices"
    echo "  --ndx_dihedrals FILE      NDX file with dihedral indices"
    echo ""
    echo "OPTIONAL ARGUMENTS:"
    echo "  --iterations NUM          Number of iterations (default: 20)"
    echo "  --workdir DIR             Working directory (default: calibration)"
    echo "  --ntomp NUM               Number of OpenMP threads (default: 12)"
    echo "  --ntmpi NUM               Number of MPI threads (default: 1)"
    echo "  --mdp_dir DIR             MDP files directory (required for MD runs)"
    echo "  --solv_ions_gro FILE      Solvated ions GRO filename (default: solv_ions_CG.gro)"
    echo "  --itp_to_optimize FILE    ITP file to optimize (default: cg.itp)"
    echo "  --topol_cg_file FILE      Topology filename (default: topol_cg.top)"
    echo "  --ff_files_dir DIR        Force field files directory (default: ff_files)"
    echo "  --index_ndx FILE          Index file for analysis (default: index.ndx)"
    echo "  --group_out STR           Group for output in analysis (default: System)"
    echo ""
    echo "MULTIMODAL OPTIONS:"
    echo "  --multimodal_mode BOOL    Use peak density mode (default: True)"
    echo "  --variance_multimodal BOOL Use variance from all data (default: True)"
    echo ""
    echo "SIMULATED ANNEALING OPTIONS:"
    echo "  --T0 FLOAT                Initial temperature (default: 10.0)"
    echo "  --alpha FLOAT             Cooling factor (default: 0.85)"
    echo ""
    echo "OUTPUT OPTIONS:"
    echo "  --output_final DIR        Final optimized output directory (default: FINAL_OPTIMIZED)"
    echo "  --keep_temp               Keep temporary files"
    echo "  --skip_md                 Skip MD simulation (use existing XTC files)"
    echo "  --verbose                 Verbose output"
    echo ""
    echo "Examples:"
    echo "  $0 --cg_md_dir MDRUN_CG \\"
    echo "     --bonds_ref STATISTICS/bond_statistics.tsv \\"
    echo "     --angles_ref STATISTICS/angle_statistics.tsv \\"
    echo "     --dihedrals_ref STATISTICS/dihedral_statistics.tsv \\"
    echo "     --aa_xvg_bonds XVG/bonds \\"
    echo "     --aa_xvg_angles XVG/angles \\"
    echo "     --aa_xvg_dihedrals XVG/dihedrals \\"
    echo "     --ndx_bonds NDX/bonds.ndx \\"
    echo "     --ndx_angles NDX/angles.ndx \\"
    echo "     --ndx_dihedrals NDX/dihedrals.ndx \\"
    echo "     --mdp_dir mdp"
    exit 1
}

# -----------------------
# Default values
# -----------------------
# Essential args (no defaults)
CG_MD_DIR=""
BONDS_REF=""
ANGLES_REF=""
DIHEDRALS_REF=""
AA_XVG_BONDS=""
AA_XVG_ANGLES=""
AA_XVG_DIHEDRALS=""
NDX_BONDS=""
NDX_ANGLES=""
NDX_DIHEDRALS=""

# Optional args with defaults
ITERATIONS=20
WORKDIR="calibration"
NTOMP=12
NTMPI=1
MDP_DIR=""
SOLV_IONS_GRO="solv_ions_CG.gro"
ITP_TO_OPTIMIZE="cg.itp"
TOPOL_CG_FILE="topol_cg.top"
FF_FILES_DIR="ff_files"
INDEX_NDX="index.ndx"
GROUP_OUT="System"

# Multimodal options
MULTIMODAL_MODE="True"
VARIANCE_MULTIMODAL="True"

# Simulated annealing options
T0=10.0
ALPHA=0.85

# Output options
OUTPUT_FINAL="FINAL_OPTIMIZED"
KEEP_TEMP="false"
SKIP_MD="false"
VERBOSE="false"

# Internal variables
ORIGINAL_DIR="$(pwd)"

# -----------------------
# Parse arguments
# -----------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --cg_md_dir) CG_MD_DIR="$2"; shift 2 ;;
        --bonds_ref) BONDS_REF="$2"; shift 2 ;;
        --angles_ref) ANGLES_REF="$2"; shift 2 ;;
        --dihedrals_ref) DIHEDRALS_REF="$2"; shift 2 ;;
        --aa_xvg_bonds) AA_XVG_BONDS="$2"; shift 2 ;;
        --aa_xvg_angles) AA_XVG_ANGLES="$2"; shift 2 ;;
        --aa_xvg_dihedrals) AA_XVG_DIHEDRALS="$2"; shift 2 ;;
        --ndx_bonds) NDX_BONDS="$2"; shift 2 ;;
        --ndx_angles) NDX_ANGLES="$2"; shift 2 ;;
        --ndx_dihedrals) NDX_DIHEDRALS="$2"; shift 2 ;;
        
        --iterations) ITERATIONS="$2"; shift 2 ;;
        --workdir) WORKDIR="$2"; shift 2 ;;
        --ntomp) NTOMP="$2"; shift 2 ;;
        --ntmpi) NTMPI="$2"; shift 2 ;;
        --mdp_dir) MDP_DIR="$2"; shift 2 ;;
        --solv_ions_gro) SOLV_IONS_GRO="$2"; shift 2 ;;
        --itp_to_optimize) ITP_TO_OPTIMIZE="$2"; shift 2 ;;
        --topol_cg_file) TOPOL_CG_FILE="$2"; shift 2 ;;
        --ff_files_dir) FF_FILES_DIR="$2"; shift 2 ;;
        --index_ndx) INDEX_NDX="$2"; shift 2 ;;
        --group_out) GROUP_OUT="$2"; shift 2 ;;
        
        --multimodal_mode) MULTIMODAL_MODE="$2"; shift 2 ;;
        --variance_multimodal) VARIANCE_MULTIMODAL="$2"; shift 2 ;;
        
        --T0) T0="$2"; shift 2 ;;
        --alpha) ALPHA="$2"; shift 2 ;;
        
        --output_final) OUTPUT_FINAL="$2"; shift 2 ;;
        --keep_temp) KEEP_TEMP="true"; shift ;;
        --skip_md) SKIP_MD="true"; shift ;;
        --verbose) VERBOSE="true"; shift ;;
        
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

# -----------------------
# Validate essential arguments
# -----------------------
if [ -z "$CG_MD_DIR" ] || [ -z "$BONDS_REF" ] || [ -z "$ANGLES_REF" ] || \
   [ -z "$DIHEDRALS_REF" ] || [ -z "$AA_XVG_BONDS" ] || [ -z "$AA_XVG_ANGLES" ] || \
   [ -z "$AA_XVG_DIHEDRALS" ] || [ -z "$NDX_BONDS" ] || [ -z "$NDX_ANGLES" ] || \
   [ -z "$NDX_DIHEDRALS" ]; then
    echo "Error: Missing essential arguments"
    echo ""
    usage
fi

# Check if MDP directory is provided when not skipping MD
if [ "$SKIP_MD" != "true" ] && [ -z "$MDP_DIR" ]; then
    echo "Error: --mdp_dir is required when running MD simulations"
    echo "Use --skip_md to skip MD simulations and use existing XTC files"
    echo ""
    usage
fi

# -----------------------
# Auto-detect paths
# -----------------------
get_abs_path() {
    if [[ "$1" = /* ]]; then
        echo "$1"
    else
        echo "$ORIGINAL_DIR/$1"
    fi
}

CG_MD_DIR_ABS=$(get_abs_path "$CG_MD_DIR")
BONDS_REF_ABS=$(get_abs_path "$BONDS_REF")
ANGLES_REF_ABS=$(get_abs_path "$ANGLES_REF")
DIHEDRALS_REF_ABS=$(get_abs_path "$DIHEDRALS_REF")
AA_XVG_BONDS_ABS=$(get_abs_path "$AA_XVG_BONDS")
AA_XVG_ANGLES_ABS=$(get_abs_path "$AA_XVG_ANGLES")
AA_XVG_DIHEDRALS_ABS=$(get_abs_path "$AA_XVG_DIHEDRALS")
NDX_BONDS_ABS=$(get_abs_path "$NDX_BONDS")
NDX_ANGLES_ABS=$(get_abs_path "$NDX_ANGLES")
NDX_DIHEDRALS_ABS=$(get_abs_path "$NDX_DIHEDRALS")

# Optional paths
if [ -n "$MDP_DIR" ]; then
    MDP_DIR_ABS=$(get_abs_path "$MDP_DIR")
else
    MDP_DIR_ABS=""
fi

INDEX_NDX_ABS=$(get_abs_path "$INDEX_NDX")

# -----------------------
# Setup working directory
# -----------------------
BASEDIR="$ORIGINAL_DIR"
WORKDIR_ABS="$BASEDIR/$WORKDIR"

# Clean previous workdir if exists
if [ -d "$WORKDIR_ABS" ]; then
    echo "Removing previous working directory: $WORKDIR_ABS"
    rm -rf "$WORKDIR_ABS"
fi

mkdir -p "$WORKDIR_ABS"

# -----------------------
# Export OpenMP settings
# -----------------------
export OMP_NUM_THREADS=$NTOMP
export GMX_GPU_DD_COMMS=true
export GMX_GPU_PME_PP_COMMS=true
export GMX_FORCE_UPDATE_DEFAULT_GPU=true
export GMX_ENABLE_DIRECT_GPU_COMM=TRUE
export GMX_DISABLE_GPU_TIMING=TRUE

# -----------------------
# Display configuration
# -----------------------
echo "=========================================="
echo "CG Parameter Optimization Pipeline"
echo "=========================================="
echo "Essential arguments:"
echo "  CG MD dir:           $CG_MD_DIR_ABS"
echo "  Bonds ref (AA):      $BONDS_REF_ABS"
echo "  Angles ref (AA):     $ANGLES_REF_ABS"
echo "  Dihedrals ref (AA):  $DIHEDRALS_REF_ABS"
echo "  AA XVG bonds:        $AA_XVG_BONDS_ABS"
echo "  AA XVG angles:       $AA_XVG_ANGLES_ABS"
echo "  AA XVG dihedrals:    $AA_XVG_DIHEDRALS_ABS"
echo "  NDX bonds:           $NDX_BONDS_ABS"
echo "  NDX angles:          $NDX_ANGLES_ABS"
echo "  NDX dihedrals:       $NDX_DIHEDRALS_ABS"
echo ""
echo "Optional arguments:"
echo "  Iterations:          $ITERATIONS"
echo "  Working dir:         $WORKDIR_ABS"
echo "  OpenMP threads:      $NTOMP"
echo "  MPI threads:         $NTMPI"
echo "  MDP dir:             $MDP_DIR_ABS"
echo "  Solv ions GRO:       $SOLV_IONS_GRO"
echo "  ITP to optimize:     $ITP_TO_OPTIMIZE"
echo "  Topology file:       $TOPOL_CG_FILE"
echo "  FF files dir:        $FF_FILES_DIR"
echo "  Index NDX:           $INDEX_NDX_ABS"
echo "  Group out:           $GROUP_OUT"
echo ""
echo "Multimodal options:"
echo "  Multimodal mode:     $MULTIMODAL_MODE"
echo "  Variance multimodal: $VARIANCE_MULTIMODAL"
echo ""
echo "Simulated annealing:"
echo "  T0:                  $T0"
echo "  Alpha:               $ALPHA"
echo ""
echo "Output options:"
echo "  Final output:        $OUTPUT_FINAL"
echo "  Keep temp:           $KEEP_TEMP"
echo "  Skip MD:             $SKIP_MD"
echo "  Verbose:             $VERBOSE"
echo "=========================================="

# -----------------------
# Function to check errors
# -----------------------
check_error() {
    if [ $? -ne 0 ]; then
        echo "Error: $1 failed"
        cd "$ORIGINAL_DIR"
        exit 1
    fi
}

# -----------------------
# Verbose output function
# -----------------------
log_verbose() {
    if [ "$VERBOSE" = "true" ]; then
        echo "$1"
    fi
}

# -----------------------
# Function to run MD simulation
# -----------------------
run_md_simulation() {
    local work_dir="$1"
    local mdp_dir="$2"
    local solv_gro="$3"
    local topol_file="$4"
    
    cd "$work_dir"
    echo "    Running MD simulation in: $(pwd)"
    
    # Check if MDP files exist
    if [ ! -d "$mdp_dir" ]; then
        echo "      Error: MDP directory not found: $mdp_dir"
        return 1
    fi
    
    # Energy minimization
    if [ -f "$mdp_dir/minimization.mdp" ]; then
        echo "      Energy minimization..."
        gmx grompp -f "$mdp_dir/minimization.mdp" -p "$topol_file" -c "$solv_gro" -o em.tpr -r "$solv_gro" -maxwarn 2
        check_error "grompp minimization"
        gmx mdrun -deffnm em -v -pin on -ntmpi $NTMPI -ntomp $NTOMP
        check_error "mdrun minimization"
    else
        echo "      Warning: minimization.mdp not found, skipping minimization"
    fi
    
    # NVT equilibration
    if [ -f "$mdp_dir/CG_nvt_1000.mdp" ]; then
        echo "      NVT equilibration..."
        gmx grompp -f "$mdp_dir/CG_nvt_1000.mdp" -p "$topol_file" -c em.gro -o nvt.tpr -r em.gro -maxwarn 2
        check_error "grompp NVT"
        gmx mdrun -deffnm nvt -ntmpi $NTMPI -ntomp $NTOMP -v
        check_error "mdrun NVT"
    else
        echo "      Warning: CG_nvt_1000.mdp not found, skipping NVT"
    fi
    
    # NPT equilibration
    if [ -f "$mdp_dir/CG_npt_1000.mdp" ]; then
        echo "      NPT equilibration..."
        gmx grompp -f "$mdp_dir/CG_npt_1000.mdp" -p "$topol_file" -c nvt.gro -o npt.tpr -r nvt.gro -t nvt.cpt -maxwarn 2
        check_error "grompp NPT"
        gmx mdrun -deffnm npt -ntmpi $NTMPI -ntomp $NTOMP -v
        check_error "mdrun NPT"
    else
        echo "      Warning: CG_npt_1000.mdp not found, skipping NPT"
    fi
    
    # Production MD
    if [ -f "$mdp_dir/CG_md.mdp" ]; then
        echo "      Production MD..."
        gmx grompp -f "$mdp_dir/CG_md.mdp" -p "$topol_file" -c npt.gro -o md.tpr -r npt.gro -t npt.cpt -maxwarn 2
        check_error "grompp MD"
        gmx mdrun -deffnm md -ntmpi $NTMPI -ntomp $NTOMP -v
        check_error "mdrun MD"
    else
        echo "      Warning: CG_md.mdp not found, skipping production MD"
        return 1
    fi
    
    echo "      ✓ MD simulation completed"
    return 0
}

# -----------------------
# Function to analyze trajectory
# -----------------------
analyze_trajectory() {
    local work_dir="$1"
    local bonds_ndx="$2"
    local angles_ndx="$3"
    local dihedrals_ndx="$4"
    local group_out="$5"
    
    cd "$work_dir"
    echo "    Analyzing trajectory..."
    
    # Check if md.xtc exists
    if [ ! -f "md.xtc" ]; then
        echo "      Error: md.xtc not found in $(pwd)"
        return 1
    fi
    
    # Step 1: Generate TPR if not exists
    if [ ! -f "solv_ions.tpr" ]; then
        echo "      Generating TPR file..."
        gmx grompp -f "$MDP_DIR_ABS/minimization.mdp" -c "$SOLV_IONS_GRO" -p "$TOPOL_CG_FILE" -o solv_ions.tpr -maxwarn 2 2>&1 | grep -v "GROMACS reminds you"
        check_error "grompp for TPR"
    fi
    
    # Step 2: Remove PBC with -pbc whole
    echo "      Removing PBC (making molecules whole)..."
    echo "0" | gmx trjconv -s solv_ions.tpr -f md.xtc -o md_whole.xtc -pbc whole 2>&1 | grep -v "GROMACS reminds you" || true
    
    if [ ! -f "md_whole.xtc" ]; then
        echo "      Warning: Failed to create md_whole.xtc, trying alternative method..."
        echo "0" | gmx trjconv -s solv_ions.tpr -f md.xtc -o md_whole.xtc -pbc mol 2>&1 | grep -v "GROMACS reminds you" || true
    fi
    
    # Step 3: Remove jumps with -pbc nojump
    if [ -f "md_whole.xtc" ]; then
        echo "      Removing jumps from trajectory..."
        echo "0" | gmx trjconv -s solv_ions.tpr -f md_whole.xtc -o md_nojump.xtc -pbc nojump 2>&1 | grep -v "GROMACS reminds you" || true
    else
        echo "      Warning: md_whole.xtc not found, trying with original md.xtc"
        echo "0" | gmx trjconv -s solv_ions.tpr -f md.xtc -o md_nojump.xtc -pbc nojump 2>&1 | grep -v "GROMACS reminds you" || true
    fi
    
    # Step 4: Center the molecule in the box
    if [ -f "md_nojump.xtc" ]; then
        echo "      Centering molecule in box..."
        printf "0\n0\n" | gmx trjconv -s solv_ions.tpr -f md_nojump.xtc -o md_center.xtc -center -pbc mol -ur compact 2>&1 | grep -v "GROMACS reminds you" || true
    else
        echo "      Warning: md_nojump.xtc not found, trying with md_whole.xtc"
        printf "0\n0\n" | gmx trjconv -s solv_ions.tpr -f md_whole.xtc -o md_center.xtc -center -pbc mol -ur compact 2>&1 | grep -v "GROMACS reminds you" || true
    fi
    
    # Step 5: Fallback
    if [ ! -f "md_center.xtc" ]; then
        echo "      Warning: Centering failed, using trajectory with PBC removed only"
        if [ -f "md_nojump.xtc" ]; then
            cp md_nojump.xtc md_center.xtc
        elif [ -f "md_whole.xtc" ]; then
            cp md_whole.xtc md_center.xtc
        else
            echo "      Error: No valid trajectory file found"
            return 1
        fi
    fi
    
    # Run analysis using bayesian-potentials analyze command
    echo "    Running bonds/angles/dihedrals analysis..."
    
    CMD="bayesian-potentials analyze \
        --bonds_ndx \"$bonds_ndx\" \
        --angles_ndx \"$angles_ndx\" \
        --dihedrals_ndx \"$dihedrals_ndx\" \
        --xtc_file \"md_center.xtc\" \
        --tpr_file \"solv_ions.tpr\" \
        --group_1 \"$group_out\" \
        --group_2 \"$group_out\""
    
    if [ -f "$INDEX_NDX_ABS" ]; then
        CMD="$CMD --index \"$INDEX_NDX_ABS\""
    fi
    
    eval $CMD
    check_error "bayesian-potentials analyze"
    
    mkdir -p bonds angles dihedrals
    mv bond_*.xvg bonds/ 2>/dev/null || true
    mv ang_*.xvg angles/ 2>/dev/null || true
    mv dih_*.xvg dihedrals/ 2>/dev/null || true
    mv distr_*.xvg bonds/ 2>/dev/null || true
    
    echo "      ✓ Analysis completed"
    
    return 0
}

# -----------------------
# Create initial iteration directory (iter_0)
# -----------------------
echo ""
echo "Setting up initial iteration (iter_0)..."

INITIAL_ITER_DIR="$WORKDIR_ABS/iter_0"
mkdir -p "$INITIAL_ITER_DIR"

if [ -d "$CG_MD_DIR_ABS" ]; then
    echo "  Copying files from $CG_MD_DIR_ABS to $INITIAL_ITER_DIR"
    cp -r "$CG_MD_DIR_ABS"/* "$INITIAL_ITER_DIR/" 2>/dev/null || true
    
    if [ -d "$CG_MD_DIR_ABS/$FF_FILES_DIR" ]; then
        cp -r "$CG_MD_DIR_ABS/$FF_FILES_DIR" "$INITIAL_ITER_DIR/" 2>/dev/null || true
    fi
    
    if [ -n "$MDP_DIR_ABS" ] && [ -d "$MDP_DIR_ABS" ]; then
        cp "$MDP_DIR_ABS"/*.mdp "$INITIAL_ITER_DIR/" 2>/dev/null || true
    fi
    
    echo "  ✓ Initial iteration setup complete"
else
    echo "  Error: CG_MD_DIR not found: $CG_MD_DIR_ABS"
    exit 1
fi

# -----------------------
# Main optimization loop
# -----------------------
echo ""
echo "=========================================="
echo "Starting optimization loop ($ITERATIONS iterations)"
echo "=========================================="

for ((i=0; i<$ITERATIONS; i++)); do
    echo ""
    echo "========================================="
    echo "Iteration $i"
    echo "========================================="

    CURRENT_ITER="iter_$i"
    NEXT_ITER="iter_$((i+1))"
    CURRENT_DIR="$WORKDIR_ABS/$CURRENT_ITER"
    NEXT_DIR="$WORKDIR_ABS/$NEXT_ITER"

    cd "$CURRENT_DIR"
    echo "  Working directory: $(pwd)"

    mkdir -p stat

    # -------------------------------
    # Run MD simulation (if not skipped)
    # -------------------------------
    if [ "$SKIP_MD" != "true" ]; then
        if [ ! -f "md.xtc" ] || [ $i -gt 0 ]; then
            run_md_simulation "$CURRENT_DIR" "$MDP_DIR_ABS" "$SOLV_IONS_GRO" "$TOPOL_CG_FILE"
        else
            echo "  Using existing md.xtc file"
        fi
    else
        echo "  Skipping MD simulation (--skip_md)"
        if [ ! -f "md.xtc" ]; then
            echo "  Error: md.xtc not found and MD simulation skipped"
            exit 1
        fi
    fi
    
    # -------------------------------
    # Analyze trajectory to generate XVG files
    # -------------------------------
    if [ ! -d "bonds" ] || [ ! -d "angles" ] || [ ! -d "dihedrals" ]; then
        analyze_trajectory "$CURRENT_DIR" "$NDX_BONDS_ABS" "$NDX_ANGLES_ABS" "$NDX_DIHEDRALS_ABS" "$GROUP_OUT"
    else
        echo "  Using existing XVG files"
    fi
    
    # -------------------------------
    # Generate TSV statistics from CG XVG files
    # -------------------------------
    echo "  Generating TSV statistics from CG simulation..."
    
    bayesian-potentials distributions \
        --bonds_dir "bonds" \
        --angles_dir "angles" \
        --dihedrals_dir "dihedrals" \
        --dir_to_output "stat" \
        --bond_out "bond_${i}.tsv" \
        --angle_out "angle_${i}.tsv" \
        --dihedral_out "dihedral_${i}.tsv"
    check_error "bayesian-potentials distributions"
    
    mv bond_${i}.tsv stat/ 2>/dev/null || true
    mv angle_${i}.tsv stat/ 2>/dev/null || true
    mv dihedral_${i}.tsv stat/ 2>/dev/null || true
    
    echo "  ✓ CG statistics saved to stat/"
    
    # -------------------------------
    # Plot distributions (AA reference vs CG simulated)
    # -------------------------------
    # Note: plot_distributions.py is not yet available as CLI command
    # This will be added in a future update
    if [ -f "$PACKAGE_DIR/scripts/plot_distributions.py" ]; then
        echo "  Plotting distributions (AA reference vs CG simulated)..."
        mkdir -p figures
        
        python3 "$PACKAGE_DIR/scripts/plot_distributions.py" \
            --bonds_ref_dir "$AA_XVG_BONDS_ABS" \
            --angles_ref_dir "$AA_XVG_ANGLES_ABS" \
            --dihedrals_ref_dir "$AA_XVG_DIHEDRALS_ABS" \
            --bonds_sim_dir "bonds" \
            --angles_sim_dir "angles" \
            --dihedrals_sim_dir "dihedrals" \
            --figures_dir "figures"
        
        echo "  ✓ Distribution plots saved to figures/"
    else
        echo "  Warning: plot_distributions.py not found, skipping plots"
    fi

    # -------------------------------
    # Update potentials (skip last iteration)
    # -------------------------------
    if [ $i -lt $((ITERATIONS-1)) ]; then
        echo "  Updating force constants..."
    
        mkdir -p "$NEXT_DIR"
    
        echo "  Copying base files to $NEXT_DIR..."
        
        # Copy all files except stat, figures, bonds, angles, dihedrals, and the ITP
        rsync -av --exclude="stat" --exclude="figures" --exclude="bonds" --exclude="angles" --exclude="dihedrals" --exclude="$ITP_TO_OPTIMIZE" \
            "$CURRENT_DIR/" "$NEXT_DIR/" 2>/dev/null || \
        find "$CURRENT_DIR" -maxdepth 1 -type f ! -name "$ITP_TO_OPTIMIZE" -exec cp {} "$NEXT_DIR/" \; 2>/dev/null || true
    
        for dir in ff_files ions_mdp; do
            if [ -d "$CURRENT_DIR/$dir" ]; then
                cp -r "$CURRENT_DIR/$dir" "$NEXT_DIR/" 2>/dev/null || true
            fi
        done
    
        echo "  Running force adjustment..."
        bayesian-potentials force-adjust \
            --bonds_ref "$BONDS_REF_ABS" \
            --angles_ref "$ANGLES_REF_ABS" \
            --dihedrals_ref "$DIHEDRALS_REF_ABS" \
            --bonds_sim "stat/bond_${i}.tsv" \
            --angles_sim "stat/angle_${i}.tsv" \
            --dihedrals_sim "stat/dihedral_${i}.tsv" \
            --itp_cg "$CURRENT_DIR/$ITP_TO_OPTIMIZE" \
            --ndx_bounds "$NDX_BONDS_ABS" \
            --ndx_angles "$NDX_ANGLES_ABS" \
            --ndx_dihedrals "$NDX_DIHEDRALS_ABS" \
            --molecule_name "molecule" \
            --multimodal_mode "$MULTIMODAL_MODE" \
            --variance_multimodal "$VARIANCE_MULTIMODAL" \
            --T0 "$T0" \
            --alpha "$ALPHA" \
            --itp_out "$NEXT_DIR/$ITP_TO_OPTIMIZE"
        
        if [ $? -eq 0 ]; then
            echo "  ✓ Updated ITP saved to: $NEXT_DIR/$ITP_TO_OPTIMIZE"
            
            if [ -f "$CURRENT_DIR/$ITP_TO_OPTIMIZE" ] && [ -f "$NEXT_DIR/$ITP_TO_OPTIMIZE" ]; then
                if diff -q "$CURRENT_DIR/$ITP_TO_OPTIMIZE" "$NEXT_DIR/$ITP_TO_OPTIMIZE" > /dev/null 2>&1; then
                    echo "  ⚠ WARNING: New ITP is identical to old ITP!"
                else
                    echo "  ✓ New ITP is different from old ITP"
                fi
            fi
        else
            echo "  Error: bayesian-potentials force-adjust failed"
            exit 1
        fi
    fi
    
    # Clean up temporary files if not keeping
    if [ "$KEEP_TEMP" != "true" ]; then
        echo "  Cleaning up temporary files..."
        rm -f *.xtc *.trr nvt* npt* md* *.cpt *.edr *.log \#* 2>/dev/null || true
        rm -rf \#* 2>/dev/null || true
    fi
    
    echo "  ✓ Iteration $i completed"

done

# -----------------------
# Save final optimized files
# -----------------------
echo ""
echo "=========================================="
echo "Saving final optimized files"
echo "=========================================="

FINAL_DIR="$BASEDIR/$OUTPUT_FINAL"
mkdir -p "$FINAL_DIR"

FINAL_ITER_DIR="$WORKDIR_ABS/iter_$((ITERATIONS-1))"
if [ -d "$FINAL_ITER_DIR" ]; then
    echo "  Copying final results from $FINAL_ITER_DIR to $FINAL_DIR"
    cp -r "$FINAL_ITER_DIR"/* "$FINAL_DIR/" 2>/dev/null || true
    
    if [ -f "$FINAL_ITER_DIR/$ITP_TO_OPTIMIZE" ]; then
        cp "$FINAL_ITER_DIR/$ITP_TO_OPTIMIZE" "$FINAL_DIR/${ITP_TO_OPTIMIZE%.itp}_optimized.itp"
        echo "  ✓ Optimized ITP saved to: $FINAL_DIR/${ITP_TO_OPTIMIZE%.itp}_optimized.itp"
    fi
fi

# Copy statistics from all iterations
mkdir -p "$FINAL_DIR/all_statistics"
for j in $(seq 0 $((ITERATIONS-1))); do
    if [ -f "$WORKDIR_ABS/iter_$j/stat/bond_${j}.tsv" ]; then
        cp "$WORKDIR_ABS/iter_$j/stat/bond_${j}.tsv" "$FINAL_DIR/all_statistics/" 2>/dev/null || true
    fi
    if [ -f "$WORKDIR_ABS/iter_$j/stat/angle_${j}.tsv" ]; then
        cp "$WORKDIR_ABS/iter_$j/stat/angle_${j}.tsv" "$FINAL_DIR/all_statistics/" 2>/dev/null || true
    fi
    if [ -f "$WORKDIR_ABS/iter_$j/stat/dihedral_${j}.tsv" ]; then
        cp "$WORKDIR_ABS/iter_$j/stat/dihedral_${j}.tsv" "$FINAL_DIR/all_statistics/" 2>/dev/null || true
    fi
done

# Create summary file
cat > "$FINAL_DIR/OPTIMIZATION_SUMMARY.txt" << EOF
==========================================
CG Parameter Optimization Summary
==========================================

Date: $(date)
Working directory: $WORKDIR_ABS
Number of iterations: $ITERATIONS

Input files:
  CG MD directory:       $CG_MD_DIR_ABS
  Bonds reference (AA):  $BONDS_REF_ABS
  Angles reference (AA): $ANGLES_REF_ABS
  Dihedrals reference (AA): $DIHEDRALS_REF_ABS
  AA XVG bonds:          $AA_XVG_BONDS_ABS
  AA XVG angles:         $AA_XVG_ANGLES_ABS
  AA XVG dihedrals:      $AA_XVG_DIHEDRALS_ABS
  MDP directory:         $MDP_DIR_ABS

Optimization parameters:
  Multimodal mode:       $MULTIMODAL_MODE
  Variance multimodal:   $VARIANCE_MULTIMODAL
  T0 (SA):              $T0
  Alpha (SA):           $ALPHA

Output files:
  Optimized ITP:         ${ITP_TO_OPTIMIZE%.itp}_optimized.itp
  Final iteration:       $FINAL_ITER_DIR
  All statistics:        $FINAL_DIR/all_statistics/

Simulation parameters:
  OpenMP threads:        $NTOMP
  MPI threads:           $NTMPI
  Skip MD:               $SKIP_MD
  Group out:             $GROUP_OUT

==========================================
EOF

echo "  ✓ Summary saved to: $FINAL_DIR/OPTIMIZATION_SUMMARY.txt"

# -----------------------
# Cleanup if not keeping temp
# -----------------------
if [ "$KEEP_TEMP" != "true" ]; then
    echo ""
    echo "Cleaning up temporary files..."
    rm -rf "$WORKDIR_ABS" 2>/dev/null || true
    echo "  ✓ Removed working directory: $WORKDIR_ABS"
fi

# -----------------------
# Final message
# -----------------------
echo ""
echo "=========================================="
echo "Optimization completed successfully!"
echo "=========================================="
echo ""
echo "Final optimized files saved in:"
echo "  $FINAL_DIR"
echo ""
echo "Optimized ITP file:"
echo "  $FINAL_DIR/${ITP_TO_OPTIMIZE%.itp}_optimized.itp"
echo ""
echo "Statistics from all iterations:"
echo "  $FINAL_DIR/all_statistics/"
echo ""
echo "To view summary:"
echo "  cat $FINAL_DIR/OPTIMIZATION_SUMMARY.txt"
echo ""
echo "=========================================="

cd "$ORIGINAL_DIR"
