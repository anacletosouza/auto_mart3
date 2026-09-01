#!/bin/bash

# =============================================================================
# SCRIPT: auto_mart_CG.sh
# DESCRIPTION: Automated parameter optimization for Coarse-Grained (CG) 
#              molecular dynamics simulations using Bayesian potentials.
#              Performs iterative optimization of bonds, angles, and dihedrals
#              force constants by comparing reference distributions from AA
#              simulations with CG simulations.
#
# USAGE:
#   ./auto_mart_CG.sh --INPUT_AUTO_MART_AA_DIR <path> --OUTPUT_AUTO_MART_CG_DIR <path> [OPTIONS]
#
# REQUIRED ARGUMENTS:
#   --INPUT_AUTO_MART_AA_DIR <path>     Directory containing folders from auto_mart_aa results (BONDS_ANGLES_DIHEDRALS_XVG_REF  CG_MARTINI3  GMX  MDRUN_CG  STATISTICS)
#   --OUTPUT_AUTO_MART_CG_DIR <path>    Directory for CG optimization output
#
# OPTIONAL ARGUMENTS:
#   Simulation Parameters:
#     --ntomp <int>                     OpenMP threads (default: 10)
#     --ntmpi <int>                     MPI threads (default: 1)
#     --ref_t <float>                   Reference temperature in K (default: 310)
#     --ref_p <float>                   Reference pressure in bar (default: 1.0)
#     --dt_nvt_ps <float>               NVT timestep in ps (default: 0.001)
#     --time_nvt_ps <float>             NVT duration in ps (default: 1000)
#     --dt_npt_ps <float>               NPT timestep in ps (default: 0.001)
#     --time_npt_ps <float>             NPT duration in ps (default: 5000)
#     --dt_md_ps <float>                MD production timestep in ps (default: 0.001)
#     --time_md_ps <float>              MD production duration in ps (default: 10000)
#
#   GROMACS MDRUN Options:
#     --pin_on                          Enable thread pinning (default: on)
#     --pin_off                         Disable thread pinning
#     --nb_gpu                          Use GPU for non-bonded interactions (default: auto)
#     --nb_cpu                          Use CPU for non-bonded interactions
#     --nb_none                         Auto-detect for non-bonded interactions
#     --pme_gpu                         Use GPU for PME (default: auto)
#     --pme_cpu                         Use CPU for PME
#     --pme_none                        Auto-detect for PME
#     --bonded_gpu                      Use GPU for bonded interactions (default: auto)
#     --bonded_cpu                      Use CPU for bonded interactions
#     --bonded_none                     Auto-detect for bonded interactions
#     --npme <int>                      Number of PME threads (default: auto)
#     --cuda_visible_devices <str>      CUDA visible devices (default: 0,1,2,3)
#     --gpu_id <str>                    Alias for CUDA_VISIBLE_DEVICES
#
#   File Names:
#     --topol_cg_file <name>            CG topology file name (default: topol_cg.top)
#     --solv_ions_gro <name>            Solvated ions GRO file (default: solv_ions_CG.gro)
#     --em_tpr <name>                   Energy minimization TPR file (default: em.tpr)
#     --itp_to_optimize <name>          ITP file to optimize (default: cg.itp)
#
#   Force Constants Limits:
#     --min_force_bond <float>          Minimum bond force constant (default: 500.0)
#     --max_force_bond <float>          Maximum bond force constant (default: 50000.0)
#     --default_force_bond <float>      Default bond force constant (default: 1250.0)
#     --min_force_angle <float>         Minimum angle force constant (default: 10.0)
#     --max_force_angle <float>         Maximum angle force constant (default: 150.0)
#     --default_force_angle <float>     Default angle force constant (default: 25.0)
#     --min_force_dihedral <float>      Minimum dihedral force constant (default: 10.0)
#     --max_force_dihedral <float>      Maximum dihedral force constant (default: 150.0)
#     --default_force_dihedral <float>  Default dihedral force constant (default: 25.0)
#
#   Optimization Parameters:
#     --T0 <float>                      Initial simulated annealing temperature (default: 10.0)
#     --alpha <float>                   Cooling factor (default: 0.85)
#     --distribution_points <int>       Points for R² calculation (default: 100)
#     --n_iter <int>                    Number of optimization iterations (default: 30)
#
#   Input/Output Paths:
#     --bonds_ndx <path>                Bonds index file (default: INPUT_AUTO_MART_AA_DIR/CG_MARTINI3/NDX/bonds.ndx)
#     --angles_ndx <path>               Angles index file (default: INPUT_AUTO_MART_AA_DIR/CG_MARTINI3/NDX/angles.ndx)
#     --dihedrals_ndx <path>            Dihedrals index file (default: INPUT_AUTO_MART_AA_DIR/CG_MARTINI3/NDX/dihedrals.ndx)
#     --bonds_ref_xvg_dir <path>        Bonds reference XVG directory (default: INPUT_AUTO_MART_AA_DIR/BONDS_ANGLES_DIHEDRALS_XVG_REF/bonds)
#     --angles_ref_xvg_dir <path>       Angles reference XVG directory (default: INPUT_AUTO_MART_AA_DIR/BONDS_ANGLES_DIHEDRALS_XVG_REF/angles)
#     --dihedrals_ref_xvg_dir <path>    Dihedrals reference XVG directory (default: INPUT_AUTO_MART_AA_DIR/BONDS_ANGLES_DIHEDRALS_XVG_REF/dihedrals)
#
#   Prefixes for XVG Files:
#     --prefix_xvg_bond_ref <str>       Prefix for bond reference XVG files (default: bond_)
#     --prefix_xvg_angle_ref <str>      Prefix for angle reference XVG files (default: ang_)
#     --prefix_xvg_dihedral_ref <str>   Prefix for dihedral reference XVG files (default: dih_)
#     --prefix_xvg_bond_sim <str>       Prefix for bond simulation XVG files (default: bond_)
#     --prefix_xvg_angle_sim <str>      Prefix for angle simulation XVG files (default: ang_)
#     --prefix_xvg_dihedral_sim <str>   Prefix for dihedral simulation XVG files (default: dih_)
#
#   GROMACS Processing:
#     --group_2 <name>                  Output group for trajectory processing (default: System)
#     --group_1 <name>                  Center group for trajectory processing (default: System)
#     --index <path>                    Index file for GROMACS selections (default: None)
#     --no_remove_gmx_files_in_iter     Preserve GMX temporary files (*.xtc, *.edr, etc.)
#
#   Other:
#     --molecule_name <name>            Molecule name for potential adjustment (default: molecule)
#
# EXAMPLES:
#   Basic usage (only required arguments):
#     ./auto_mart_CG.sh \
#         --INPUT_AUTO_MART_AA_DIR "/path/to/AA_results" \
#         --OUTPUT_AUTO_MART_CG_DIR "/path/to/CG_results"
#
#   GPU-accelerated simulation with specific devices:
#     ./auto_mart_CG.sh \
#         --INPUT_AUTO_MART_AA_DIR "/path/to/AA_results" \
#         --OUTPUT_AUTO_MART_CG_DIR "/path/to/CG_results" \
#         --cuda_visible_devices "0,1,2,3" \
#         --nb_gpu --pme_gpu --bonded_gpu \
#         --pin_on --npme 4
#
#   CPU-only simulation:
#     ./auto_mart_CG.sh \
#         --INPUT_AUTO_MART_AA_DIR "/path/to/AA_results" \
#         --OUTPUT_AUTO_MART_CG_DIR "/path/to/CG_results" \
#         --nb_cpu --pme_cpu --bonded_cpu \
#         --pin_on --ntomp 20
#
#   Full customization:
#     ./auto_mart_CG.sh \
#         --INPUT_AUTO_MART_AA_DIR "/path/to/AA_results" \
#         --OUTPUT_AUTO_MART_CG_DIR "/path/to/CG_results" \
#         --ntomp 20 --ntmpi 2 \
#         --ref_t 300 --ref_p 1.01325 \
#         --time_nvt_ps 2000 --time_npt_ps 10000 --time_md_ps 20000 \
#         --n_iter 50 --T0 20.0 --alpha 0.9 \
#         --group_2 "Protein" --group_1 "Protein" \
#         --index "custom_index.ndx" \
#         --no_remove_gmx_files_in_iter \
#         --molecule_name "my_protein"
#
# NOTES:
#   - The script requires GROMACS and bayesian_potentials to be installed
#   - Reference data from AA simulations must be present in INPUT_AUTO_MART_AA_DIR
#   - The optimization loop runs iteratively to match CG distributions to AA references
#   - Use --no_remove_gmx_files_in_iter for debugging trajectory issues
#   - The --index option allows using custom index files for group selections in GROMACS
#
# AUTHOR: Anacleto Silva de Souza
# DATE: 22-04-2026
# =============================================================================

# =========================
# Helper function to convert relative path to absolute path
# =========================
to_absolute_path() {
    local path="$1"
    if [[ -z "$path" ]]; then
        echo ""
        return
    fi
    # If path is already absolute (starts with /), return as is
    if [[ "$path" =~ ^/ ]]; then
        echo "$path"
    else
        # Convert relative path to absolute using current directory
        echo "$(cd "$(dirname "$path")" 2>/dev/null && pwd)/$(basename "$path")" 2>/dev/null || echo "$(pwd)/$path"
    fi
}

# =========================
# Function to build mdrun command with GPU/CPU options
# =========================
build_mdrun_cmd() {
    local deffnm="$1"
    local extra_opts=""
    
    # Add pinning option
    if [ "$PIN_ON" = true ]; then
        extra_opts="$extra_opts -pin on"
    fi
    
    # Add non-bonded option
    if [ -n "$NB_OPTION" ]; then
        extra_opts="$extra_opts -nb $NB_OPTION"
    fi
    
    # Add PME option
    if [ -n "$PME_OPTION" ]; then
        extra_opts="$extra_opts -pme $PME_OPTION"
    fi
    
    # Add bonded option
    if [ -n "$BONDED_OPTION" ]; then
        extra_opts="$extra_opts -bonded $BONDED_OPTION"
    fi
    
    # Add NPME option
    if [ -n "$NPME" ] && [ "$NPME" != "auto" ] && [ "$NPME" != "None" ]; then
        extra_opts="$extra_opts -npme $NPME"
    fi
    
    # Build the complete command
    if [ -n "$CUDA_VISIBLE_DEVICES_SET" ]; then
        echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES_SET gmx mdrun -deffnm $deffnm -v -ntmpi $NTMPI -ntomp $NTOMP $extra_opts"
    else
        echo "gmx mdrun -deffnm $deffnm -v -ntmpi $NTMPI -ntomp $NTOMP $extra_opts"
    fi
}

# =========================
# Parse command line arguments (all arguments use underscores)
# =========================

# Default values
NO_REMOVE_GMX_FILES=false
GROUP_2="System"
GROUP_1="System"
INDEX_FILE="None"

# GROMACS mdrun default values
PIN_ON=true
NB_OPTION=""
PME_OPTION=""
BONDED_OPTION=""
NPME=""
CUDA_VISIBLE_DEVICES_SET=""

# Parse arguments - ALL ARGUMENTS MUST USE UNDERSCORES
while [[ $# -gt 0 ]]; do
    case $1 in
        --no_remove_gmx_files_in_iter)
            NO_REMOVE_GMX_FILES=true
            shift
            ;;
        --group_2)
            GROUP_2="$2"
            shift 2
            ;;
        --group_1)
            GROUP_1="$2"
            shift 2
            ;;
        --index)
            INDEX_FILE="$2"
            shift 2
            ;;
        --INPUT_AUTO_MART_AA_DIR)
            INPUT_AUTO_MART_AA_DIR="$2"
            shift 2
            ;;
        --OUTPUT_AUTO_MART_CG_DIR)
            OUTPUT_AUTO_MART_CG_DIR="$2"
            shift 2
            ;;
        --bonds_ndx)
            BONDS_NDX="$2"
            shift 2
            ;;
        --angles_ndx)
            ANGLES_NDX="$2"
            shift 2
            ;;
        --dihedrals_ndx)
            DIHEDRALS_NDX="$2"
            shift 2
            ;;
        --bonds_ref_xvg_dir)
            BONDS_REF_XVG_DIR="$2"
            shift 2
            ;;
        --angles_ref_xvg_dir)
            ANGLES_REF_XVG_DIR="$2"
            shift 2
            ;;
        --dihedrals_ref_xvg_dir)
            DIHEDRALS_REF_XVG_DIR="$2"
            shift 2
            ;;
        --molecule_name)
            MOLECULE_NAME="$2"
            shift 2
            ;;
        --T0)
            T0="$2"
            shift 2
            ;;
        --alpha)
            ALPHA="$2"
            shift 2
            ;;
        --distribution_points)
            DISTRIBUTION_POINTS="$2"
            shift 2
            ;;
        --ntomp)
            NTOMP="$2"
            shift 2
            ;;
        --ntmpi)
            NTMPI="$2"
            shift 2
            ;;
        --ref_t)
            REF_T="$2"
            shift 2
            ;;
        --ref_p)
            REF_P="$2"
            shift 2
            ;;
        --dt_nvt_ps)
            DT_NVT_PS="$2"
            shift 2
            ;;
        --time_nvt_ps)
            TIME_NVT_PS="$2"
            shift 2
            ;;
        --dt_npt_ps)
            DT_NPT_PS="$2"
            shift 2
            ;;
        --time_npt_ps)
            TIME_NPT_PS="$2"
            shift 2
            ;;
        --dt_md_ps)
            DT_MD_PS="$2"
            shift 2
            ;;
        --time_md_ps)
            TIME_MD_PS="$2"
            shift 2
            ;;
        --topol_cg_file)
            TOPOL_CG_FILE="$2"
            shift 2
            ;;
        --solv_ions_gro)
            SOLV_IONS_GRO="$2"
            shift 2
            ;;
        --em_tpr)
            EM_TPR="$2"
            shift 2
            ;;
        --min_force_bond)
            MIN_FORCE_BOND="$2"
            shift 2
            ;;
        --max_force_bond)
            MAX_FORCE_BOND="$2"
            shift 2
            ;;
        --default_force_bond)
            DEFAULT_FORCE_BOND="$2"
            shift 2
            ;;
        --min_force_angle)
            MIN_FORCE_ANGLE="$2"
            shift 2
            ;;
        --max_force_angle)
            MAX_FORCE_ANGLE="$2"
            shift 2
            ;;
        --default_force_angle)
            DEFAULT_FORCE_ANGLE="$2"
            shift 2
            ;;
        --min_force_dihedral)
            MIN_FORCE_DIHEDRAL="$2"
            shift 2
            ;;
        --max_force_dihedral)
            MAX_FORCE_DIHEDRAL="$2"
            shift 2
            ;;
        --default_force_dihedral)
            DEFAULT_FORCE_DIHEDRAL="$2"
            shift 2
            ;;
        --prefix_xvg_bond_ref)
            PREFIX_XVG_BOND_REF="$2"
            shift 2
            ;;
        --prefix_xvg_angle_ref)
            PREFIX_XVG_ANGLE_REF="$2"
            shift 2
            ;;
        --prefix_xvg_dihedral_ref)
            PREFIX_XVG_DIHEDRAL_REF="$2"
            shift 2
            ;;
        --prefix_xvg_bond_sim)
            PREFIX_XVG_BOND_SIM="$2"
            shift 2
            ;;
        --prefix_xvg_angle_sim)
            PREFIX_XVG_ANGLE_SIM="$2"
            shift 2
            ;;
        --prefix_xvg_dihedral_sim)
            PREFIX_XVG_DIHEDRAL_SIM="$2"
            shift 2
            ;;
        --itp_to_optimize)
            ITP_TO_OPTIMIZE="$2"
            shift 2
            ;;
        --n_iter)
            N_ITER="$2"
            shift 2
            ;;
        # GROMACS mdrun options
        --pin_on)
            PIN_ON=true
            shift
            ;;
        --pin_off)
            PIN_ON=false
            shift
            ;;
        --nb_gpu)
            NB_OPTION="gpu"
            shift
            ;;
        --nb_cpu)
            NB_OPTION="cpu"
            shift
            ;;
        --nb_none)
            NB_OPTION=""
            shift
            ;;
        --pme_gpu)
            PME_OPTION="gpu"
            shift
            ;;
        --pme_cpu)
            PME_OPTION="cpu"
            shift
            ;;
        --pme_none)
            PME_OPTION=""
            shift
            ;;
        --bonded_gpu)
            BONDED_OPTION="gpu"
            shift
            ;;
        --bonded_cpu)
            BONDED_OPTION="cpu"
            shift
            ;;
        --bonded_none)
            BONDED_OPTION=""
            shift
            ;;
        --npme)
            NPME="$2"
            shift 2
            ;;
        --cuda_visible_devices)
            CUDA_VISIBLE_DEVICES_SET="$2"
            shift 2
            ;;
        --gpu_id)
            CUDA_VISIBLE_DEVICES_SET="$2"
            shift 2
            ;;
        --help|-h)
            sed -n '2,/^# =============================================================================/p' "$0" | sed 's/^# //'
            exit 0
            ;;
        *)
            echo "ERROR: Unknown option: $1"
            echo "All options must use underscores (_) instead of hyphens (-)"
            echo "Example: --INPUT_AUTO_MART_AA_DIR instead of --INPUT-AUTO-MART-AA-DIR"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [ -z "$INPUT_AUTO_MART_AA_DIR" ] || [ -z "$OUTPUT_AUTO_MART_CG_DIR" ]; then
    echo "ERROR: --INPUT_AUTO_MART_AA_DIR and --OUTPUT_AUTO_MART_CG_DIR are required"
    echo "Use --help for usage information"
    exit 1
fi

# Convert paths to absolute paths
INPUT_AUTO_MART_AA_DIR=$(to_absolute_path "$INPUT_AUTO_MART_AA_DIR")
OUTPUT_AUTO_MART_CG_DIR=$(to_absolute_path "$OUTPUT_AUTO_MART_CG_DIR")

echo "=== Path Information ==="
echo "INPUT_AUTO_MART_AA_DIR (absolute): $INPUT_AUTO_MART_AA_DIR"
echo "OUTPUT_AUTO_MART_CG_DIR (absolute): $OUTPUT_AUTO_MART_CG_DIR"
echo ""

# Set default values for optional parameters
: ${NTOMP:=12}
: ${NTMPI:=1}
: ${REF_T:=310}
: ${REF_P:=1.0}
: ${DT_NVT_PS:=0.0005}
: ${TIME_NVT_PS:=10000}
: ${DT_NPT_PS:=0.0005}
: ${TIME_NPT_PS:=5000}
: ${DT_MD_PS:=0.0005}
: ${TIME_MD_PS:=100000}
: ${TOPOL_CG_FILE:="topol_cg.top"}
: ${SOLV_IONS_GRO:="solv_ions_CG.gro"}
: ${EM_TPR:="em.tpr"}
: ${MIN_FORCE_BOND:=500.0}
: ${MAX_FORCE_BOND:=50000.0}
: ${DEFAULT_FORCE_BOND:=1250.0}
: ${MIN_FORCE_ANGLE:=10.0}
: ${MAX_FORCE_ANGLE:=1000.0}
: ${DEFAULT_FORCE_ANGLE:=25.0}
: ${MIN_FORCE_DIHEDRAL:=10.0}
: ${MAX_FORCE_DIHEDRAL:=1000.0}
: ${DEFAULT_FORCE_DIHEDRAL:=25.0}
: ${T0:=10.0}
: ${ALPHA:=0.85}
: ${DISTRIBUTION_POINTS:=200}
: ${PREFIX_XVG_BOND_REF:="bond_"}
: ${PREFIX_XVG_ANGLE_REF:="ang_"}
: ${PREFIX_XVG_DIHEDRAL_REF:="dih_"}
: ${PREFIX_XVG_BOND_SIM:="bond_"}
: ${PREFIX_XVG_ANGLE_SIM:="ang_"}
: ${PREFIX_XVG_DIHEDRAL_SIM:="dih_"}
: ${ITP_TO_OPTIMIZE:="cg.itp"}
: ${N_ITER:=50}
: ${BONDS_NDX:="$INPUT_AUTO_MART_AA_DIR/CG_MARTINI3/NDX/bonds.ndx"}
: ${ANGLES_NDX:="$INPUT_AUTO_MART_AA_DIR/CG_MARTINI3/NDX/angles.ndx"}
: ${DIHEDRALS_NDX:="$INPUT_AUTO_MART_AA_DIR/CG_MARTINI3/NDX/dihedrals.ndx"}
: ${BONDS_REF_XVG_DIR:="$INPUT_AUTO_MART_AA_DIR/BONDS_ANGLES_DIHEDRALS_XVG_REF/bonds"}
: ${ANGLES_REF_XVG_DIR:="$INPUT_AUTO_MART_AA_DIR/BONDS_ANGLES_DIHEDRALS_XVG_REF/angles"}
: ${DIHEDRALS_REF_XVG_DIR:="$INPUT_AUTO_MART_AA_DIR/BONDS_ANGLES_DIHEDRALS_XVG_REF/dihedrals"}
: ${MOLECULE_NAME:="molecule"}

# Convert additional paths to absolute paths
BONDS_NDX=$(to_absolute_path "$BONDS_NDX")
ANGLES_NDX=$(to_absolute_path "$ANGLES_NDX")
DIHEDRALS_NDX=$(to_absolute_path "$DIHEDRALS_NDX")
BONDS_REF_XVG_DIR=$(to_absolute_path "$BONDS_REF_XVG_DIR")
ANGLES_REF_XVG_DIR=$(to_absolute_path "$ANGLES_REF_XVG_DIR")
DIHEDRALS_REF_XVG_DIR=$(to_absolute_path "$DIHEDRALS_REF_XVG_DIR")

if [ "$INDEX_FILE" != "None" ]; then
    INDEX_FILE=$(to_absolute_path "$INDEX_FILE")
fi

# Display GROMACS mdrun configuration
echo "=== GROMACS mdrun Configuration ==="
echo "Thread pinning: $PIN_ON"
echo "Non-bonded (nb): ${NB_OPTION:-auto}"
echo "PME: ${PME_OPTION:-auto}"
echo "Bonded: ${BONDED_OPTION:-auto}"
echo "PME threads (npme): ${NPME:-auto}"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES_SET:-not set}"
echo "OpenMP threads (ntomp): $NTOMP"
echo "MPI threads (ntmpi): $NTMPI"
echo ""

# Function to build GROMACS command with optional index file
build_gmx_cmd() {
    local cmd="$1"
    local input="$2"
    local output="$3"
    local options="$4"
    
    if [ "$INDEX_FILE" != "None" ] && [ -f "$INDEX_FILE" ]; then
        eval "$cmd -n $INDEX_FILE $options $input $output"
    else
        eval "$cmd $options $input $output"
    fi
}

# Function to run mdrun with GPU/CPU options
run_mdrun() {
    local deffnm="$1"
    local mdrun_cmd=$(build_mdrun_cmd "$deffnm")
    
    echo "Running: $mdrun_cmd"
    eval "$mdrun_cmd"
}

rm -fr AUTO_MART_CG_RESULTS

# ------------------------------------------------------------------------------
# STEP 1: Creating folder for output pipeline
# ------------------------------------------------------------------------------

i=0

mkdir -p "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/MD_CG"
mkdir -p "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/RESULTS/BONDS_ANGLES_DIHEDRALS_XVG_SIM"/{bonds,angles,dihedrals}
mkdir -p "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/RESULTS/FIGURES_DISTR_SERIES_REF_VS_SIM"
   
# Copy initial CG system files
if [ -d "$INPUT_AUTO_MART_AA_DIR/MDRUN_CG" ]; then
    cp -r "$INPUT_AUTO_MART_AA_DIR/MDRUN_CG"/* "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/MD_CG" 2>/dev/null || true
fi
    
cd "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/MD_CG"
echo "Working in: $(pwd)"
    
# ------------------------------------------------------------------------------
# STEP 1: CREATE MDP FILES FOR SIMULATION
# ------------------------------------------------------------------------------
echo ""
echo "=== STEP 1: Creating MDP files ==="
    
# Create minimization.mdp
cat > minimization.mdp << 'EOF'
; Energy minimization for Coarse-Grained system
integrator               = steep
emtol                    = 100.0
emstep                   = 0.01
nsteps                   = 50000
nstcgsteep               = 1000
nstlog                   = 1000
nstenergy                = 1000
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstxout-compressed       = 0
cutoff-scheme            = Verlet
nstlist                  = 20
rlist                    = 1.35
coulombtype              = reaction-field
rcoulomb                 = 1.1
epsilon_r                = 15
epsilon_rf               = 0
vdw_type                 = cutoff
vdw-modifier             = Potential-shift-verlet
rvdw                     = 1.1
pbc                      = xyz
tcoupl                   = no
pcoupl                   = no
constraints              = none
EOF
    
# Create nvt.mdp
NSTEPS_NVT=$(echo "$TIME_NVT_PS / $DT_NVT_PS" | bc -l | awk '{print int($1)}')
    
cat > nvt.mdp << EOF
integrator               = md
dt                       = $DT_NVT_PS
nsteps                   = $NSTEPS_NVT

refcoord_scaling         = com

; Center of mass removal
comm-mode                = Linear
nstcomm                  = 100
comm-grps                = System

; Output control
nstxout                  = 5000
nstvout                  = 0
nstfout                  = 0
nstlog                   = 5000
nstenergy                = 5000
nstxout-compressed       = 5000
compressed-x-precision   = 5000

; Neighborsearching
cutoff-scheme            = Verlet
nstlist                  = 20
rlist                    = 1.1

; Electrostatics
coulombtype              = reaction-field
rcoulomb                 = 1.1
epsilon_r                = 15
epsilon_rf               = 0

vdw_type                 = cutoff
vdw-modifier             = Potential-shift-verlet
rvdw                     = 1.1

pbc                      = xyz

; Temperature coupling
tcoupl                   = v-rescale
tc-grps                  = System
tau_t                    = 1.0 
ref_t                    = $REF_T

; NO PRESSURE COUPLING 

Pcoupl                   = no
Pcoupltype               = isotropic
tau_p                    = 20.0
compressibility          = 4e-5
ref_p                    =  1.0

; Velocity generation
gen_vel                  = yes
gen_temp                 = $REF_T
continuation             = no

constraints              = h-bonds
constraint_algorithm     = Lincs
EOF
    
# Create npt.mdp
NSTEPS_NPT=$(echo "$TIME_NPT_PS / $DT_NPT_PS" | bc -l | awk '{print int($1)}')
    
cat > npt.mdp << EOF
; NPT equilibration for Coarse-Grained MEMBRANE system

integrator               = md
dt                       = $DT_NPT_PS
nsteps                   = $NSTEPS_NPT


; Center of mass removal
comm-mode                = Linear
nstcomm                  = 100
comm-grps                = System
refcoord_scaling         = com
; Output control
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstlog                   = 5000
nstenergy                = 5000
nstxout-compressed       = 5000
compressed-x-precision   = 1000

; Neighborsearching
cutoff-scheme            = Verlet
nstlist                  = 20
rlist                    = 1.1

; Electrostatics
coulombtype              = reaction-field
rcoulomb                 = 1.1
epsilon_r                = 15
epsilon_rf               = 0

vdw_type                 = cutoff
vdw-modifier             = Potential-shift-verlet
rvdw                     = 1.1

pbc                      = xyz

; Temperature coupling (EXACTLY as your working NVT)
tcoupl                   = v-rescale
tc-grps                  = System
tau_t                    = 1.0 
ref_t                    = $REF_T

; ============= PRESSURE COUPLING ADDED FOR NPT =============

Pcoupl                   = Parrinello-Rahman
Pcoupltype               = isotropic
tau_p                    = 12.0               
compressibility          = 4.5e-5      
ref_p                    = $REF_P
; ============================================================

; Velocity generation
gen_vel                  = no               ; Continue from NVT
gen_temp                 = 310
continuation             = yes              ; Continuing from NVT

constraints              = h-bonds
lincs_order              = 8
lincs_iter               = 2
lincs_warnangle          = 90
EOF
    
# Create md.mdp
NSTEPS_MD=$(echo "$TIME_MD_PS / $DT_MD_PS" | bc -l | awk '{print int($1)}')
    
cat > md.mdp << EOF
; NPT equilibration for Coarse-Grained MD simulation

integrator               = md
dt                       = $DT_MD_PS
nsteps                   = $NSTEPS_MD

; Center of mass removal
comm-mode                = Linear
nstcomm                  = 100
comm-grps                = System

; Output control
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstlog                   = 1000
nstenergy                = 1000
nstxout-compressed       = 5000
compressed-x-precision   = 1000

; Neighborsearching
cutoff-scheme            = Verlet
nstlist                  = 20
rlist                    = 1.1

; Electrostatics
coulombtype              = reaction-field
rcoulomb                 = 1.1
epsilon_r                = 15
epsilon_rf               = 0

vdw_type                 = cutoff
vdw-modifier             = Potential-shift-verlet
rvdw                     = 1.1

pbc                      = xyz

; Temperature coupling (EXACTLY as your working NVT)
tcoupl                   = v-rescale
tc-grps                  = System
tau_t                    = 1.0 
ref_t                    = $REF_T

; ============= PRESSURE COUPLING ADDED FOR NPT =============

Pcoupl                   = Parrinello-Rahman
Pcoupltype               = isotropic
tau_p                    = 12.0               
compressibility          = 4.5e-5      
ref_p                    = $REF_P 
; ============================================================

; Velocity generation (CHANGED for continuation)
gen_vel                  = no               ; Continue from NVT
gen_temp                 = 310
continuation             = yes              ; Continuing from NVT

constraints              = h-bonds
lincs_order              = 8
lincs_iter               = 2
lincs_warnangle          = 90
EOF
    
echo "  ✓ MDP files created"
    
# ------------------------------------------------------------------------------
# STEP 2: RUN ENERGY MINIMIZATION
# ------------------------------------------------------------------------------
echo ""
echo "=== STEP 2: Energy Minimization ==="
    
export OMP_NUM_THREADS=$NTOMP
    
gmx grompp -f minimization.mdp -p "$TOPOL_CG_FILE" -c "$SOLV_IONS_GRO" -o em.tpr -r "$SOLV_IONS_GRO" -maxwarn 2
gmx mdrun -v -deffnm em -ntmpi $NTMPI -ntomp $NTOMP
    
echo "  ✓ Minimization completed"
    
# ------------------------------------------------------------------------------
# STEP 3: RUN NVT EQUILIBRATION
# ------------------------------------------------------------------------------
echo ""
echo "=== STEP 3: NVT Equilibration ==="
    
gmx grompp -f nvt.mdp -p "$TOPOL_CG_FILE" -c em.gro -o nvt.tpr -r em.gro -maxwarn 2
run_mdrun "nvt"
    
echo "  ✓ NVT completed"
    
# ------------------------------------------------------------------------------
# STEP 4: RUN NPT EQUILIBRATION
# ------------------------------------------------------------------------------
echo ""
echo "=== STEP 4: NPT Equilibration ==="
    
gmx grompp -f npt.mdp -p "$TOPOL_CG_FILE" -c nvt.gro -o npt.tpr -r nvt.gro -t nvt.cpt -maxwarn 2
run_mdrun "npt"
    
echo "  ✓ NPT completed"
    
# ------------------------------------------------------------------------------
# STEP 5: RUN PRODUCTION MD
# ------------------------------------------------------------------------------
echo ""
echo "=== STEP 5: Production MD ==="
    
gmx grompp -f md.mdp -p "$TOPOL_CG_FILE" -c npt.gro -o md.tpr -r npt.gro -t npt.cpt -maxwarn 2
run_mdrun "md"
    
echo "  ✓ Production MD completed"
    
# ------------------------------------------------------------------------------
# STEP 6: TRAJECTORY PROCESSING (Remove PBC, center molecule)
# ------------------------------------------------------------------------------
echo ""
echo "=== STEP 6: Trajectory Processing ==="
    
# Generate TPR for processing
gmx grompp -f minimization.mdp -c "$SOLV_IONS_GRO" -p "$TOPOL_CG_FILE" -o solv_ions.tpr -maxwarn 2

# Function to run trjconv with optional index
run_trjconv() {
    local cmd="$1"
    local options="$2"
    
    if [ "$INDEX_FILE" != "None" ] && [ -f "$INDEX_FILE" ]; then
        eval "$cmd -n $INDEX_FILE $options"
    else
        eval "$cmd $options"
    fi
}
    
# Remove PBC (make molecules whole)
echo "$GROUP_2" | run_trjconv "gmx trjconv -s solv_ions.tpr -f md.xtc -o md_whole.xtc" "-pbc whole"
    
# Remove jumps
echo "$GROUP_2" | run_trjconv "gmx trjconv -s solv_ions.tpr -f md_whole.xtc -o md_nojump.xtc" "-pbc nojump"
    
# Center molecule in box
echo "$GROUP_1" "$GROUP_2" | run_trjconv "gmx trjconv -s solv_ions.tpr -f md_nojump.xtc -o md_center.xtc" "-center -pbc mol -ur compact"
    
echo "  ✓ Trajectory processed: md_center.xtc"
    
# ------------------------------------------------------------------------------
# STEP 7: ANALYZE BONDS, ANGLES, DIHEDRALS (auto_analyze)
# ------------------------------------------------------------------------------
echo ""
echo "=== STEP 7: Running auto_analyze ==="
    
auto_mart3 auto_analyze \
        --bonds_ndx "$BONDS_NDX" \
        --angles_ndx "$ANGLES_NDX" \
        --dihedrals_ndx "$DIHEDRALS_NDX" \
        --xtc_file "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/MD_CG/md_center.xtc" \
        --tpr_file "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/MD_CG/solv_ions.tpr" \
        --output_all_files "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/RESULTS/BONDS_ANGLES_DIHEDRALS_XVG_SIM"
    
echo "  ✓ Analysis completed"
    
# ------------------------------------------------------------------------------
# STEP 8: PLOT DISTRIBUTIONS (auto_plot_distributions)
# ------------------------------------------------------------------------------
echo ""
echo "=== STEP 8: Running auto_plot_distributions ==="
    
auto_mart3 auto_plot_distributions \
        --bonds_ref_dir "$BONDS_REF_XVG_DIR" \
        --angles_ref_dir "$ANGLES_REF_XVG_DIR" \
        --dihedrals_ref_dir "$DIHEDRALS_REF_XVG_DIR" \
        --bonds_sim_dir "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/RESULTS/BONDS_ANGLES_DIHEDRALS_XVG_SIM/bonds" \
        --angles_sim_dir "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/RESULTS/BONDS_ANGLES_DIHEDRALS_XVG_SIM/angles" \
        --dihedrals_sim_dir "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/RESULTS/BONDS_ANGLES_DIHEDRALS_XVG_SIM/dihedrals" \
        --figures_dir "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/RESULTS/FIGURES_DISTR_SERIES_REF_VS_SIM"
    
echo "  ✓ Plots saved to $OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/RESULTS/FIGURES_DISTR_SERIES_REF_VS_SIM"
    
# ------------------------------------------------------------------------------
# STEP 9: UPDATE FORCE CONSTANTS (bayes_potential_adjust)
# ------------------------------------------------------------------------------
echo ""
echo "=== STEP 9: Updating force constants ==="
    
auto_mart3 bayes_potential_adjust \
        --bonds_ref_xvg_dir "$BONDS_REF_XVG_DIR" \
        --angles_ref_xvg_dir "$ANGLES_REF_XVG_DIR" \
        --dihedrals_ref_xvg_dir "$DIHEDRALS_REF_XVG_DIR" \
        --bonds_sim_xvg_dir "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/RESULTS/BONDS_ANGLES_DIHEDRALS_XVG_SIM/bonds" \
        --angles_sim_xvg_dir "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/RESULTS/BONDS_ANGLES_DIHEDRALS_XVG_SIM/angles" \
        --dihedrals_sim_xvg_dir "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/RESULTS/BONDS_ANGLES_DIHEDRALS_XVG_SIM/dihedrals" \
        --itp_cg "$ITP_TO_OPTIMIZE" \
        --ndx_bounds "$BONDS_NDX" \
        --ndx_angles "$ANGLES_NDX" \
        --ndx_dihedrals "$DIHEDRALS_NDX" \
        --molecule_name "$MOLECULE_NAME" \
        --T0 $T0 \
        --alpha $ALPHA \
        --itp_out "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/MD_CG/cg_updated_iter_${i}.itp" \
        --distribution_points $DISTRIBUTION_POINTS

# Remove GMX files unless --no_remove_gmx_files_in_iter was specified
if [ "$NO_REMOVE_GMX_FILES" = false ]; then
    rm -fr *.xtc *.edr *.trr *.log *.cpt \#*
    echo "  ✓ GMX temporary files removed"
else
    echo "  GMX temporary files preserved (--no_remove_gmx_files_in_iter specified)"
fi
    
echo "  ✓ Force constants updated in iteration ${i}"
echo " Topology was updated in $OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/cg_updated_iter_${i}.itp"
    
# ------------------------------------------------------------------------------
# STEP 10: STARTING OPTIMIZATION LOOP
# ------------------------------------------------------------------------------
echo ""
echo "=== STEP 10: Starting optimization loop ==="

for i in $(seq 1 $N_ITER); do
    
    # ------------------------------------------------------------------------------
    # STEP 11: Loop of optimization
    # ------------------------------------------------------------------------------
    
    mkdir -p "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/MD_CG"
    mkdir -p "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/RESULTS/BONDS_ANGLES_DIHEDRALS_XVG_SIM"/{bonds,angles,dihedrals}
    mkdir -p "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/RESULTS/FIGURES_DISTR_SERIES_REF_VS_SIM"
    
    
    # Copy previous iteration CG system files
    if [ -d "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_$((i-1))/MD_CG" ]; then
        cp -r "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_$((i-1))/MD_CG"/* "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/MD_CG" 2>/dev/null || true
    fi
    
    cd "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/MD_CG"
    
    echo "Working in: $(pwd)"
    
    # Copy the updated ITP file from previous iteration
    cp "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_$((i-1))/MD_CG/cg_updated_iter_$((i-1)).itp" "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/MD_CG/cg.itp"
    ITP_TO_OPTIMIZE="cg.itp"
    
    # ------------------------------------------------------------------------------
    # STEP 12: PREPARING FILES FOR SIMULATION
    # ------------------------------------------------------------------------------
    echo ""
    echo "=== Files transferred to $OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/MD_CG ==="
    echo ""
    echo "    $ITP_TO_OPTIMIZE was optimized according to iteration $((i-1))"
    
    # ------------------------------------------------------------------------------
    # STEP 2: RUN ENERGY MINIMIZATION
    # ------------------------------------------------------------------------------
    echo ""
    echo "=== STEP 2: Energy Minimization ==="
    
    export OMP_NUM_THREADS=$NTOMP
    
    gmx grompp -f minimization.mdp -p "$TOPOL_CG_FILE" -c md.gro -o em.tpr -r md.gro -maxwarn 2
    gmx mdrun -v -deffnm em -ntmpi $NTMPI -ntomp $NTOMP
    
    echo "  ✓ Minimization completed"
    
    # ------------------------------------------------------------------------------
    # STEP 3: RUN NVT EQUILIBRATION
    # ------------------------------------------------------------------------------
    echo ""
    echo "=== STEP 3: NVT Equilibration ==="
    
    gmx grompp -f nvt.mdp -p "$TOPOL_CG_FILE" -c em.gro -o nvt.tpr -r em.gro -maxwarn 2
    run_mdrun "nvt"
    
    echo "  ✓ NVT completed"
    
    # ------------------------------------------------------------------------------
    # STEP 4: RUN NPT EQUILIBRATION
    # ------------------------------------------------------------------------------
    echo ""
    echo "=== STEP 4: NPT Equilibration ==="
    
    gmx grompp -f npt.mdp -p "$TOPOL_CG_FILE" -c nvt.gro -o npt.tpr -r nvt.gro -t nvt.cpt -maxwarn 2
    run_mdrun "npt"
    
    echo "  ✓ NPT completed"
    
    # ------------------------------------------------------------------------------
    # STEP 5: RUN PRODUCTION MD
    # ------------------------------------------------------------------------------
    echo ""
    echo "=== STEP 5: Production MD ==="
    
    gmx grompp -f md.mdp -p "$TOPOL_CG_FILE" -c npt.gro -o md.tpr -r npt.gro -t npt.cpt -maxwarn 2
    run_mdrun "md"
    
    echo "  ✓ Production MD completed"
    
    # ------------------------------------------------------------------------------
    # STEP 6: TRAJECTORY PROCESSING (Remove PBC, center molecule)
    # ------------------------------------------------------------------------------
    echo ""
    echo "=== STEP 6: Trajectory Processing ==="
    
    # Generate TPR for processing
    gmx grompp -f minimization.mdp -c "$SOLV_IONS_GRO" -p "$TOPOL_CG_FILE" -o solv_ions.tpr -maxwarn 2
    
    # Remove PBC (make molecules whole)
    echo "$GROUP_2" | run_trjconv "gmx trjconv -s solv_ions.tpr -f md.xtc -o md_whole.xtc" "-pbc whole"
    
    # Remove jumps
    echo "$GROUP_2" | run_trjconv "gmx trjconv -s solv_ions.tpr -f md_whole.xtc -o md_nojump.xtc" "-pbc nojump"
    
    # Center molecule in box
    echo "$GROUP_1" "$GROUP_2" | run_trjconv "gmx trjconv -s solv_ions.tpr -f md_nojump.xtc -o md_center.xtc" "-center -pbc mol -ur compact"
    
    echo "  ✓ Trajectory processed: md_center.xtc"
    
    # ------------------------------------------------------------------------------
    # STEP 7: ANALYZE BONDS, ANGLES, DIHEDRALS (auto_analyze)
    # ------------------------------------------------------------------------------
    echo ""
    echo "=== STEP 7: Running auto_analyze ==="
    
    auto_mart3 auto_analyze \
        --bonds_ndx "$BONDS_NDX" \
        --angles_ndx "$ANGLES_NDX" \
        --dihedrals_ndx "$DIHEDRALS_NDX" \
        --xtc_file "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/MD_CG/md_center.xtc" \
        --tpr_file "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/MD_CG/solv_ions.tpr" \
        --output_all_files "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/RESULTS/BONDS_ANGLES_DIHEDRALS_XVG_SIM"
    
    echo "  ✓ Analysis completed"
    
    # ------------------------------------------------------------------------------
    # STEP 8: PLOT DISTRIBUTIONS (auto_plot_distributions)
    # ------------------------------------------------------------------------------
    echo ""
    echo "=== STEP 8: Running auto_plot_distributions ==="
    
    auto_mart3 auto_plot_distributions \
        --bonds_ref_dir "$BONDS_REF_XVG_DIR" \
        --angles_ref_dir "$ANGLES_REF_XVG_DIR" \
        --dihedrals_ref_dir "$DIHEDRALS_REF_XVG_DIR" \
        --bonds_sim_dir "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/RESULTS/BONDS_ANGLES_DIHEDRALS_XVG_SIM/bonds" \
        --angles_sim_dir "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/RESULTS/BONDS_ANGLES_DIHEDRALS_XVG_SIM/angles" \
        --dihedrals_sim_dir "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/RESULTS/BONDS_ANGLES_DIHEDRALS_XVG_SIM/dihedrals" \
        --figures_dir "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/RESULTS/FIGURES_DISTR_SERIES_REF_VS_SIM"
    
    echo "  ✓ Plots saved to $OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/RESULTS/FIGURES_DISTR_SERIES_REF_VS_SIM"
    
    # ------------------------------------------------------------------------------
    # STEP 9: UPDATE FORCE CONSTANTS (bayes_potential_adjust)
    # ------------------------------------------------------------------------------
    echo ""
    echo "=== STEP 9: Updating force constants ==="
    
    auto_mart3 bayes_potential_adjust \
        --bonds_ref_xvg_dir "$BONDS_REF_XVG_DIR" \
        --angles_ref_xvg_dir "$ANGLES_REF_XVG_DIR" \
        --dihedrals_ref_xvg_dir "$DIHEDRALS_REF_XVG_DIR" \
        --bonds_sim_xvg_dir "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/RESULTS/BONDS_ANGLES_DIHEDRALS_XVG_SIM/bonds" \
        --angles_sim_xvg_dir "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/RESULTS/BONDS_ANGLES_DIHEDRALS_XVG_SIM/angles" \
        --dihedrals_sim_xvg_dir "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/RESULTS/BONDS_ANGLES_DIHEDRALS_XVG_SIM/dihedrals" \
        --itp_cg "$ITP_TO_OPTIMIZE" \
        --ndx_bounds "$BONDS_NDX" \
        --ndx_angles "$ANGLES_NDX" \
        --ndx_dihedrals "$DIHEDRALS_NDX" \
        --molecule_name "$MOLECULE_NAME" \
        --T0 $T0 \
        --alpha $ALPHA \
        --itp_out "$OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/MD_CG/cg_updated_iter_${i}.itp" \
        --distribution_points $DISTRIBUTION_POINTS

    ITP_TO_OPTIMIZE="cg_updated_iter_${i}.itp"
    
    # Remove GMX files unless --no_remove_gmx_files_in_iter was specified
    if [ "$NO_REMOVE_GMX_FILES" = false ]; then
        rm -fr *.xtc *.edr *.trr *.log *.cpt \#*
        echo "  ✓ GMX temporary files removed"
    else
        echo "  GMX temporary files preserved (--no_remove_gmx_files_in_iter specified)"
    fi
    
    echo "  ✓ Force constants updated in iteration ${i}"
    echo " Topology was updated in $OUTPUT_AUTO_MART_CG_DIR/CALIBRATION/iter_${i}/MD_CG/cg_updated_iter_${i}.itp"
    
done

echo ""
echo "=== Optimization completed ==="
