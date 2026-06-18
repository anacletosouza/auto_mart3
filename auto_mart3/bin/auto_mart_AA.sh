#!/bin/bash

# ================================================================
# AUTO_MART_AA: COMPLETE CG MAPPING AND PARAMETERIZATION PIPELINE
# ================================================================
# This pipeline performs:
#   1. Atomistic to Coarse-Grained (CG) mapping using cg-martini3
#   2. Trajectory mapping from AA to CG representation
#   3. Topology adaptation and preparation
#   4. Bond, angle, and dihedral analysis
#   5. Statistical distribution calculations
#   6. Final system preparation for CG simulations
#
# REQUIRED ARGUMENTS (must be provided by user):
#   --aa_tpr       : Path to AA simulation .tpr file (run input file)
#   --aa_xtc       : Path to AA trajectory .xtc file (compressed trajectory)
#   --aa_gro       : Path to AA structure .gro file (coordinates)
#   --aa_itp       : Path to AA filtered .itp file (topology for molecule)
#   --beads_json   : Path to JSON config file defining CG beads mapping
#   --input_mdp    : Path to GROMACS .mdp file (minimization/equilibration params)
#   --path_ff      : Path to force field directory (contains .itp files)
#   --output_dir   : Path to output directory for all results
#   --name_molecule: Name of the molecule (e.g., "protein", "ligand")
#
# OPTIONAL ARGUMENTS (default values provided):
#   --force_application : Force application method (default: "random=[1250,30;30,1]")
#   --beads_position    : Beads positioning method (default: "geom" = geometric center)
#   --cycle_restr       : Cycle restraints configuration (default: "fix=3,mode=cycle")
#   --maxwarn           : Max warnings for gmx grompp (default: 2)
#   --distance_from_atom: Distance for solvent placement (default: 2.0 nm)
#   --salt              : Salt concentration in M (default: 0.15 M)
# =====================================================

# Author: Anacleto Silva de Souza

# =====================================================
# Helper function to convert relative path to absolute path
# =====================================================
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
        if [[ -f "$path" ]] || [[ -d "$path" ]]; then
            echo "$(cd "$(dirname "$path")" 2>/dev/null && pwd)/$(basename "$path")" 2>/dev/null || echo "$(pwd)/$path"
        else
            # If file/directory doesn't exist yet, convert based on current directory
            echo "$(pwd)/$path"
        fi
    fi
}

# =====================================================
# PARSE COMMAND LINE ARGUMENTS (mantém hífens nos argumentos)
# =====================================================
while [[ $# -gt 0 ]]; do
    case $1 in
        --aa_tpr)
            AA_TPR="$2"
            shift 2
            ;;
        --aa_xtc)
            AA_XTC="$2"
            shift 2
            ;;
        --aa_gro)
            AA_GRO="$2"
            shift 2
            ;;
        --aa_itp)
            AA_ITP="$2"
            shift 2
            ;;
        --beads_json)
            BEADS_JSON="$2"
            shift 2
            ;;
        --input_mdp)
            INPUT_MDP="$2"
            shift 2
            ;;
        --path_ff)
            PATH_FF="$2"
            shift 2
            ;;
        --output_dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --name_molecule)
            NAME_MOLECULE="$2"
            shift 2
            ;;
        --force_application)
            FORCE_APPLICATION="$2"
            shift 2
            ;;
        --beads_position)
            BEADS_POSITION="$2"
            shift 2
            ;;
        --cycle_restr)
            CYCLE_RESTR="$2"
            shift 2
            ;;
        --maxwarn)
            MAXWARN="$2"
            shift 2
            ;;
        --distance_from_atom)
            DISTANCE_FROM_ATOM="$2"
            shift 2
            ;;
        --salt)
            SALT="$2"
            shift 2
            ;;
        *)
            echo "ERROR: Unknown option $1"
            echo "Usage: $0 --aa_tpr <path> --aa_xtc <path> --aa_gro <path> --aa_itp <path> --beads_json <path> --input_mdp <path> --path_ff <path> --output_dir <path> --name_molecule <name> [optional args]"
            exit 1
            ;;
    esac
done

# =====================================================
# VALIDATE REQUIRED ARGUMENTS
# =====================================================
if [ -z "$AA_TPR" ] || [ -z "$AA_XTC" ] || [ -z "$AA_GRO" ] || [ -z "$AA_ITP" ] || \
   [ -z "$BEADS_JSON" ] || [ -z "$INPUT_MDP" ] || [ -z "$PATH_FF" ] || \
   [ -z "$OUTPUT_DIR" ] || [ -z "$NAME_MOLECULE" ]; then
    echo "ERROR: Missing required arguments!"
    echo "Required arguments:"
    echo "  --aa_tpr, --aa_xtc, --aa_gro, --aa_itp, --beads_json,"
    echo "  --input_mdp, --path_ff, --output_dir, --name_molecule"
    exit 1
fi

# =====================================================
# CONVERT ALL PATHS TO ABSOLUTE PATHS
# =====================================================
echo "=========================================="
echo "CONVERTING PATHS TO ABSOLUTE"
echo "=========================================="

AA_TPR=$(to_absolute_path "$AA_TPR")
AA_XTC=$(to_absolute_path "$AA_XTC")
AA_GRO=$(to_absolute_path "$AA_GRO")
AA_ITP=$(to_absolute_path "$AA_ITP")
BEADS_JSON=$(to_absolute_path "$BEADS_JSON")
INPUT_MDP=$(to_absolute_path "$INPUT_MDP")
PATH_FF=$(to_absolute_path "$PATH_FF")
OUTPUT_DIR=$(to_absolute_path "$OUTPUT_DIR")

echo "AA_TPR (absolute): $AA_TPR"
echo "AA_XTC (absolute): $AA_XTC"
echo "AA_GRO (absolute): $AA_GRO"
echo "AA_ITP (absolute): $AA_ITP"
echo "BEADS_JSON (absolute): $BEADS_JSON"
echo "INPUT_MDP (absolute): $INPUT_MDP"
echo "PATH_FF (absolute): $PATH_FF"
echo "OUTPUT_DIR (absolute): $OUTPUT_DIR"
echo ""

# =====================================================
# SET DEFAULT VALUES FOR OPTIONAL ARGUMENTS
# =====================================================
FORCE_APPLICATION=${FORCE_APPLICATION:-"random=[1250,30;30,1]"}
BEADS_POSITION=${BEADS_POSITION:-"geom"}
CYCLE_RESTR=${CYCLE_RESTR:-"fix=3,mode=cycle"}
MAXWARN=${MAXWARN:-2}
DISTANCE_FROM_ATOM=${DISTANCE_FROM_ATOM:-2.0}
SALT=${SALT:-0.15}

# =====================================================
# CREATE OUTPUT DIRECTORY
# =====================================================
echo "=========================================="
echo "CG MAPPING AND PARAMETERIZATION PIPELINE"
echo "=========================================="
echo "Output directory: $OUTPUT_DIR"
echo "Molecule name: $NAME_MOLECULE"
echo ""

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || exit 1

# =====================================================
# STEP 0: cg-martini3 (AA -> CG mapping)
# =====================================================
echo "=== STEP 0: Running cg-martini3 (AA -> CG mapping) ==="
echo "   Mapping atomistic structure to coarse-grained beads..."
cg-martini3 \
    --input_gro "$AA_GRO" \
    --input_beads_definitions "$BEADS_JSON" \
    --output_dir CG_MARTINI3 \
    --beads_position "$BEADS_POSITION" \
    --aa_itp "$AA_ITP" \
    --force_application "$FORCE_APPLICATION" \
    --cycle_restr "$CYCLE_RESTR" \
    --name_molecule "$NAME_MOLECULE"

if [ $? -ne 0 ]; then
    echo "ERROR: cg-martini3 failed!"
    exit 1
fi

mkdir -p "$OUTPUT_DIR/GMX"
echo "   ✓ Step 0 completed"
echo ""

# =====================================================
# STEP 1: auto_map (trajectory mapping) - usando UNDERSCORES
# =====================================================
echo "=== STEP 1: Running auto_map (trajectory mapping) ==="
echo "   Mapping atomistic trajectory to CG representation..."
auto_mart3 auto_map \
    --index_cg "CG_MARTINI3/NDX/cg.ndx" \
    --aa_tpr "$AA_TPR" \
    --aa_xtc "$AA_XTC" \
    --output_mapped "GMX/mapped.xtc" \
    --output_cg_gro "GMX/cg.gro" \
    --remove_pbc

if [ $? -ne 0 ]; then
    echo "ERROR: auto_map failed!"
    exit 1
fi
echo "   ✓ Step 1 completed"
echo ""

# =====================================================
# STEP 2: Copy generated ITP
# =====================================================
echo "=== STEP 2: Copying CG ITP ==="
cp "CG_MARTINI3/ITP/final_cg.itp" "GMX/cg.itp"
echo "   ✓ Step 2 completed"
echo ""

# =====================================================
# STEP 2.5: auto_adapt_itp (adapt atom names) - usando UNDERSCORES
# =====================================================
echo "=== STEP 2.5: Running auto_adapt_itp (adapt atom names) ==="
echo "   Adapting atom names in topology to match structure file..."
auto_mart3 auto_adapt_itp \
    --input_itp "GMX/cg.itp" \
    --input_gro_ref "GMX/cg.gro" \
    --output_itp_adapted "GMX/cg.itp.adapted"

if [ $? -ne 0 ]; then
    echo "ERROR: auto_adapt_itp failed!"
    exit 1
fi

mv "GMX/cg.itp.adapted" "GMX/cg.itp"
echo "   ✓ Step 2.5 completed"
echo ""

# =====================================================
# STEP 3: auto_gen_top + grompp - usando UNDERSCORES
# =====================================================
echo "=== STEP 3: Generating topology and running grompp ==="
echo "   Creating complete system topology..."

# Create temporary directory
mkdir -p topology_temp
cd topology_temp || exit 1

# Copy required files
cp "$PATH_FF"/*.itp .
cp "$OUTPUT_DIR/GMX/cg.itp" .
cp "$OUTPUT_DIR/GMX/cg.gro" .

# Generate topology
auto_mart3 auto_gen_top \
    --path_ff . \
    --ff martini_v3.0.0.itp \
    --ions martini_v3.0.0_ions_v1.itp \
    --solvent martini_v3.0.0_solvents_v1.itp \
    --itp_ligand cg.itp \
    --name_molecule "$NAME_MOLECULE" \
    --number_molecule 1 \
    --output_topol topol_cg.top

if [ $? -ne 0 ]; then
    echo "ERROR: auto_gen_top failed!"
    exit 1
fi

# Fix includes (remove absolute paths)
sed -i 's|#include ".*/|#include "|g' topol_cg.top

# Run grompp
echo "   Running grompp to create TPR file..."
gmx grompp \
    -f "$INPUT_MDP" \
    -c cg.gro \
    -p topol_cg.top \
    -o CG.tpr \
    -maxwarn "$MAXWARN"

if [ $? -ne 0 ]; then
    echo "ERROR: gmx grompp failed!"
    exit 1
fi

# Copy results back
cp CG.tpr "$OUTPUT_DIR/GMX/"
cp topol_cg.top "$OUTPUT_DIR/GMX/"

cd "$OUTPUT_DIR" || exit 1
rm -rf topology_temp
echo "   ✓ Step 3 completed"
echo ""

# =====================================================
# STEP 4: auto_analyze (bonds, angles, dihedrals analysis) - usando UNDERSCORES
# =====================================================
echo "=== STEP 4: Running auto_analyze (bonds, angles, dihedrals) ==="
echo "   Analyzing CG trajectory for internal coordinates..."

mkdir -p "$OUTPUT_DIR/BONDS_ANGLES_DIHEDRALS_XVG_REF"

auto_mart3 auto_analyze \
    --bonds_ndx "$OUTPUT_DIR/CG_MARTINI3/NDX/bonds.ndx" \
    --angles_ndx "$OUTPUT_DIR/CG_MARTINI3/NDX/angles.ndx" \
    --dihedrals_ndx "$OUTPUT_DIR/CG_MARTINI3/NDX/dihedrals.ndx" \
    --xtc_file "$OUTPUT_DIR/GMX/mapped.xtc" \
    --tpr_file "$OUTPUT_DIR/GMX/CG.tpr" \
    --output_all_files "$OUTPUT_DIR/BONDS_ANGLES_DIHEDRALS_XVG_REF"

if [ $? -ne 0 ]; then
    echo "ERROR: auto_analyze failed!"
    exit 1
fi
echo "   ✓ Step 4 completed"
echo ""

# =====================================================
# STEP 5: auto_distributions (distribution statistics) - usando UNDERSCORES
# =====================================================
echo "=== STEP 5: Running auto_distributions (statistics) ==="
echo "   Calculating statistical distributions..."

auto_mart3 auto_distributions \
    --bonds_dir "$OUTPUT_DIR/BONDS_ANGLES_DIHEDRALS_XVG_REF/bonds" \
    --angles_dir "$OUTPUT_DIR/BONDS_ANGLES_DIHEDRALS_XVG_REF/angles" \
    --dihedrals_dir "$OUTPUT_DIR/BONDS_ANGLES_DIHEDRALS_XVG_REF/dihedrals" \
    --dir_to_output "$OUTPUT_DIR/STATISTICS" \
    --bond_out bond_statistics.tsv \
    --angle_out angle_statistics.tsv \
    --dihedral_out dihedral_statistics.tsv

if [ $? -ne 0 ]; then
    echo "ERROR: auto_distributions failed!"
    exit 1
fi
echo "   ✓ Step 5 completed"
echo ""

# =====================================================
# STEP 6: auto_prep (prepare CG system for simulation) - usando UNDERSCORES
# =====================================================
echo "=== STEP 6: Running auto_prep (prepare CG system) ==="
echo "   Preparing solvated CG system for simulation..."

auto_mart3 auto_prep \
    --input_ref_dir "$OUTPUT_DIR/GMX" \
    --input_gro cg.gro \
    --input_itp cg.itp \
    --input_topol topol_cg.top \
    --output_dir "$OUTPUT_DIR/MDRUN_CG" \
    --ff martini_v3.0.0.itp \
    --ions martini_v3.0.0_ions_v1.itp \
    --solvent martini_v3.0.0_solvents_v1.itp \
    --input_ff_dir "$PATH_FF" \
    --use_distance_from_atom \
    --distance_from_atom "$DISTANCE_FROM_ATOM" \
    --salt "$SALT"

if [ $? -ne 0 ]; then
    echo "ERROR: auto_prep failed!"
    exit 1
fi
echo "   ✓ Step 6 completed"
echo ""

# =====================================================
# SUMMARY
# =====================================================
echo "=========================================="
echo "PIPELINE COMPLETED SUCCESSFULLY"
echo "=========================================="
echo "Output directory: $OUTPUT_DIR"
echo "Molecule: $NAME_MOLECULE"
echo ""
echo "=== GENERATED FILES AND DIRECTORIES ==="
echo ""
echo "   CG_MARTINI3/           # Initial CG mapping results"
echo "   ├── NDX/"
echo "   │   ├── cg.ndx         # CG atom indices"
echo "   │   ├── bonds.ndx      # Bond connectivity indices"
echo "   │   ├── angles.ndx     # Angle indices"
echo "   │   └── dihedrals.ndx  # Dihedral indices"
echo "   ├── ITP/"
echo "   │   └── final_cg.itp   # Original CG topology"
echo "   ├── GRO/"
echo "   │   └── cg.gro         # Original CG structure"
echo "   └── PDB/"
echo "       └── cg.pdb         # CG structure in PDB format"
echo ""
echo "   GMX/                   # GROMACS-compatible files"
echo "   ├── mapped.xtc         # CG trajectory (from AA mapping)"
echo "   ├── cg.gro             # CG structure (frame 0)"
echo "   ├── cg.itp             # Adapted CG topology (names fixed)"
echo "   ├── CG.tpr             # GROMACS run input file"
echo "   └── topol_cg.top       # Complete system topology (no solvent)"
echo ""
echo "   BONDS_ANGLES_DIHEDRALS_XVG_REF/  # Analysis output"
echo "   ├── bonds/             # Bond distance vs time (XVG files)"
echo "   ├── angles/            # Angle values vs time (XVG files)"
echo "   ├── dihedrals/         # Dihedral angles vs time (XVG files)"
echo "   ├── bonds_histograms/  # Bond distance histograms"
echo "   ├── angles_histograms/ # Angle histograms"
echo "   └── dihedrals_histograms/ # Dihedral histograms"
echo ""
echo "   STATISTICS/            # Distribution statistics"
echo "   ├── bond_statistics.tsv     # Bond: mean, std, bins, counts"
echo "   ├── angle_statistics.tsv    # Angle: mean, std, bins, counts"
echo "   └── dihedral_statistics.tsv # Dihedral: mean, std, bins, counts"
echo ""
echo "   MDRUN_CG/              # Simulation-ready system"
echo "   ├── topol.top          # Final topology (with solvent and ions)"
echo "   ├── conf.gro           # Solvated CG structure"
echo "   ├── em.mdp             # Energy minimization parameters"
echo "   ├── em.tpr             # Energy minimization TPR"
echo "   ├── em.gro             # Energy-minimized structure"
echo "   ├── md.mdp             # MD simulation parameters"
echo "   ├── md.tpr             # MD simulation TPR"
echo "   └── run_simulation.sh  # Script to run the simulation"
echo ""
echo "=========================================="
echo "NEXT STEPS:"
echo "1. Review the generated statistics in STATISTICS/"
echo "2. Check the solvated system in MDRUN_CG/conf.gro"
echo "3. Run energy minimization: cd MDRUN_CG && gmx mdrun -deffnm em"
echo "4. Run production MD: cd MDRUN_CG && bash run_simulation.sh"
echo "=========================================="
