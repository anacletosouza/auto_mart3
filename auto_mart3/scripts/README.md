# Auto_Mart3: Automated Coarse-Grained Mapping and Parameterization Pipeline for Martini 3

Auto_Mart3 is a comprehensive Python pipeline for mapping atomistic molecular dynamics trajectories to coarse-grained (CG) resolution, generating CG topologies, analyzing structural properties (bonds, angles, dihedrals), and preparing full simulation systems using the Martini 3 force field.

## Table of Contents

- [Overview](#overview)
- [Pipeline Workflow](#pipeline-workflow)
- [Module Descriptions](#module-descriptions)
  - [1. map_aa_traj_to_cg.py](#1-map_aa_traj_to_cgpy)
  - [2. generate_cg_top.py](#2-generate_cg_toppy)
  - [3. generate_bonds_angles_dihedrals.py](#3-generate_bonds_angles_dihedralspy)
  - [4. bp_distributions.py](#4-bp_distributionspy)
  - [5. adaptation_gro_itp.py](#5-adaptation_gro_itppy)
  - [6. bp_prep.py](#6-bp_preppy)
  - [7. plot_distributions.py](#7-plot_distributionspy)
- [Quick Start](#quick-start)
- [Dependencies](#dependencies)

---

## Overview

Auto_Mart3 automates the conversion of atomistic molecular dynamics simulations to coarse-grained representations compatible with the Martini 3 force field. The pipeline handles:

- **Mapping** atomistic trajectories to CG beads using GROMACS
- **Topology generation** with proper force field includes
- **Structural analysis** of bonds, angles, and dihedrals
- **Statistical processing** of distribution data
- **System preparation** for CG simulations (box definition, solvation, ions)
- **Visualization** of reference vs simulated distributions

---

## Pipeline Workflow

```
┌─────────────────┐
│   AA Trajectory │
│   (md.xtc/tpr)  │
└────────┬────────┘
         │
         ▼
┌─────────────────────────────┐
│ 1. map_aa_traj_to_cg.py     │
│    Map AA → CG beads         │
└────────┬────────────────────┘
         │
         ▼
┌─────────────────────────────┐
│ 2. generate_cg_top.py       │
│    Create CG topology        │
└────────┬────────────────────┘
         │
         ▼
┌─────────────────────────────┐
│ 3. generate_bonds_angles_   │
│    dihedrals.py              │
│    Extract structural data   │
└────────┬────────────────────┘
         │
         ▼
┌─────────────────────────────┐
│ 4. bp_distributions.py      │
│    Compute statistics        │
└────────┬────────────────────┘
         │
         ▼
┌─────────────────────────────┐
│ 5. adaptation_gro_itp.py    │
│    Sync atom names           │
└────────┬────────────────────┘
         │
         ▼
┌─────────────────────────────┐
│ 6. bp_prep.py               │
│    Prepare simulation system│
└────────┬────────────────────┘
         │
         ▼
┌─────────────────────────────┐
│ 7. plot_distributions.py    │
│    Visualize results         │
└─────────────────────────────┘
```

---

## Module Descriptions

### 1. map_aa_traj_to_cg.py

**Description:** Maps atomistic (AA) trajectories to coarse-grained (CG) resolution using GROMACS. Automatically generates CG `.gro` and `.xtc` files from AA MD simulations. The script handles PBC correction and computes centers of mass for each CG bead defined in the index file.

**Required files:**
- `md.xtc` - AA trajectory (with or without PBC)
- `cg.ndx` - CG bead definitions (atoms grouped per bead)
- `md.tpr` - AA topology/run input file

**Key features:**
- Automatic PBC removal with `gmx trjconv -pbc whole`
- CG bead mapping using center-of-mass calculation
- Automatic CG `.gro` extraction from first frame
- Flexible output naming

**Examples:**

```bash
# Basic mapping with PBC removal (default)
python map_aa_traj_to_cg.py \
    --index_cg cg.ndx \
    --aa_tpr md.tpr \
    --aa_xtc md.xtc \
    --output_mapped mapped.xtc \
    --output_cg_gro molecule.gro

# Using trajectory with PBC already corrected
python map_aa_traj_to_cg.py \
    --index_cg cg.ndx \
    --aa_tpr md.tpr \
    --aa_xtc md_no_pbc.xtc \
    --output_mapped mapped.xtc \
    --corrected_pbc \
    --output_cg_gro molecule.gro

# Dry run to check commands without executing
python map_aa_traj_to_cg.py \
    --index_cg cg.ndx \
    --aa_tpr md.tpr \
    --aa_xtc md.xtc \
    --dry-run \
    --verbose

# Keep temporary files for debugging
python map_aa_traj_to_cg.py \
    --index_cg cg.ndx \
    --aa_tpr md.tpr \
    --aa_xtc md.xtc \
    --keep-temp
```

**Arguments:**

| Argument | Description | Default |
|----------|-------------|---------|
| `--index_cg`, `-i` | CG index file (.ndx) | **Required** |
| `--aa_tpr`, `-t` | Atomistic .tpr file | **Required** |
| `--aa_xtc`, `-x` | Atomistic trajectory (.xtc) | **Required** |
| `--output_mapped`, `-o` | Output CG trajectory | mapped.xtc |
| `--output_cg_gro` | Output CG .gro file (first frame) | molecule.gro |
| `--remove_pbc` | Remove PBC with gmx trjconv | True |
| `--corrected_pbc` | Skip PBC correction (use when already corrected) | False |
| `--verbose`, `-v` | Verbose output | False |
| `--dry-run` | Print commands only | False |
| `--keep-temp` | Keep temporary files | False |
| `--gmx_cmd` | GROMACS command | gmx |

---

### 2. generate_cg_top.py

**Description:** Generates a GROMACS topology (`.top`) file for CG simulations, assembling necessary `#include` directives for force field, ions, solvent, and ligand definitions. The output is ready for use with GROMACS `grompp`.

**Key features:**
- Flexible force field file naming
- Automatic path handling for ITP files
- Optional system title and comments
- Support for ligand ITP inclusion

**Examples:**

```bash
# Basic topology with all Martini 3 components
python generate_cg_top.py \
    --path_ff ../ff_files/ \
    --ff martini_v3.0.0.itp \
    --ions martini_v3.0.0_ions_v1.itp \
    --solvent martini_v3.0.0_solvents_v1.itp \
    --itp_ligand ../setup/cg.itp \
    --name_molecule MOL \
    --number_molecule 1 \
    --output_topol topol_cg.top

# Skip ligand (e.g., for pure solvent system)
python generate_cg_top.py \
    --path_ff ff_files/ \
    --itp_ligand None \
    --name_molecule SOL \
    --number_molecule 1000 \
    --output_topol topol_cg.top \
    --title_system "Pure Martini water box" \
    --title_comments "Generated by Auto_Mart3"

# Multiple molecules (e.g., 2 ligands in solution)
python generate_cg_top.py \
    --path_ff ff_files/ \
    --itp_ligand molecule.itp \
    --name_molecule LIG \
    --number_molecule 2 \
    --output_topol topol_cg.top \
    --title_system "Ligand in aqueous solution"
```

**Arguments:**

| Argument | Description | Default |
|----------|-------------|---------|
| `--path_ff` | Directory containing force field ITP files | **Required** |
| `--ff` | Force field ITP filename | martini_v3.0.0.itp |
| `--ions` | Ions ITP filename | martini_v3.0.0_ions_v1.itp |
| `--solvent` | Solvent ITP filename | martini_v3.0.0_solvents_v1.itp |
| `--itp_ligand` | Ligand ITP file (use "None" to skip) | **Required** |
| `--name_molecule` | Molecule name for topology | **Required** |
| `--number_molecule` | Number of molecules | **Required** |
| `--title_comments` | Optional comment header | "" |
| `--title_system` | Optional system name | "" |
| `--output_topol` | Output topology filename | **Required** |

---

### 3. generate_bonds_angles_dihedrals.py

**Description:** Computes bond distances, angles, and dihedral distributions from molecular dynamics trajectories using GROMACS tools. Reads index files defining interactions and extracts structural properties including averages, standard deviations, and distributions.

**Key features:**
- Automatic PBC removal and trajectory alignment
- Distribution generation for each interaction
- Organized output directory structure
- Support for time series and histogram data
- Automatic cleanup of GROMACS backup files

**Examples:**

```bash
# Basic usage with all three interaction types
python generate_bonds_angles_dihedrals.py \
    --bonds_ndx bonds.ndx \
    --angles_ndx angles.ndx \
    --dihedrals_ndx dihedrals.ndx \
    --xtc_file md.xtc \
    --tpr_file md.tpr \
    --output_all_files analysis_results/

# With PBC removal and alignment (recommended)
python generate_bonds_angles_dihedrals.py \
    --bonds_ndx bonds.ndx \
    --angles_ndx angles.ndx \
    --dihedrals_ndx dihedrals.ndx \
    --xtc_file md.xtc \
    --tpr_file md.tpr \
    --output_all_files analysis_results/ \
    --remove_pbc \
    --index index.ndx \
    --group_1 "Backbone" \
    --group_2 "System"

# Keep intermediate trajectory files for inspection
python generate_bonds_angles_dihedrals.py \
    --bonds_ndx bonds.ndx \
    --angles_ndx angles.ndx \
    --dihedrals_ndx dihedrals.ndx \
    --xtc_file md.xtc \
    --tpr_file md.tpr \
    --output_all_files analysis_results/ \
    --remove_pbc \
    --index index.ndx \
    --group_1 "C-alpha" \
    --keep_intermediate
```

**Index file formats:**

*Bonds index file (bonds.ndx):*
```
[ bonds ]
1 2
2 3
3 4
4 5
```

*Angles index file (angles.ndx):*
```
[ angles ]
1 2 3
2 3 4
3 4 5
```

*Dihedrals index file (dihedrals.ndx):*
```
[ dihedrals ]
1 2 3 4
2 3 4 5
3 4 5 6
```

**Arguments:**

| Argument | Description | Default |
|----------|-------------|---------|
| `--bonds_ndx` | Index file defining bonds | **Required** |
| `--angles_ndx` | Index file defining angles | **Required** |
| `--dihedrals_ndx` | Index file defining dihedrals | **Required** |
| `--xtc_file` | Input trajectory file (.xtc) | **Required** |
| `--tpr_file` | Input topology file (.tpr) | **Required** |
| `--output_all_files` | Base output directory path | **Required** |
| `--index` | Index file for group selection | None |
| `--remove_pbc` | Enable PBC removal and alignment | False |
| `--group_1` | Group used for fitting | None |
| `--group_2` | Group for output | "System" |
| `--keep_intermediate` | Keep intermediate trajectory files | False |

**Output structure:**
```
output_all_files/
├── bonds/
│   ├── report_bonds.txt
│   ├── data_bonds.txt
│   ├── bond_0.xvg
│   ├── distr_bond_0.xvg
│   └── errors.log
├── angles/
│   ├── report_angles.txt
│   ├── data_angles.txt
│   ├── ang_0.xvg
│   ├── distr_ang_0.xvg
│   └── errors.log
└── dihedrals/
    ├── report_dihedrals.txt
    ├── data_dihedrals.txt
    ├── dih_0.xvg
    ├── distr_dih_0.xvg
    └── errors.log
```

---

### 4. bp_distributions.py

**Description:** Computes statistical properties (mean and standard deviation) of bond, angle, and dihedral distributions from GROMACS XVG output files. Processes both time series data and precomputed distributions, saving results as tab-separated values (TSV) files.

**Key features:**
- Automatic detection of time series vs distribution data
- Support for histogram-based distributions
- Parallel processing of bonds, angles, and dihedrals
- TSV output for easy import into spreadsheets or analysis tools

**Examples:**

```bash
# Process all distributions from analysis_results directory
python bp_distributions.py \
    --bonds_dir analysis_results/bonds/ \
    --angles_dir analysis_results/angles/ \
    --dihedrals_dir analysis_results/dihedrals/ \
    --dir_to_output statistics/

# Custom output filenames
python bp_distributions.py \
    --bonds_dir analysis_results/bonds/ \
    --angles_dir analysis_results/angles/ \
    --dihedrals_dir analysis_results/dihedrals/ \
    --bond_out bond_stats.tsv \
    --angle_out angle_stats.tsv \
    --dihedral_out dihedral_stats.tsv \
    --dir_to_output .

# Using default directory names (if structure matches)
python bp_distributions.py \
    --bonds_dir bonds/ \
    --angles_dir angles/ \
    --dihedrals_dir dihedrals/
```

**Arguments:**

| Argument | Description | Default |
|----------|-------------|---------|
| `--bonds_dir` | Directory containing bond XVG files | bonds |
| `--angles_dir` | Directory containing angle XVG files | angles |
| `--dihedrals_dir` | Directory containing dihedral XVG files | dihedrals |
| `--dir_to_output` | Output directory for TSV files | TSV_statistics |
| `--bond_out` | Output filename for bond statistics | bond_statistics.tsv |
| `--angle_out` | Output filename for angle statistics | angle_statistics.tsv |
| `--dihedral_out` | Output filename for dihedral statistics | dihedral_statistics.tsv |

**Output format (TSV):**

*bond_statistics.tsv:*
```
index    mean       sd
0        0.1542     0.0021
1        0.1538     0.0019
2        0.1551     0.0023
```

---

### 5. adaptation_gro_itp.py

**Description:** Updates atom names in a GROMACS ITP file to match those from a reference GRO file. Ensures consistency between topology and coordinate files when atom naming differs between parametrization and structural files.

**Key features:**
- Sequential mapping preserving atom order
- Preserves ITP file formatting
- Verbose mode for debugging
- Error checking for atom count mismatches

**Examples:**

```bash
# Basic adaptation
python adaptation_gro_itp.py \
    --input_itp cg.itp \
    --input_gro_ref molecule.gro \
    --output_itp_adapted cg_adapted.itp

# Verbose mode to see name changes
python adaptation_gro_itp.py \
    --input_itp cg.itp \
    --input_gro_ref molecule.gro \
    --output_itp_adapted cg_adapted.itp \
    --verbose
```

**Arguments:**

| Argument | Description | Default |
|----------|-------------|---------|
| `--input_itp` | Input ITP file to be modified | **Required** |
| `--input_gro_ref` | Reference GRO file with correct atom names | **Required** |
| `--output_itp_adapted` | Output ITP file with updated atom names | **Required** |
| `--verbose` | Print detailed debug information | False |

---

### 6. bp_prep.py

**Description:** Prepares a complete coarse-grained molecular dynamics system using GROMACS and the Martini 3 force field. Automates box definition, solvation with Martini water beads, ion addition for charge neutralization and target salt concentration, and generates energy minimization input files.

**Workflow:**
1. Copies input structure, topology, force field, and MDP files
2. Centers the molecule in the simulation box
3. Defines simulation box (fixed size or distance-based)
4. Solvates the system with Martini water beads
5. Adds ions (Na+/Cl−) for neutralization and target concentration
6. Updates topology with correct molecule counts
7. Generates `.tpr` for energy minimization

**Examples:**

```bash
# Standard Martini setup with distance-based box (recommended)
python bp_prep.py \
    --use_distance_from_atom \
    --distance_from_atom 1.2

# Fixed cubic box of 12 nm
python bp_prep.py \
    --box_size 12.0

# Limit solvent beads (useful for smaller systems)
python bp_prep.py \
    --use_distance_from_atom \
    --distance_from_atom 1.2 \
    --max_solvent 768 \
    --solvent_radius 0.21

# Custom salt concentration (physiological)
python bp_prep.py \
    --use_distance_from_atom \
    --distance_from_atom 1.2 \
    --salt 0.15

# Complete example with custom paths
python bp_prep.py \
    --input_ref_dir /path/to/input \
    --input_gro cg.gro \
    --input_itp cg.itp \
    --input_topol topol_cg.top \
    --input_ff_dir /path/to/ff_files \
    --ff martini_v3.0.0.itp \
    --ions martini_v3.0.0_ions_v1.itp \
    --solvent martini_v3.0.0_solvents_v1.itp \
    --output_dir /path/to/output \
    --water_dir /path/to/water \
    --water_file_gro water.gro \
    --use_distance_from_atom \
    --distance_from_atom 1.2 \
    --max_solvent 768 \
    --salt 0.15
```

**Arguments:**

| Argument | Description | Default |
|----------|-------------|---------|
| `--input_ref_dir` | Directory with input GRO, ITP, TOP | *Custom path* |
| `--input_gro` | Input CG structure (.gro) | cg.gro |
| `--input_itp` | Molecule topology (.itp) | cg.itp |
| `--input_topol` | System topology (.top) | topol_cg.top |
| `--input_ff_dir` | Directory with Martini force field files | None |
| `--output_dir` | Directory for final system files | *Custom path* |
| `--water_dir` | Directory with water bead GRO file | *Custom path* |
| `--water_file_gro` | Water structure file | water.gro |
| `--ff` | Force field file | martini_v3.0.0.itp |
| `--ions` | Ions file | martini_v3.0.0_ions_v1.itp |
| `--solvent` | Solvent file | martini_v3.0.0_solvents_v1.itp |
| `--distance_from_atom` | Distance (nm) between solute and box edge | 2.0 |
| `--use_distance_from_atom` | Use distance-based box calculation | False |
| `--box_size` | Fixed box size (nm) | None |
| `--max_solvent` | Maximum number of solvent beads | 2000 |
| `--solvent_radius` | Martini bead radius (nm) | 0.21 |
| `--salt` | Salt concentration (mol/L) | 0.15 |
| `--input_mdp_dir` | Directory with MDP parameter files | ../../../../mdp/ |
| `--input_name_file_mdp` | MDP filename for minimization | minimization.mdp |
| `--ions_mdp_dir` | Directory with MDP for ion addition | *Custom path* |
| `--ions_file_mdp` | MDP filename for ion addition | ions.mdp |

---

### 7. plot_distributions.py

**Description:** Processes XVG files containing molecular simulation data for bonds, angles, and dihedrals. Compares reference and simulated datasets by generating distribution plots (using Kernel Density Estimation) and time series plots. Automatically organizes outputs into structured directories.

**Key features:**
- Gaussian KDE for smooth probability density curves
- Combined and separate distribution plots
- Time series visualization with mean lines
- Automatic directory creation
- PNG and SVG output formats
- Memory-efficient figure handling

**Examples:**

```bash
# Compare reference and simulated distributions
python plot_distributions.py \
    --bonds_ref_dir reference/bonds/ \
    --bonds_sim_dir simulated/bonds/ \
    --angles_ref_dir reference/angles/ \
    --angles_sim_dir simulated/angles/ \
    --dihedrals_ref_dir reference/dihedrals/ \
    --dihedrals_sim_dir simulated/dihedrals/ \
    --figures_dir comparison_figures/

# Using default directory structure
python plot_distributions.py \
    --bonds_ref_dir bonds_ref/ \
    --bonds_sim_dir bonds_sim/
```

**Output directory structure:**
```
figures_dir/
├── bonds/
│   ├── distr_ref_simulated_joined/
│   │   ├── bond_0.png
│   │   ├── bond_0.svg
│   │   └── ...
│   ├── distr_ref_simulated_separated/
│   │   ├── ref/
│   │   └── simulated/
│   ├── time_series_ref_simulated_joined/
│   └── time_series_ref_simulated_separated/
│       ├── ref/
│       └── simulated/
├── angles/
└── dihedrals/
```

**Arguments:**

| Argument | Description | Default |
|----------|-------------|---------|
| `--bonds_ref_dir` | Directory with reference bond XVG files | bonds |
| `--angles_ref_dir` | Directory with reference angle XVG files | angles |
| `--dihedrals_ref_dir` | Directory with reference dihedral XVG files | dihedrals |
| `--bonds_sim_dir` | Directory with simulated bond XVG files | bonds |
| `--angles_sim_dir` | Directory with simulated angle XVG files | angles |
| `--dihedrals_sim_dir` | Directory with simulated dihedral XVG files | dihedrals |
| `--figures_dir` | Output directory for figures | figures |

---

## Quick Start

Complete pipeline example from AA trajectory to CG simulation ready system:

```bash
# Step 1: Map AA trajectory to CG
python map_aa_traj_to_cg.py \
    --index_cg cg.ndx \
    --aa_tpr md.tpr \
    --aa_xtc md.xtc \
    --output_mapped mapped.xtc \
    --output_cg_gro molecule.gro

# Step 2: Generate CG topology
python generate_cg_top.py \
    --path_ff ff_files/ \
    --itp_ligand cg.itp \
    --name_molecule MOL \
    --number_molecule 1 \
    --output_topol topol_cg.top

# Step 3: Analyze structural properties
python generate_bonds_angles_dihedrals.py \
    --bonds_ndx bonds.ndx \
    --angles_ndx angles.ndx \
    --dihedrals_ndx dihedrals.ndx \
    --xtc_file mapped.xtc \
    --tpr_file md.tpr \
    --output_all_files analysis/

# Step 4: Compute statistics
python bp_distributions.py \
    --bonds_dir analysis/bonds/ \
    --angles_dir analysis/angles/ \
    --dihedrals_dir analysis/dihedrals/ \
    --dir_to_output stats/

# Step 5: Adapt ITP names if needed
python adaptation_gro_itp.py \
    --input_itp cg.itp \
    --input_gro_ref molecule.gro \
    --output_itp_adapted cg_adapted.itp

# Step 6: Prepare simulation system
python bp_prep.py \
    --use_distance_from_atom \
    --distance_from_atom 1.2 \
    --salt 0.15

# Step 7: Generate comparison plots
python plot_distributions.py \
    --bonds_ref_dir reference/bonds/ \
    --bonds_sim_dir analysis/bonds/ \
    --figures_dir figures/
```

---

## Dependencies

- **Python 3.7+** with packages:
  - `numpy`
  - `pandas`
  - `scipy`
  - `matplotlib`
  - `seaborn`
- **GROMACS 2020+** (for modules 1, 3, 6)
- **Standard Unix utilities** (printf, echo, etc.)



