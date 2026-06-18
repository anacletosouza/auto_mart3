# Auto_Mart3 - Coarse-Grained Mapping Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![GROMACS](https://img.shields.io/badge/GROMACS-2021+-green.svg)](https://www.gromacs.org/)

A complete pipeline to map atomistic (AA) molecular dynamics trajectories to coarse-grained (CG) representations, generate CG topologies, prepare GROMACS input files, and perform Bayesian optimization of force constants for Martini 3 force field parameters.

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Required Files](#required-files)
- [Commands](#commands)
  - [auto-map](#auto-map)
  - [auto-gen-top](#auto-gen-top)
  - [auto-analyze](#auto-analyze)
  - [auto-distributions](#auto-distributions)
  - [auto-adapt-itp](#auto-adapt-itp)
  - [auto-prep](#auto-prep)
  - [auto-plot-distributions](#auto-plot-distributions)
  - [bayes-potential-adjust](#bayes-potential-adjust)
  - [autoparam-aa](#autoparam-aa)
  - [autoparam-cg](#autoparam-cg)
- [Complete Workflow Examples](#complete-workflow-examples)
- [Output Directory Structure](#output-directory-structure)

## Overview

**Auto_Mart3** is a comprehensive toolkit for coarse-grained molecular dynamics simulations using the Martini 3 force field. It bridges the gap between atomistic and coarse-grained resolutions by providing:

1. **AA → CG Mapping**: Project atomistic trajectories onto CG beads
2. **Topology Generation**: Create complete GROMACS topology files
3. **Trajectory Analysis**: Extract bonds, angles, and dihedrals from CG simulations
4. **Statistical Analysis**: Compute distribution statistics (mean, std, histograms)
5. **Bayesian Optimization**: Iteratively adjust force constants to match reference distributions
6. **Full Automation**: End-to-end pipelines for parameterization and optimization

## Installation

```bash
# Clone the repository
git clone https://github.com/anacletosouza/auto_mart3.git
cd auto_mart3

# Install using pip
pip install -e .

# Or install dependencies manually
pip install numpy scipy matplotlib pandas


General workflow

cd /grain/anacleto/projects/project_1_virion_simulation/CG_model/spike_protein/CG-carbohydrate_types

(1) Generate the initial file using pipeline auto_mart_aa. If you want to use different CPU affinity (cpu 127), we suggest:

export OMP_NUM_THREADS=1
export OMP_PROC_BIND=true
export OMP_PLACES="{127}"
export GOMP_CPU_AFFINITY="127"

numactl --physcpubind=127 \
  auto_mart3 auto_mart_aa \
  --aa_tpr AA_simulations_3/O_glycan_SER_A/md.tpr \
  --aa_xtc AA_simulations_3/O_glycan_SER_A/md.xtc \
  --aa_gro AA_simulations_3/O_glycan_SER_A/carb.gro \
  --aa_itp AA_simulations_3/O_glycan_SER_A/carb_filtered.itp \
  --beads_json AA_simulations/json/beads_config.json \
  --output_dir AUTO_MART3_AA/results_O_glycan_SER_A \
  --name_molecule molecule


(2) If you want to use different CPU affinity, we suggest:

export OMP_NUM_THREADS=32
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export GOMP_CPU_AFFINITY="$(seq -s ' ' 1 2 63)"

numactl --physcpubind=$(seq -s, 1 2 63) auto_mart3 auto_mart_cg   --INPUT_AUTO_MART_AA_DIR AUTO_MART3_AA/results_M8   --OUTPUT_AUTO_MART_CG_DIR AUTO_MART3_CG/results_M8   --dt_nvt_ps 0.0005 --time_nvt_ps 10000   --dt_npt_ps 0.0005 --time_npt_ps 10000   --dt_md_ps 0.0005 --time_md_ps 100000   --pin_on --ntomp $OMP_NUM_THREADS --cuda_visible_devices 0

export OMP_NUM_THREADS=32
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export GOMP_CPU_AFFINITY="$(seq -s ' ' 64 2 126)"

numactl --physcpubind=$(seq -s, 64 2 126) auto_mart3 auto_mart_cg   --INPUT_AUTO_MART_AA_DIR AUTO_MART3_AA/results_O_glycan_SER_A   --OUTPUT_AUTO_MART_CG_DIR AUTO_MART3_CG/results_O_glycan_SER_A   --dt_nvt_ps 0.0005 --time_nvt_ps 10000   --dt_npt_ps 0.0005 --time_npt_ps 10000   --dt_md_ps 0.0005 --time_md_ps 100000   --pin_on --ntomp $OMP_NUM_THREADS --cuda_visible_devices 0

export OMP_NUM_THREADS=32
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export GOMP_CPU_AFFINITY="$(seq -s ' ' 65 2 127)"

numactl --physcpubind=$(seq -s, 65 2 127) auto_mart3 auto_mart_cg   --INPUT_AUTO_MART_AA_DIR AUTO_MART3_AA/results_O_glycan_THR_A   --OUTPUT_AUTO_MART_CG_DIR AUTO_MART3_CG/results_O_glycan_THR_A   --dt_nvt_ps 0.0005 --time_nvt_ps 10000   --dt_npt_ps 0.0005 --time_npt_ps 10000   --dt_md_ps 0.0005 --time_md_ps 100000   --pin_on --ntomp $OMP_NUM_THREADS --cuda_visible_devices 1
  

```

**Dependencies:**
- Python 3.8+
- GROMACS (2021 or later)
- NumPy, SciPy, Matplotlib, Pandas

## Required Files

Before running Auto_Mart3, ensure you have the following files:

| File Type | Description | Example |
|-----------|-------------|---------|
| `.tpr` | GROMACS run input file (solvent+molecule+ions) | `md.tpr` |
| `.xtc` | Compressed trajectory (solvent+molecule+ions) | `traj.xtc` |
| `.gro` | Structure file (aa model for molecule (no solvent and no ions) | `conf.gro` |
| `.itp` | Molecular topology (aa model for molecule- no solvent and no ions) | `molecule.itp` |
| `.ndx` | Index file (bonds, angles, dihedrals) | `bonds.ndx` |
| `.json` | Beads definition for mapping | `beads.json` |
| `.mdp` | GROMACS parameters | `minimization.mdp` |

## Commands

### auto-map

Map atomistic trajectory to coarse-grained (CG) resolution using GROMACS.

**Description:** Projects atomistic coordinates onto CG beads based on index definitions, creating a CG trajectory and structure file.

**Inputs:**
| Argument | Required | Description |
|----------|----------|-------------|
| `--index_cg` | Yes | CG index file (.ndx) containing bead definitions |
| `--aa_tpr` | Yes | Atomistic topology file (.tpr) from AA simulation (with solvent+ions) |
| `--aa_xtc` | Yes | Atomistic trajectory file (.xtc) from AA simulation (with solvent+ions) |
| `--output_mapped` | No | Output CG trajectory filename (default: mapped.xtc) (with solvent+ions) |
| `--output_cg_gro` | No | Output CG .gro file (default: molecule.gro) (no solvent+ions) |
| `--remove_pbc` | No | Remove PBC using `gmx trjconv -pbc whole` (default: True) (with solvent+ions) |
| `--gmx_cmd` | No | GROMACS command to use (default: gmx) |

**Outputs:**
- CG trajectory file (`.xtc`)
- CG structure file (`.gro`)

**Example:**
```bash
auto_mart3 auto-map \
    --index_cg CG_MARTINI3/NDX/cg.ndx \
    --aa_tpr ../AA_sim/md.tpr \
    --aa_xtc ../AA_sim/traj.xtc \
    --output_mapped GMX/mapped.xtc \
    --output_cg_gro GMX/cg.gro \
    --remove_pbc
```

---

### auto-gen-top

Generate GROMACS topology (.top) file for CG simulations.

**Description:** Creates a complete system topology by combining force field files (Martini 3, ions, solvents) with the CG molecule topology.

**Inputs:**
| Argument | Required | Description |
|----------|----------|-------------|
| `--path_ff` | Yes | Directory containing force field ITP files |
| `--itp_ligand` | Yes | Ligand ITP file (use 'None' to skip) |
| `--name_molecule` | Yes | Molecule name for the ligand |
| `--number_molecule` | Yes | Number of molecules in the system |
| `--output_topol` | Yes | Output topology filename |
| `--ff` | No | Force field ITP (default: martini_v3.0.0.itp) |
| `--ions` | No | Ions ITP (default: martini_v3.0.0_ions_v1.itp) |
| `--solvent` | No | Solvent ITP (default: martini_v3.0.0_solvents_v1.itp) |

**Outputs:**
- Complete system topology file (`.top`)

**Example:**
```bash
auto_mart3 auto-gen-top \
    --path_ff ./ff_files \
    --itp_ligand cg.itp \
    --name_molecule protein \
    --number_molecule 1 \
    --output_topol topol_cg.top
```

---

### auto-analyze

Calculate bonds, angles, and dihedrals from MD trajectories.

**Description:** Extracts bond distances, bond angles, and dihedral angles from CG trajectories, saving time-series data in XVG format.

**Inputs:**
| Argument | Required | Description |
|----------|----------|-------------|
| `--bonds_ndx` | Yes | Index file defining bonds |
| `--angles_ndx` | Yes | Index file defining angles |
| `--dihedrals_ndx` | Yes | Index file defining dihedrals |
| `--xtc_file` | Yes | Input trajectory file (.xtc) |
| `--tpr_file` | Yes | Input topology file (.tpr) |
| `--output_all_files` | Yes | Base output directory path |
| `--index` | No | Index file for group selection |
| `--remove_pbc` | No | Remove PBC and align trajectory |
| `--group_1` | No | Group for fitting/alignment |
| `--group_2` | No | Group for output (default: System) |

**Outputs:**
- `bonds/` - Bond distance vs time (XVG files)
- `angles/` - Angle values vs time (XVG files)
- `dihedrals/` - Dihedral angles vs time (XVG files)
- `*_histograms/` - Histogram data files

**Example:**
```bash
auto_mart3 auto-analyze \
    --bonds_ndx CG_MARTINI3/NDX/bonds.ndx \
    --angles_ndx CG_MARTINI3/NDX/angles.ndx \
    --dihedrals_ndx CG_MARTINI3/NDX/dihedrals.ndx \
    --xtc_file GMX/mapped.xtc \
    --tpr_file GMX/CG.tpr \
    --output_all_files BONDS_ANGLES_DIHEDRALS_XVG_REF
```

---

### auto-distributions

Calculate statistical properties from XVG distribution files.

**Description:** Computes mean, standard deviation, histogram bins, and counts for bond, angle, and dihedral distributions.

**Inputs:**
| Argument | Required | Description |
|----------|----------|-------------|
| `--bonds_dir` | No | Directory with bond XVG files (default: bonds) |
| `--angles_dir` | No | Directory with angle XVG files (default: angles) |
| `--dihedrals_dir` | No | Directory with dihedral XVG files (default: dihedrals) |
| `--dir_to_output` | No | Output directory (default: TSV_statistics) |
| `--bond_out` | No | Bond statistics filename (default: bond_statistics.tsv) |
| `--angle_out` | No | Angle statistics filename (default: angle_statistics.tsv) |
| `--dihedral_out` | No | Dihedral statistics filename (default: dihedral_statistics.tsv) |

**Outputs:**
- TSV files with statistics (mean, std, bins, counts)

**Example:**
```bash
auto_mart3 auto-distributions \
    --bonds_dir BONDS_ANGLES_DIHEDRALS_XVG_REF/bonds \
    --angles_dir BONDS_ANGLES_DIHEDRALS_XVG_REF/angles \
    --dihedrals_dir BONDS_ANGLES_DIHEDRALS_XVG_REF/dihedrals \
    --dir_to_output STATISTICS
```

---

### auto-adapt-itp

Update atom names in ITP file to match GRO reference.

**Description:** Ensures consistency between topology and structure files by renaming atoms in the ITP to match those in the GRO file.

**Inputs:**
| Argument | Required | Description |
|----------|----------|-------------|
| `--input_itp` | Yes | Input ITP file to modify |
| `--input_gro_ref` | Yes | Reference GRO file with correct atom names |
| `--output_itp_adapted` | Yes | Output ITP file with updated names |

**Outputs:**
- Adapted ITP file with matching atom names

**Example:**
```bash
auto_mart3 auto-adapt-itp \
    --input_itp GMX/cg.itp \
    --input_gro_ref GMX/cg.gro \
    --output_itp_adapted GMX/cg.itp.adapted
```

---

### auto-prep

Setup CG system with solvent and ions using Martini 3 force field.

**Description:** Prepares a complete solvated CG system ready for MD simulation, including adding water, ions, energy minimization, and creating run scripts.

**Inputs:**
| Argument | Required | Description |
|----------|----------|-------------|
| `--input_ref_dir` | Yes | Directory with input files |
| `--output_dir` | Yes | Output directory |
| `--input_gro` | No | CG structure filename (default: cg.gro) |
| `--input_itp` | No | CG topology ITP (default: cg.itp) |
| `--input_topol` | No | System topology (default: topol_cg.top) |
| `--input_ff_dir` | No | Force field directory |
| `--use_distance_from_atom` | No | Use distance to calculate box size |
| `--distance_from_atom` | No | Distance (nm) from molecule (default: 2.0) |
| `--salt` | No | Salt concentration in M (default: 0.15) |

**Outputs:**
- `topol.top` - Final topology with solvent/ions
- `conf.gro` - Solvated CG structure
- `em.mdp`, `em.tpr`, `em.gro` - Energy minimization files
- `md.mdp`, `md.tpr` - MD simulation files
- `run_simulation.sh` - Simulation script

**Example:**
```bash
auto_mart3 auto-prep \
    --input_ref_dir GMX \
    --output_dir MDRUN_CG \
    --input_gro cg.gro \
    --input_itp cg.itp \
    --input_topol topol_cg.top \
    --input_ff_dir ./ff_files \
    --use_distance_from_atom \
    --distance_from_atom 2.0 \
    --salt 0.15
```

---

### auto-plot-distributions

Plot reference vs simulated distributions from XVG files.

**Description:** Generates comparison plots between reference (AA-derived) and simulated (CG) distributions for bonds, angles, and dihedrals.

**Inputs:**
| Argument | Required | Description |
|----------|----------|-------------|
| `--bonds_ref_dir` | No | Reference bond XVG directory (default: bonds) |
| `--angles_ref_dir` | No | Reference angle XVG directory (default: angles) |
| `--dihedrals_ref_dir` | No | Reference dihedral XVG directory (default: dihedrals) |
| `--bonds_sim_dir` | No | Simulated bond XVG directory (default: bonds) |
| `--angles_sim_dir` | No | Simulated angle XVG directory (default: angles) |
| `--dihedrals_sim_dir` | No | Simulated dihedral XVG directory (default: dihedrals) |
| `--figures_dir` | No | Output directory for figures (default: figures) |

**Outputs:**
- PNG/SVG figures comparing reference vs simulated distributions

**Example:**
```bash
auto_mart3 auto-plot-distributions \
    --bonds_ref_dir INPUT/BONDS_REF/bonds \
    --bonds_sim_dir OUTPUT/BONDS_SIM/bonds \
    --figures_dir FIGURES_COMPARISON
```

---

### bayes-potential-adjust

Adjust CG force constants using Bayesian update with R² correction.

**Description:** Implements simulated annealing with Bayesian inference to optimize force constants, minimizing the difference between reference and simulated distributions using R² as the objective function.

**Inputs:**
| Argument | Required | Description |
|----------|----------|-------------|
| `--bonds_ref_xvg_dir` | Yes | Reference bond XVG directory |
| `--angles_ref_xvg_dir` | Yes | Reference angle XVG directory |
| `--dihedrals_ref_xvg_dir` | Yes | Reference dihedral XVG directory |
| `--bonds_sim_xvg_dir` | Yes | Simulated bond XVG directory |
| `--angles_sim_xvg_dir` | Yes | Simulated angle XVG directory |
| `--dihedrals_sim_xvg_dir` | Yes | Simulated dihedral XVG directory |
| `--itp_cg` | Yes | CG ITP file with topology |
| `--ndx_bounds` | Yes | NDX file with bond indices |
| `--ndx_angles` | Yes | NDX file with angle indices |
| `--ndx_dihedrals` | Yes | NDX file with dihedral indices |
| `--molecule_name` | No | Molecule name (default: molecule) |
| `--T0` | No | Initial SA temperature (default: 10.0) |
| `--alpha` | No | Cooling factor (default: 0.85) |
| `--distribution_points` | No | Points for R² (default: 200) |
| `--itp_out` | No | Output ITP filename (default: cg.itp) |

**Outputs:**
- Updated ITP file with optimized force constants
- Console output with R² values per iteration

**Example:**
```bash
auto_mart3 bayes-potential-adjust \
    --bonds_ref_xvg_dir REF/bonds \
    --angles_ref_xvg_dir REF/angles \
    --dihedrals_ref_xvg_dir REF/dihedrals \
    --bonds_sim_xvg_dir SIM/bonds \
    --angles_sim_xvg_dir SIM/angles \
    --dihedrals_sim_xvg_dir SIM/dihedrals \
    --itp_cg cg.itp \
    --ndx_bounds bonds.ndx \
    --ndx_angles angles.ndx \
    --ndx_dihedrals dihedrals.ndx \
    --molecule_name protein \
    --T0 10.0 --alpha 0.85 \
    --itp_out cg_optimized.itp
```

---

### auto_mart_aa

Run full AA to CG parametrization pipeline.

**Description:** Complete end-to-end pipeline for mapping atomistic simulations to CG representation, including trajectory mapping, topology generation, and distribution analysis.

**Required Arguments:**
| Argument | Description |
|----------|-------------|
| `--aa_tpr` | AA TPR file (with solvent+ions) |
| `--aa_xtc` | AA trajectory XTC file (with solvent+ions) |
| `--aa_gro` | AA GRO file (ONLY molecule of interest)|
| `--aa_itp` | AA ITP file (from molecule of interest)|
| `--beads_json` | Beads definition JSON file |
| `--input_mdp` | MDP file for grompp (mypathway/em.mdp for example) |
| `--path_ff` | Force field directory |
| `--output_dir` | Output directory |
| `--name_molecule` | Molecule name |

**Optional Arguments:**
| Argument | Default | Description |
|----------|---------|-------------|
| `--force_application` | random=[1250,30;30,1] | Force application method |
| `--beads_position` | com | Beads positioning (com/cog/min) |
| `--cycle_restr` | fix=3,mode=cycle | Cycle restraints |
| `--maxwarn` | 2 | Max warnings for grompp |
| `--distance_from_atom` | 2.0 | Distance for solvent (nm) |
| `--salt` | 0.15 | Salt concentration (M) |

**Output Directory Structure:**
```
output_dir/
├── CG_MARTINI3/          # Initial CG mapping
├── GMX/                  # GROMACS-compatible files
├── BONDS_ANGLES_DIHEDRALS_XVG_REF/  # Analysis output
├── STATISTICS/           # Distribution statistics
└── MDRUN_CG/            # Simulation-ready system
```

**Example:**
```bash
auto_mart3 autoparam-aa \
                         --aa_tpr ../AA_sim/md.tpr \
                         --aa_xtc ../AA_sim/traj.xtc \
                         --aa_gro ../AA_sim/conf.gro \
                         --aa_itp molecule.itp \
                         --beads_json beads_definition.json \
                         --input_mdp minimization.mdp \
                         --path_ff ./martini_3_ff \
                         --output_dir ./results_AA2CG \
                         --name_molecule protein \
                         --distance_from_atom 2.5 \
                         --salt 0.15
    
auto_mart3 auto_mart_aa \
                        --aa_tpr AA_simulations_1/A2/md.tpr \
                        --aa_xtc AA_simulations_1/A2/md.xtc \
                        --aa_gro AA_simulations_1/A2/carb.gro \
                        --aa_itp AA_simulations_1/A2/carb_filtered.itp \
                        --beads_json AA_simulations/json/beads_config.json \
                        --input_mdp AA_simulations/mdp/minimization.mdp \
                        --path_ff AA_simulations/ff_files \
                        --output_dir AUTOPARAM_AA_2/results_A2 \
                        --name_molecule molecule
```

---

### autoparam-cg

Run CG parameter optimization pipeline (iterative Bayesian optimization).

**Description:** Performs iterative optimization of bond, angle, and dihedral force constants by comparing reference distributions from AA simulations with CG simulations using Bayesian inference and simulated annealing.

**Required Arguments:**
| Argument | Description |
|----------|-------------|
| `--INPUT_AUTOPARAM_AA_DIR` | Directory with autoparam-aa results |
| `--OUTPUT_AUTOPARAM_CG_DIR` | Directory for CG optimization output |

**Optional Arguments:**
| Category | Argument | Default | Description |
|----------|----------|---------|-------------|
| **Simulation** | `--ntomp` | 10 | OpenMP threads |
| | `--ntmpi` | 1 | MPI threads |
| | `--ref_t` | 310 | Reference temperature (K) |
| | `--ref_p` | 1.0 | Reference pressure (bar) |
| | `--dt_nvt_ps` | 0.001 | NVT timestep (ps) |
| | `--time_nvt_ps` | 1000 | NVT duration (ps) |
| | `--dt_npt_ps` | 0.001 | NPT timestep (ps) |
| | `--time_npt_ps` | 5000 | NPT duration (ps) |
| | `--dt_md_ps` | 0.002 | MD timestep (ps) |
| | `--time_md_ps` | 10000 | MD duration (ps) |
| **Force Limits** | `--min_force_bond` | 500.0 | Min bond force constant |
| | `--max_force_bond` | 50000.0 | Max bond force constant |
| | `--default_force_bond` | 1250.0 | Default bond force |
| | `--min_force_angle` | 10.0 | Min angle force |
| | `--max_force_angle` | 150.0 | Max angle force |
| | `--default_force_angle` | 25.0 | Default angle force |
| | `--min_force_dihedral` | 10.0 | Min dihedral force |
| | `--max_force_dihedral` | 150.0 | Max dihedral force |
| | `--default_force_dihedral` | 25.0 | Default dihedral force |
| **Optimization** | `--T0` | 10.0 | Initial SA temperature |
| | `--alpha` | 0.85 | Cooling factor |
| | `--distribution_points` | 100 | Points for R² |
| | `--n_iter` | 30 | Number of iterations |
| **GROMACS** | `--group_2` | System | Output group |
| | `--group_1` | System | Center group |
| | `--index` | None | Custom index file |
| | `--no_remove_gmx_files_in_iter` | False | Keep GMX files |

**Output Directory Structure:**
```
OUTPUT_AUTOPARAM_CG_DIR/
└── CALIBRATION/
    ├── iter_0/
    │   ├── MD_CG/              # Initial simulation
    │   └── RESULTS/            # Analysis results
    ├── iter_1/
    ├── iter_2/
    └── ...
```

**Example:**
```bash
# Basic usage
auto_mart3 autoparam-cg \
    --INPUT_AUTOPARAM_AA_DIR ./results_AA2CG \
    --OUTPUT_AUTOPARAM_CG_DIR ./results_CG_opt

# Full customization
auto_mart3 autoparam-cg \
    --INPUT_AUTOPARAM_AA_DIR ./results_AA2CG \
    --OUTPUT_AUTOPARAM_CG_DIR ./results_CG_opt \
    --ntomp 20 --ntmpi 2 \
    --ref_t 300 --ref_p 1.01325 \
    --time_nvt_ps 2000 --time_npt_ps 10000 --time_md_ps 20000 \
    --n_iter 50 --T0 20.0 --alpha 0.9 \
    --min_force_bond 1000 --max_force_bond 100000 \
    --group_2 "Protein" --group_1 "Protein" \
    --index custom_index.ndx \
    --no_remove_gmx_files_in_iter \
    --molecule_name "my_protein"
```

---

## Complete Workflow Examples

### Example 1: Full AA to CG Parameterization

```bash
# Step 1: Run AA to CG mapping and initial setup
auto_mart3 autoparam-aa \
    --aa_tpr AA_simulation/md.tpr \
    --aa_xtc AA_simulation/traj.xtc \
    --aa_gro AA_simulation/conf.gro \
    --aa_itp molecule.itp \
    --beads_json beads_definition.json \
    --input_mdp minimization.mdp \
    --path_ff ./martini_3_ff \
    --output_dir ./results_AA2CG \
    --name_molecule protein

# Step 2: Run CG optimization
auto_mart3 autoparam-cg \
    --INPUT_AUTOPARAM_AA_DIR ./results_AA2CG \
    --OUTPUT_AUTOPARAM_CG_DIR ./results_CG_opt
```

### Example 2: Individual Module Usage

```bash
# Map AA trajectory to CG
auto_mart3 auto-map \
    --index_cg cg.ndx \
    --aa_tpr md.tpr \
    --aa_xtc traj.xtc \
    --output_mapped cg_traj.xtc

# Analyze CG trajectory
auto_mart3 auto-analyze \
    --bonds_ndx bonds.ndx \
    --angles_ndx angles.ndx \
    --dihedrals_ndx dihedrals.ndx \
    --xtc_file cg_traj.xtc \
    --tpr_file cg.tpr \
    --output_all_files analysis/

# Optimize force constants
auto_mart3 bayes-potential-adjust \
    --bonds_ref_xvg_dir ref/bonds \
    --angles_ref_xvg_dir ref/angles \
    --dihedrals_ref_xvg_dir ref/dihedrals \
    --bonds_sim_xvg_dir sim/bonds \
    --angles_sim_xvg_dir sim/angles \
    --dihedrals_sim_xvg_dir sim/dihedrals \
    --itp_cg molecule.itp \
    --ndx_bounds bonds.ndx \
    --ndx_angles angles.ndx \
    --ndx_dihedrals dihedrals.ndx \
    --itp_out molecule_optimized.itp
```

---

## Output Directory Structure

### autoparam-aa Output
```
results_AA2CG/
├── CG_MARTINI3/
│   ├── NDX/
│   │   ├── cg.ndx
│   │   ├── bonds.ndx
│   │   ├── angles.ndx
│   │   └── dihedrals.ndx
│   ├── ITP/
│   │   └── final_cg.itp
│   ├── GRO/
│   │   └── cg.gro
│   └── PDB/
│       └── cg.pdb
├── GMX/
│   ├── mapped.xtc
│   ├── cg.gro
│   ├── cg.itp
│   ├── CG.tpr
│   └── topol_cg.top
├── BONDS_ANGLES_DIHEDRALS_XVG_REF/
│   ├── bonds/
│   ├── angles/
│   ├── dihedrals/
│   ├── bonds_histograms/
│   ├── angles_histograms/
│   └── dihedrals_histograms/
├── STATISTICS/
│   ├── bond_statistics.tsv
│   ├── angle_statistics.tsv
│   └── dihedral_statistics.tsv
└── MDRUN_CG/
    ├── topol.top
    ├── conf.gro
    ├── em.mdp, em.tpr, em.gro
    ├── md.mdp, md.tpr
    └── run_simulation.sh
```

### autoparam-cg Output
```
results_CG_opt/
└── CALIBRATION/
    ├── iter_0/
    │   ├── MD_CG/           # Initial simulation files
    │   └── RESULTS/
    │       ├── BONDS_ANGLES_DIHEDRALS_XVG_SIM/
    │       │   ├── bonds/
    │       │   ├── angles/
    │       │   └── dihedrals/
    │       └── FIGURES_DISTR_SERIES_REF_VS_SIM/
    ├── iter_1/
    ├── iter_2/
    └── ...
```

---

## License

MIT License - see LICENSE file for details.

## Citation

If you use Auto_Mart3 in your research, please cite:

```
Anacleto Silva de Souza, Auto_Mart3: Automated Coarse-Grained Mapping and Parameterization Pipeline,
https://github.com/anacletosouza/auto_mart3, 2026.
```
