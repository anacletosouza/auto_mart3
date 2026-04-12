# Bayesian Potentials - Coarse-Grained Mapping Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

A complete pipeline to map atomistic (AA) molecular dynamics trajectories to coarse-grained (CG) representations, generate CG topologies, prepare GROMACS input files, and analyze bonds/angles/dihedrals from CG trajectories.

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Required Files](#required-files)
- [Commands](#commands)
  - [`bayesian-potentials map`](#bayesian-potentials-map)
  - [`bayesian-potentials gen-top`](#bayesian-potentials-gen-top)
  - [`bayesian-potentials analyze`](#bayesian-potentials-analyze)
  - [`bayesian-potentials pipeline`](#bayesian-potentials-pipeline)
- [Complete Examples](#complete-examples)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)

## Overview

This tool converts atomistic MD simulations to coarse-grained representations (Martini 3 compatible) and provides comprehensive analysis of CG trajectories. The pipeline:

1. **Generates CG mapping** using `cg-martini3` (beads definition, index files, topology)
2. **Maps** AA trajectories to CG beads based on index definitions
3. **Creates** complete GROMACS topology and runs `grompp` for simulation readiness
4. **Analyzes** bonds, angles, and dihedrals from CG trajectories

## Installation

```bash
# Clone the repository
git clone https://github.com/anacletosouza/bayesian-potentials.git
cd bayesian-potentials

# Install the package
pip install -e .

# Or install directly
pip install bayesian-potentials
```

### Additional Requirements

The package utilizes modules from `cg-martini3`. Install it as follows:

```bash
# Clone the repository
git clone https://github.com/anacletosouza/carb_param_automated_martini3.git
cd carb_param_automated_martini3

# Install the package
pip install -e .

# Or install directly
pip install carb_param_automated_martini3

# Or install via git
pip install git+ssh://git@github.com/anacletosouza/carb_param_automated_martini3.git
```

## Required Files

Before running the pipeline, prepare these files:

| File | Description | How to obtain |
|------|-------------|----------------|
| `md.tpr` | AA GROMACS input file | `gmx grompp -f md.mdp -c conf.gro -p topol.top -o md.tpr` |
| `md.xtc` | AA trajectory (XTC/TRR) | Output from MD simulation |
| `md.gro` | AA structure file | Initial structure or last frame |
| `carb.itp` | AA molecule ITP | Original force field ITP |
| `beads_config.json` | Beads definition JSON | Define bead groups from AA atoms |
| `minimization.mdp` | GROMACS MDP parameters | Included in package data |
| `ff_files/` | Martini3 force field files | Included in package data |

### Beads Configuration File (`beads_config.json`)

Example configuration for a molecule with multiple beads:

```json
{
  "bead_groups": [
    {
      "name": "BB",
      "atoms": [1, 2, 3, 4, 12, 13],
      "description": "Backbone bead"
    },
    {
      "name": "SC1",
      "atoms": [5, 6, 7, 8, 9, 10, 11],
      "description": "Side chain bead 1"
    },
    {
      "name": "BG1",
      "atoms": [14, 15, 16, 17, 18],
      "description": "Bead group 1"
    },
    {
      "name": "BG2",
      "atoms": [19, 20, 21, 22, 23, 24, 25, 26, 27, 28],
      "description": "Bead group 2"
    }
  ]
}
```

## Commands

### `bayesian-potentials map`

Maps an atomistic trajectory to coarse-grained coordinates using an existing CG index file.

```bash
bayesian-potentials map \
    --index_cg NDX/cg.ndx \
    --aa_tpr setup/md.tpr \
    --aa_xtc setup/md.xtc \
    --output_mapped mapped.xtc \
    --output_cg_gro cg.gro \
    --remove_pbc \
    --verbose
```

#### Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--index_cg, -i` | ✓ | CG index file with bead definitions |
| `--aa_tpr, -t` | ✓ | Atomistic TPR file |
| `--aa_xtc, -x` | ✓ | Atomistic trajectory |
| `--output_mapped, -o` | | Output mapped trajectory (default: mapped.xtc) |
| `--output_cg_gro` | | Output CG GRO file (default: molecule.gro) |
| `--remove_pbc` | | Remove periodic boundary conditions (default: True) |
| `--verbose, -v` | | Verbose output |

### `bayesian-potentials gen-top`

Generates a complete GROMACS topology file including force field and molecule.

```bash
bayesian-potentials gen-top \
    --path_ff ff_files \
    --ff martini_v3.0.0.itp \
    --ions martini_v3.0.0_ions_v1.itp \
    --solvent martini_v3.0.0_solvents_v1.itp \
    --itp_ligand cg.itp \
    --name_molecule "MOL" \
    --number_molecule 1 \
    --output_topol topol.top
```

#### Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--path_ff` | ✓ | Directory with force field ITP files |
| `--ff` | | FF filename (default: martini_v3.0.0.itp) |
| `--ions` | | Ions ITP (default: martini_v3.0.0_ions_v1.itp) |
| `--solvent` | | Solvent ITP (default: martini_v3.0.0_solvents_v1.itp) |
| `--itp_ligand` | ✓ | CG ITP file for your molecule |
| `--name_molecule` | ✓ | Molecule name in topology |
| `--number_molecule` | ✓ | Number of molecules |
| `--output_topol` | ✓ | Output topology file |

### `bayesian-potentials analyze`

Analyzes bonds, angles, and dihedrals from CG trajectories.

```bash
bayesian-potentials analyze \
    --bonds_ndx NDX/bonds.ndx \
    --angles_ndx NDX/angles.ndx \
    --dihedrals_ndx NDX/dihedrals.ndx \
    --xtc_file mapped.xtc \
    --tpr_file CG.tpr \
    --remove_pbc \
    --group_1 "System" \
    --group_2 "System"
```

#### Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--bonds_ndx` | ✓ | - | Bonds index file |
| `--angles_ndx` | ✓ | - | Angles index file |
| `--dihedrals_ndx` | ✓ | - | Dihedrals index file |
| `--xtc_file` | ✓ | - | XTC trajectory file |
| `--tpr_file` | ✓ | - | TPR topology file |
| `--index` | | None | Index file for group selection |
| `--remove_pbc` | | False | Remove PBC and align trajectory |
| `--group_1` | | "System" | Group for fitting |
| `--group_2` | | "System" | Group for output |
| `--keep_intermediate` | | False | Keep intermediate files |

### `bayesian-potentials pipeline`

**Complete pipeline** - generates CG mapping, maps trajectory, creates topology, runs grompp, and analyzes bonds/angles/dihedrals in one command.

```bash
bayesian-potentials pipeline \
    --aa_tpr setup/md.tpr \
    --aa_xtc setup/md.xtc \
    --aa_gro setup/md.gro \
    --aa_itp setup/carb.itp \
    --beads_json json/beads_config.json \
    --force_application "random=[1250,30;30,4]" \
    --beads_position com \
    --cycle_restr "fix=3,mode=cycle" \
    --input_mdp mdp/minimization.mdp \
    --path_ff ff_files/
```

#### Arguments

| Category | Argument | Default | Description |
|----------|----------|---------|-------------|
| **Required** | `--aa_tpr` | - | AA TPR file |
| | `--aa_xtc` | - | AA trajectory XTC file |
| | `--aa_gro` | - | AA GRO file |
| | `--aa_itp` | - | AA ITP file |
| | `--beads_json` | - | Beads definition JSON file |
| | `--force_application` | - | Force application (random or fixed) |
| | `--beads_position` | - | Beads position ('com' or 'geom') |
| | `--input_mdp` | - | MDP file for grompp |
| | `--path_ff` | - | Force field directory |
| **Optional** | `--output_dir` | `results` | Output directory |
| | `--name_molecule` | `molecule` | Molecule name |
| | `--number_molecule` | `1` | Number of molecules |
| | `--ff` | `martini_v3.0.0.itp` | Force field ITP filename |
| | `--ions` | `martini_v3.0.0_ions_v1.itp` | Ions ITP filename |
| | `--solvent` | `martini_v3.0.0_solvents_v1.itp` | Solvent ITP filename |
| | `--title_comments` | `Topology system in Martini 3` | Topology comments |
| | `--title_system` | `molecule in aqueous solution` | System title |
| | `--output_topol` | `topol_cg.top` | Output topology filename |
| | `--cycle_restr` | `fix=3,mode=cycle` | Cycle constraints |
| | `--default_martini` | `False` | Use default Martini3 parameters |
| | `--maxwarn` | `1` | Max warnings for grompp |
| | `--remove_pbc` | `True` | Remove PBC during mapping |
| | `--skip_grompp` | `False` | Skip grompp step |
| | `--skip_analysis` | `False` | Skip bonds/angles/dihedrals analysis |
| | `--analyze_remove_pbc` | `False` | Remove PBC before analysis |
| | `--analyze_group_1` | `System` | Group for fitting in analysis |
| | `--analyze_group_2` | `System` | Group for output in analysis |
| | `--keep_intermediate` | `False` | Keep intermediate analysis files |
| | `--keep_temp` | `False` | Keep temporary files |
| | `--verbose` | `False` | Verbose output |

#### Cycle Constraints (`--cycle_restr`)

The `--cycle_restr` parameter controls how internal bonds within each residue are determined:

| Value | Description |
|-------|-------------|
| `"none"` | Preserve original ITP connections (no modifications) |
| `"fix=3,mode=cycle"` | Create 3-member cycles (triangles) - **default** |
| `"fix=4,mode=cycle"` | Create 4-member cycles (squares/rings) |
| `"fix=5,mode=cycle"` | Create 5-member cycles |
| `"mode=linear"` | Force linear chains (Minimum Spanning Tree) |

## Complete Examples

### Example 1: Default Pipeline (with cycle constraints)

```bash
bayesian-potentials pipeline \
    --aa_tpr setup/md.tpr \
    --aa_xtc setup/md.xtc \
    --aa_gro setup/md.gro \
    --aa_itp setup/carb.itp \
    --beads_json json/beads_config.json \
    --force_application "random=[1250,30;30,4]" \
    --beads_position com \
    --input_mdp mdp/minimization.mdp \
    --path_ff ff_files/
```

### Example 2: Preserve Original ITP Connections

```bash
bayesian-potentials pipeline \
    --aa_tpr setup/md.tpr \
    --aa_xtc setup/md.xtc \
    --aa_gro setup/md.gro \
    --aa_itp setup/carb.itp \
    --beads_json json/beads_config.json \
    --force_application "random=[1250,30;30,4]" \
    --beads_position com \
    --cycle_restr "none" \
    --input_mdp mdp/minimization.mdp \
    --path_ff ff_files/
```

### Example 3: Linear Mode (no cycles)

```bash
bayesian-potentials pipeline \
    --aa_tpr setup/md.tpr \
    --aa_xtc setup/md.xtc \
    --aa_gro setup/md.gro \
    --aa_itp setup/carb.itp \
    --beads_json json/beads_config.json \
    --force_application "fix=[1250;25]" \
    --beads_position geom \
    --cycle_restr "mode=linear" \
    --input_mdp mdp/minimization.mdp \
    --path_ff ff_files/
```

### Example 4: Full Pipeline with Verbose Output

```bash
bayesian-potentials pipeline \
    --aa_tpr setup/md.tpr \
    --aa_xtc setup/md.xtc \
    --aa_gro setup/md.gro \
    --aa_itp setup/carb.itp \
    --beads_json json/beads_config.json \
    --force_application "random=[1250,30;30,4]" \
    --beads_position com \
    --cycle_restr "fix=3,mode=cycle" \
    --input_mdp mdp/minimization.mdp \
    --path_ff ff_files/ \
    --output_dir my_results \
    --name_molecule "GLYCAN" \
    --number_molecule 2 \
    --analyze_remove_pbc \
    --verbose \
    --keep_temp
```

### Example 5: Skip Analysis (only mapping and topology)

```bash
bayesian-potentials pipeline \
    --aa_tpr setup/md.tpr \
    --aa_xtc setup/md.xtc \
    --aa_gro setup/md.gro \
    --aa_itp setup/carb.itp \
    --beads_json json/beads_config.json \
    --force_application "random=[1250,30;30,4]" \
    --beads_position com \
    --input_mdp mdp/minimization.mdp \
    --path_ff ff_files/ \
    --skip_analysis
```

## Output Files

After successful pipeline execution, you'll find:

```
results/
├── GMX/                           # GROMACS files
│   ├── mapped.xtc                 # Coarse-grained trajectory
│   ├── cg.gro                     # CG coordinates
│   ├── cg.itp                     # CG molecule topology
│   ├── topol_cg.top               # Complete GROMACS topology
│   └── CG.tpr                     # GROMACS run input file
├── CG_MARTINI3/                   # cg-martini3 outputs
│   ├── GRO/                       # CG structure files
│   ├── ITP/                       # CG topology files
│   ├── JSON/                      # Beads definition files
│   └── NDX/                       # Index files (cg.ndx, bonds.ndx, etc.)
├── XVG/                           # Analysis results
│   ├── bonds/                     # Bond analysis
│   ├── angles/                    # Angle analysis
│   └── dihedrals/                 # Dihedral analysis
├── NDX/                           # Index files copy
│   ├── cg.ndx
│   ├── bonds.ndx
│   ├── angles.ndx
│   ├── dihedrals.ndx
│   └── cg.map
└── SUMMARY.txt                    # Execution summary
```

### Analysis Output Format

**bonds/ directory:**
- `bond_0.xvg` - Distance vs time
- `distr_bond_0.xvg` - Distance distribution
- `data_bonds.txt` - Average and standard deviation
- `report_bonds.txt` - Processed bonds list

**angles/ directory:**
- `ang_0.xvg` - Angle vs time
- `distr_ang_0.xvg` - Angle distribution
- `data_angles.txt` - Average and standard deviation

**dihedrals/ directory:**
- `dih_0.xvg` - Dihedral vs time
- `distr_dih_0.xvg` - Dihedral distribution
- `data_dihedrals.txt` - Average and standard deviation

## Troubleshooting

### Common Issues

**1. "cg-martini3 not found"**
```bash
# Install cg-martini3
pip install carb_param_automated_martini3

# Or install from source
git clone https://github.com/anacletosouza/carb_param_automated_martini3.git
cd carb_param_automated_martini3
pip install -e .
```

**2. "Missing force field files"**
```bash
# The package includes default ff_files in data directory
# They should be auto-detected. If not, specify path:
--path_ff /path/to/ff_files
```

**3. "grompp fails with warnings"**
```bash
# Increase max warnings
--maxwarn 10

# Or check your MDP parameters
gmx check -f minimization.mdp
```

**4. "PBC issues in trajectory"**
```bash
# Force PBC removal during mapping
--remove_pbc

# Or during analysis
--analyze_remove_pbc
```

**5. "No bonds/angles/dihedrals found"**
```bash
# Check that index files were generated correctly
ls -la results/NDX/
# Should show: bonds.ndx, angles.ndx, dihedrals.ndx

# Check connectivity report
cat results/NDX/connectivity_report.txt
```

**6. "Module not found"**
```bash
# Install in development mode
pip install -e .

# Check installation
python -c "import bayesian_potentials; print(bayesian_potentials.__file__)"
```

## License

MIT License - see LICENSE file for details.

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{bayesian_potentials_2026,
  author = {Anacleto Silva de Souza},
  title = {Bayesian Potentials: AA to CG Mapping Pipeline with Analysis Tools},
  year = {2026},
  url = {https://github.com/anacletosouza/bayesian-potentials}
}
```

