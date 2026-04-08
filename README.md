```markdown
# Bayesian Potentials - Coarse-Grained Mapping Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)

A complete pipeline to map atomistic (AA) molecular dynamics trajectories to coarse-grained (CG) representations, generate CG topologies, and prepare GROMACS input files.

##   Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Required Files](#required-files)
- [Commands](#commands)
  - [`bayesian-potentials map`](#bayesian-potentials-map)
  - [`bayesian-potentials gen-top`](#bayesian-potentials-gen-top)
  - [`bayesian-potentials pipeline`](#bayesian-potentials-pipeline)
- [Complete Examples](#complete-examples)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)

## Overview

This tool converts atomistic MD simulations to coarse-grained representations (Martini 3 compatible). The pipeline:

1. **Maps** AA trajectories to CG beads based on index definitions
2. **Generates** CG topology (ITP) with proper masses and charges
3. **Creates** complete GROMACS topology and runs `grompp` for simulation readiness

## Installation

```bash
# Clone the repository
git clone https://git@github.com/anacletosouza/bayesian-potentials.git
cd bayesian-potentials

# Install the package
pip install -e .

# Or install directly
pip install bayesian-potentials

In addition, the package utilize modules of the cg-martini3. You need to install cg-martini3:

# Clone the repository
git clone https://git@github.com/anacletosouza/carb_param_automated_martini3.git
cd carb_param_automated_martini3

# Install the package
pip install -e .

# Or install directly
pip install carb_param_automated_martini3

or

pip install git+ssh://git@github.com/anacletosouza/carb_param_automated_martini3.git
```

## Required Files

Before running the pipeline, prepare these files:

| File | Description | Example |
|------|-------------|---------|
| `cg.ndx` | CG index file mapping AA atoms to CG beads | See as follows in the [example](#index-file-format) |
| `md.tpr` | AA GROMACS input file | `gmx grompp -f md.mdp -c conf.gro -p topol.top -o md.tpr` |
| `md.xtc` | AA trajectory (XTC/TRR) | Output from MD simulation |
| `carb.itp` | AA molecule ITP | Original force field ITP |
| `minimization.mdp` | GROMACS MDP parameters | Included in package |

### Index File Format (`cg.ndx`)

It is necessary to simulate your system, which will have the files:

- md.tpr
- md.xtc (without removing PBC condition)
- minimization.mdp (MDP file)

Furthermore, it is necessary

You can obtain this cg.ndx from cgbuilder webserver or running:

```
bash
cg-martini3 --input_gro carb.gro \                          #AA model containing all information about your molecule
            --input_beads_definitions beads_config.json \   #define the beads from AA model
            --output_dir OUTPUT_GEOM \                      #define the folder to output files in each GRO, ITP, NDX, JSON
            --beads_position geom \                         #or com
            --aa_itp carb.itp                               #topology from AA model
            --force_application "random=[1250,30;30,4]"     # or you can use "fix=[1250,25]"
```

In OUTPUT_GEOM/NDX, you will obtain cg.ndx file, which contains the bead mapping from AA model.

```
[ BB ]
1 2 3 4 12 13

[ SC1 ]
5 6 7 8 9 10 11

[ BG1 ]
14 15 16 17 18

[ BG2 ]
19 20 21 22 23 24 25 26 27 28

...
```

Each `[ Group_Name ]` defines a CG bead with the atom indices that belong to it.

## Commands

### `bayesian-potentials map`

Maps an atomistic trajectory to coarse-grained coordinates.

```bash
bayesian-potentials map \
                        --index_cg       ../ndx/cg.ndx   \
                        --aa_tpr         ../setup/md.tpr \
                        --aa_xtc         ../setup/md.xtc \
                        --output_mapped  mapped.xtc      \
                        --output_cg_gro  cg.gro          \
                        --remove_pbc                     \
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

#### Example Output
```
✓ Loaded index: 30 CG beads defined
✓ Mapped 1000 frames
✓ Saved: mapped.xtc
✓ Saved: cg.gro
```

### `bayesian-potentials gen-top`

Generates a CG topology file (ITP) for GROMACS.

```bash
bayesian-potentials gen-top \
                            --path_ff           ../ff_files                    \
                            --ff                martini_v3.0.0.itp             \
                            --ions              martini_v3.0.0_ions_v1.itp     \
                            --solvent           martini_v3.0.0_solvents_v1.itp \
                            --itp_ligand        cg.itp                         \
                            --name_molecule     "molecule"                     \
                            --number_molecule   1                              \
                            --output_topol      topol.top
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

#### Example Output
```
✓ Generated: topol.top
  #include "martini_v3.0.0.itp"
  #include "martini_v3.0.0_solvents_v1.itp"
  #include "cg.itp"
  
  [ system ]
  Molecule in aqueous solution
  
  [ molecules ]
  molecule     1
```

### `bayesian-potentials pipeline`

**Complete pipeline** - runs mapping, ITP generation, topology creation, and `grompp` in one command.

```bash
bayesian-potentials pipeline \
                             --cg_ndx             ../ndx/cg.ndx \
                             --aa_tpr             ../setup/md.tpr \
                             --aa_xtc             ../setup/md.xtc \
                             --output_mapped      mapped.xtc \
                             --output_cg_gro      cg.gro \
                             --aa_itp             ../setup/carb.itp \
                             --output_cg_itp      cg.itp \
                             --input_mdp          ../mdp/minimization.mdp \
                             --output_dir         ../results/ \
                             --output_topol       topol.top \
                             --path_ff            ../ff_files \
                             --name_molecule      "molecule" \
                             --number_molecule    1
```

#### Arguments

| Category | Argument | Default | Description |
|----------|----------|---------|-------------|
| **Required** | `--cg_ndx` | - | CG index file |
| | `--aa_tpr` | - | AA TPR file |
| | `--aa_xtc` | - | AA trajectory |
| | `--output_mapped` | - | Output mapped trajectory |
| | `--output_cg_gro` | - | Output CG GRO file |
| | `--aa_itp` | - | AA ITP file |
| | `--output_cg_itp` | - | Output CG ITP |
| **Optional** | `--input_mdp` | auto-detected | MDP for grompp |
| | `--path_ff` | auto-detected | Force field directory |
| | `--output_dir` | `results` | Output directory |
| | `--output_topol` | `topol_cg.top` | Topology filename |
| | `--name_molecule` | `molecule` | Molecule name |
| | `--number_molecule` | `1` | Number of molecules |
| | `--maxwarn` | `1` | Max grompp warnings |
| | `--remove_pbc` | `True` | Remove PBC |
| | `--skip_grompp` | `True` | Skip grompp step |
| | `--default_martini` | `False` | Use default Martini3 parameters |
| | `--verbose` | `False` | Verbose output |

## Complete Examples

### Example 1: Skip grompp (Generate only topology)

```bash
bayesian-potentials pipeline \
                            --cg_ndx indices/system.ndx \
                            --aa_tpr md.tpr \
                            --aa_xtc md.xtc \
                            --output_mapped mapped.xtc \
                            --output_cg_gro system.gro \
                            --aa_itp molecule.itp \
                            --output_cg_itp molecule_cg.itp \
                            --skip_grompp \
                            --output_dir topol_only/
```

### Example 4: Step-by-step (Manual mode)

```bash
# Step 1: Map trajectory
bayesian-potentials map \
                        --index_cg cg.ndx \
                        --aa_tpr md.tpr \
                        --aa_xtc md.xtc \
                        --output_mapped mapped.xtc \
                        --output_cg_gro cg.gro

# Step 2: Generate CG ITP
cg-generate-itp \
                 --cg_gro cg.gro \
                 --cg_ndx cg.ndx \
                 --aa_itp molecule.itp \
                 --output molecule_cg.itp \
                 --name_molecule molecule

# Step 3: Create topology
bayesian-potentials gen-top \
                            --path_ff ff_files/ \
                            --itp_ligand molecule_cg.itp \
                            --name_molecule MY_MOL \
                            --number_molecule 1 \
                            --output_topol topol.top

# Step 4: Run grompp manually
gmx grompp -f minimization.mdp -c cg.gro -p topol.top -o CG.tpr
```

## Output Files

After successful pipeline execution, you'll find:

```
results/
├── mapped.xtc          # Coarse-grained trajectory
├── cg.gro              # CG coordinates (last frame)
├── cg.itp              # CG molecule topology
├── topol.top           # Complete GROMACS topology
└── CG.tpr              # GROMACS run input file (if grompp ran)
```

### File Descriptions

| File | Purpose |
|------|---------|
| `mapped.xtc` | CG trajectory for analysis or simulation restart |
| `cg.gro` | CG coordinates for visualization or as input |
| `cg.itp` | Molecule definition with [atoms], [bonds], etc. |
| `topol.top` | Complete topology including force field, solvent, ions |
| `CG.tpr` | Binary run file ready for `gmx mdrun` |

## Troubleshooting

### Common Issues

**1. "Index file not found"**
```bash
# Use absolute paths or ensure relative paths are correct
--cg_ndx $(pwd)/ndx/cg.ndx
```

**2. "Missing force field files"**
```bash
# Set path explicitly
--path_ff /usr/local/share/martini/ff_files

# Or install default files
cp -r data/ff_files ~/.bayesian_potentials/
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
# Force PBC removal
--remove_pbc

# Or pre-process trajectory
echo "0" | gmx trjconv -f md.xtc -o md_noPBC.xtc -pbc mol -center
```

**5. "Module not found"**
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
```
@software{bayesian_potentials_2026,
  author = {Anacleto Silva de Souza},
  title = {Bayesian Potentials: AA to CG Mapping Pipeline},
  year = {2026},
  url = {https://github.com/anacletosouza/bayesian-potentials}
}
```

