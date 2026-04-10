```markdown
# CG Trajectory Mapper

A Python tool to map atomistic molecular dynamics trajectories to coarse-grained (CG) resolution using GROMACS.

## Method

The tool maps atomistic trajectories to CG resolution using **center-of-mass** (COM) mapping:

$$\mathbf{R}_{\text{bead}} = \frac{\sum m_i \mathbf{r}_i}{\sum m_i}$$

where:
- $\mathbf{R}_{\text{bead}}$ is the position of the CG bead
- $m_i$ and $\mathbf{r}_i$ are the mass and position of atom $i$
- The sum runs over all atoms assigned to the bead

The tool wraps two GROMACS commands:
1. `gmx trjconv -pbc whole` - fixes periodic boundary conditions
2. `gmx traj -com` - computes COM positions for each bead group

## Installation

### Requirements
- Python 3.6+
- GROMACS (2019.x or later)
- NumPy (optional, for analysis)

### Download
```bash
git clone https://github.com/yourusername/cg-mapper.git
cd cg-mapper
chmod +x map_trajectory_to_cg.py
```

## Input Files

| File | Format | Description |
|------|--------|-------------|
| **Index file** (`.ndx`) | GROMACS index | Defines which atoms belong to each CG bead |
| **TPR file** (`.tpr`) | GROMACS binary | Run input with topology and coordinates |
| **XTC file** (`.xtc`) | GROMACS trajectory | Atomistic trajectory to map |

### Index File Format Example
```bash
[ bead1 ]
1 2 3 4 5 6
[ bead2 ]
7 8 9 10 11 12
[ bead3 ]
13 14 15 16
```

## Usage

### Basic Usage

Map a full trajectory to CG resolution:
```bash
python map_trajectory_to_cg.py -i cg.ndx -t setup/md.tpr -x setup/md.xtc
```

### Common Examples

#### 1. Map trajectory with custom output name
```bash
python map_trajectory_to_cg.py -i cg.ndx -t setup/md.tpr -x setup/md.xtc -o mapped.xtc
```

#### 2. Extract first frame as CG structure (.gro)
```bash
python map_trajectory_to_cg.py -i cg.ndx -t setup/md.tpr -x setup/md.xtc --output_cg_gro molecule.gro
```

#### 3. Generate CG .tpr for simulation
```bash
python map_trajectory_to_cg.py \
    -i cg.ndx \
    -t setup/md.tpr \
    -x setup/md.xtc \
    --output_cg_gro molecule.gro \
    --output_cg_tpr CG.tpr \
    --cg_top system_CG.top \
    --input_mdp martini.mdp
```

#### 4. Verbose mode (see GROMACS commands)
```bash
python map_trajectory_to_cg.py -i cg.ndx -t setup/md.tpr -x setup/md.xtc --verbose
```

#### 5. Dry run (preview commands without executing)
```bash
python map_trajectory_to_cg.py -i cg.ndx -t setup/md.tpr -x setup/md.xtc --dry-run
```

#### 6. Skip PBC correction
```bash
python map_trajectory_to_cg.py -i cg.ndx -t setup/md.tpr -x setup/md.xtc --no-pbc
```

#### 7. Use MPI version of GROMACS
```bash
python map_trajectory_to_cg.py -i cg.ndx -t setup/md.tpr -x setup/md.xtc --gmx-cmd gmx_mpi
```

## Command Line Options

### Required Arguments

| Argument | Short | Description |
|----------|-------|-------------|
| `--index_cg` | `-i` | GROMACS index file with bead definitions |
| `--aa_tpr` | `-t` | Atomistic GROMACS .tpr file |
| `--aa_xtc` | `-x` | Atomistic trajectory (.xtc) file |

### Optional Arguments

| Argument | Short | Default | Description |
|----------|-------|---------|-------------|
| `--output_mapped` | `-o` | `mapped.xtc` | Output CG trajectory file |
| `--output_cg_gro` | - | `None` | Extract first frame as CG .gro |
| `--output_cg_tpr` | - | `None` | Generate CG .tpr for simulation |
| `--cg_top` | - | `None` | CG topology file (.top) for grompp |
| `--input_mdp` | - | `None` | GROMACS parameter file (.mdp) |
| `--maxwarn` | - | `1` | Max warnings for grompp |
| `--pbc` | `-p` | `True` | Apply PBC correction |
| `--no-pbc` | - | - | Skip PBC correction |
| `--verbose` | `-v` | `False` | Print detailed execution info |
| `--dry-run` | - | `False` | Show commands without executing |
| `--keep-temp` | - | `False` | Preserve temporary files |
| `--gmx-cmd` | - | `gmx` | GROMACS command |

## Workflow Example

Complete workflow for parameterizing a new molecule:

```bash
# 1. Generate atomistic reference data
cd AA_REF/MDRUN
gmx mdrun -deffnm md -cpi md.cpt

# 2. Create CG mapping using CGbuilder
# http://www.cgmartini.nl/index.php/tools/cgbuilder

# 3. Map trajectory to CG resolution
python map_trajectory_to_cg.py \
    -i ndx/cg.ndx \
    -t setup/md.tpr \
    -x setup/md_center_fit.xtc \
    -o mapped.xtc \
    --output_cg_gro molecule.gro \
    --verbose

# 4. Generate CG .tpr for simulation
python map_trajectory_to_cg.py \
    -i ndx/cg.ndx \
    -t setup/md.tpr \
    -x setup/md.xtc \
    --output_cg_tpr CG.tpr \
    --cg_top system_CG.top \
    --input_mdp martini.mdp \
    --maxwarn 2

# 5. To generate CG.tpr, molecule_CG.gro, mapped.xtc    
python map_trajectory_to_cg.py \
    -i ndx/cg.ndx \
    -t setup/md.tpr \
    -x setup/md_center_fit.xtc \
    -o setup/mapped.xtc \
    --output_cg_gro setup/molecule.gro \
    --output_cg_tpr setup/CG.tpr \
    --cg_top setup/system_CG.top \
    --input_mdp mdp/minimization.mdp \
    --maxwarn 2 \
    --verbose

# 5. Generate bonded parameter distributions
# Basic analysis (without PBC removal)
bayesian-potentials analyze \
    --bonds_ndx OUTPUT/NDX/bonds.ndx \
    --angles_ndx OUTPUT/NDX/angles.ndx \
    --dihedrals_ndx OUTPUT/NDX/dihedrals.ndx \
    --xtc_file mapped.xtc \
    --tpr_file CG.tpr

# With PBC removal (default groups are "System" for both)
bayesian-potentials analyze \
    --bonds_ndx OUTPUT/NDX/bonds.ndx \
    --angles_ndx OUTPUT/NDX/angles.ndx \
    --dihedrals_ndx OUTPUT/NDX/dihedrals.ndx \
    --xtc_file mapped.xtc \
    --tpr_file CG.tpr \
    --remove_pbc

# With custom groups and keeping intermediate files
bayesian-potentials analyze \
    --bonds_ndx OUTPUT/NDX/bonds.ndx \
    --angles_ndx OUTPUT/NDX/angles.ndx \
    --dihedrals_ndx OUTPUT/NDX/dihedrals.ndx \
    --xtc_file mapped.xtc \
    --tpr_file CG.tpr \
    --remove_pbc \
    --group_1 "Backbone" \
    --group_2 "System" \
    --keep_intermediate

```

## Output Verification

Verify the generated files:

```bash
# Check CG trajectory integrity
gmx check -f mapped.xtc

# Extract first frame from CG trajectory
echo 0 | gmx trjconv -f mapped.xtc -s setup/md.tpr -o first_frame.gro -dump 0

# Count number of beads in output
grep -v "^;" first_frame.gro | tail -n +2 | head -1 | awk '{print (NF-2)}'

# Check CG .tpr
gmx check -t CG.tpr
```

## Example Output

```
Found 7 bead groups in cg.ndx

>>> Correcting periodic boundary conditions
Command: echo 0 | gmx trjconv -f setup/md.xtc -s setup/md.tpr -o /tmp/gmx_map_xxxxx/whole.xtc -pbc whole

>>> Mapping 7 beads to CG trajectory
Command: echo 0 1 2 3 4 5 6 | gmx traj -f /tmp/gmx_map_xxxxx/whole.xtc -s setup/md.tpr -oxt mapped.xtc -n cg.ndx -ng 7 -com

✓ Success! CG trajectory saved to: mapped.xtc

>>> Generating CG .gro file: molecule.gro
Command: echo 0 1 2 3 4 5 6 | gmx traj -f mapped.xtc -s setup/md.tpr -oxt molecule.gro -n cg.ndx -ng 7 -com -b 0 -e 0
✓ CG .gro saved to: molecule.gro

==================================================
SUMMARY
==================================================
Input files:
  - Index file:      cg.ndx
  - AA .tpr:         setup/md.tpr
  - AA .xtc:         setup/md.xtc

Parameters:
  - Number of beads: 7
  - PBC correction:  Yes

Output files:
  - CG trajectory: mapped.xtc (45.23 MB)
  - CG .gro: molecule.gro (0.15 MB)
```

## Error Handling

The tool provides clear error messages when inputs are missing:

```bash
# Missing required arguments for .tpr generation
python map_trajectory_to_cg.py --output_cg_tpr CG.tpr

ERROR: --output_cg_tpr requires the following arguments: --cg_top, --input_mdp

Example:
  map_trajectory_to_cg.py --output_cg_tpr CG.tpr --cg_top system_CG.top --input_mdp martini.mdp
```

## Notes

- **Center-of-mass mapping** is used (requires proper masses in TPR file)
- For Martini 3 **center-of-geometry** mapping, modify masses to 1.0 in topology first
- The index file must contain only bead group definitions
- PBC correction is applied by default for better quality mapping
- Temporary files are automatically cleaned up unless `--keep-temp` is used

## Citation

If you use this tool, please cite:
- GROMACS: Abraham et al., SoftwareX, 2015, 1-2, 19-25
- Martini 3: Souza et al., Nat. Methods, 2021, 18, 382-388

## License

MIT License

## Contact

For issues or questions, please open an issue on GitHub.
```




