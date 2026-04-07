# carb_param_automated_martini3

[![PyPI version](https://badge.fury.io/py/carb-param-automated-martini3.svg)](https://badge.fury.io/py/carb-param-automated-martini3)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Automated coarse-grained topology generation for glycans using Martini 3 force field.

## Installation
```
conda create -n cg_martini3 python=3.10 -y
conda activate cg_martini3
```
```
pip install git+ssh://git@github.com/anacletosouza/carb_param_automated_martini3.git
```
```
cg-martini3 -h (it contain all instrutions for usage)
```

### Option 1: Install via pip (recommended)

```bash
pip install git+https://github.com/anacletosouza/carb_param_automated_martini3.git
```

# Complete Procedure for CG Topology Generation

## Overview

This pipeline generates a complete coarse-grained (CG) topology for glycans by mapping all-atom (AA) molecular dynamics data onto a CG representation. The process combines structural information from CG coordinate files, index-based atom grouping, all-atom parameter files, and Martini 3 force field definitions to produce GROMACS-compatible topology files including bonds, angles, and dihedrals.

## Pipeline Workflow

```
INPUT FILES                    STEP 1-3                    OUTPUT FILES
├── carb.gro          ──►     Generate CG               ├── JSON/
├── beads_config.json          Structure & Index        │   └── beads_definition_com.json
└── carb.itp                   Files                    ├── GRO/
                                                        ├── carb.gro
                                                 │      │   └── cg.gro
                                                     ▼                   ├── NDX/
                              STEP 4-7                  │   ├── cg.ndx
                              Generate Topology         │   ├── cg.map
                              & Geometric               │   ├── bonds.ndx
                              Parameters                │   ├── angles.ndx
                                                        │   └── dihedrals.ndx
                                                        └── ITP/
                                                            └── carb_cg.itp
                                                            └── cg_final.itp
```

## Quick Start

### Basic Usage

## Individual Scripts Usage

After installation, you can use each script independently:

### 1. Convert GRO to JSON bead definitions

```bash
gro-to-json --input_gro carb.gro --output_json beads.json \
            --config_json beads_config.json --position_beads com
```

### 2. Convert JSON to coarse-grained GRO

```bash
json-to-gro --json_beads_definition beads.json --input_gro carb.gro \
            --output_beads_gro cg.gro --coordinate_option com
```

### 3. Generate index and map files

```bash
generate-index --json_beads_definition beads.json --input_gro carb.gro \
               --output_ndx cg.ndx --output_map cg.map
```

### 4. Generate ITP topology

```bash
generate-itp --cg_gro cg.gro --cg_ndx cg.ndx --aa_itp carb.itp \
             --def_json martini3_map.json --output carb_cg.itp
```

### 5. Generate bonds

```bash
generate-bonds --cg_gro cg.gro --output_ndx bonds.ndx \
               --max_distance 5.0 --tolerance 0.5
```

### 6. Generate angles

```bash
generate-angles --gro cg.gro --bonds bonds.ndx --output angles.ndx
```

### 7. Generate dihedrals

```bash
generate-dihedrals -g cg.gro -b bonds.ndx -o dihedrals.ndx
```

### 8. Generate final ITP

```bash
generate-final --input_bonds bonds.ndx --input_angles angles.ndx \
               --input_dihedrals dihedrals.ndx --itp_incomplete carb_cg.itp \
               --cg_gro cg.gro --itp_complete carb_final.itp \
               --exclusion all --force_application "fix=[1250;25]"
```

## Input Files Description

### 1. carb.gro - All-Atom Structure File
**Format:** GROMACS GRO file with atomic coordinates

```
Glycan AA structure              # Title line
   123                           # Number of atoms
    1ASN     N    1   6.339   0.630  25.238
    1ASN    HN    2   6.224   0.872  25.367
    1ASN    CA    3   5.946   1.068  25.447
    ...
```


### 2. beads_config.json - Beads Definition File
**Format:** JSON file defining CG bead composition

```json
{
  "protein": {
    "BB": ["N", "HN", "CA", "HA", "C", "O"]
  },
  "link_residues": {
    "ASN": {
      "SC1": ["CB", "HB1", "HB2", "CG", "OD1", "ND2", "HD22"]
    },
    "SER": {
      "SC1": ["CB", "HB1", "HB2", "OG", "HG"]
    },
    "THR": {
      "SC1": ["CB", "HB", "OG1", "HG1", "CG2", "HG21", "HG22", "HG23"]
    }
  },
  "carbohydrates": {
    "BG": {
      "1": ["C1", "H1", "O5", "C5", "H5"],
      "2": ["C2", "H2", "N", "HN", "C", "O", "CT", "HT1", "HT2", "HT3"],
      "3": ["C3", "H3", "O3", "HO3"],
      "4": ["C4", "H4", "C6", "H61", "H62", "O6"]
    },
    "AF": {
      "1": ["C1", "H1", "O5", "C5", "H5"],
      "2": ["C2", "H2", "O2", "HO2", "C3", "H3", "O3", "HO3"],
      "3": ["C4", "H4", "O4", "HO4", "C6", "H61", "H62", "H63"]
    },
    "BM": {
      "1": ["C1", "H1", "O5", "C5", "H5"],
      "2": ["C2", "H2", "O2", "HO2", "C3", "H3", "O3"],
      "3": ["C4", "H4", "O4", "HO4", "C6", "H61", "H62", "O6"]
    },
    "AM": {
      "1": ["C1", "H1", "O5", "C5", "H5"],
      "2": ["C2", "H2", "O2", "C3", "H3", "O3", "HO3"],
      "3": ["C4", "H4", "O4", "HO4", "C6", "H61", "H62", "O6", "HO6"]
    }
  }
}

```

**Structure:**
- Key: Bead identifier (must match cg.gro column 2)
- Value: [Residue name, [List of atom names to include]]

### 3. carb.itp - All-Atom Topology File
**Format:** GROMACS ITP file with atom parameters (CHARMM Force Field)

```
[ atoms ]
1      NH1      17     ASN      N      1      -0.470000  14.0070   
2      H        17     ASN      HN     2      0.310000   1.0080    
3      CT1      17     ASN      CA     3      0.070000   12.0110   
...
```

**Purpose:** Provides reference masses and charges for AA atoms

### 4. definitions_atoms_ff_martini3.json - Martini 3 Mapping
**Format:** JSON mapping CG bead types to Martini 3 particle types

```json
{
    "BB": "SN4a",
    "SC1": "P3",
    "BG": "TP5",
    "AF": "TP5",
    "BM": "TP5",
    "AM": "TP5"
}
```

## Pipeline Steps

## Script 1: `1-gro_to_beads_json.py`
**Purpose:** Converts an all-atom GRO file into a JSON bead definition file for coarse-graining.

**What it does:**
- Parses a GRO file with atomic coordinates
- Groups atoms by residue
- Identifies protein residues and carbohydrate monomers based on a configuration file
- Assigns atoms to coarse-grained beads (e.g., BB for backbone, SC1 for sidechain, BG1/BG2 for carbohydrates)
- Computes bead positions (center of mass or geometric center)
- Outputs a JSON file with bead definitions (atoms per bead + positions)

---

## Script 2: `2-json_to_beads_gro.py`
**Purpose:** Converts a JSON bead definition back into a GROMACS GRO file with CG beads.

**What it does:**
- Reads a JSON bead definition file (from Script 1)
- Parses the original GRO file to get atomic coordinates
- For each bead, identifies constituent atoms and computes bead position (COM or geometric center)
- Chooses a representative atom name for each bead (prioritizing C4 > C3 > C2 > C1 > C > O)
- Writes a new GRO file containing only CG beads

---

## Script 3: `3-index_map.py`
**Purpose:** Generates index (NDX) and mapping (MAP) files for GROMACS simulations.

**What it does:**
- Reads the JSON bead definition and original GRO file
- Creates atom index groups for each CG bead in NDX format (`[ bead_name ]` followed by atom indices)
- Generates a MAP file with:
  - `[ martini ]` section listing all bead names
  - `[ atoms ]` section mapping each atom type to its primary bead and all beads containing it
- Orders beads with protein residues first, then carbohydrates

---

## Script 4: `4-defining_atoms_type.py`
**Purpose:** Builds a coarse-grained topology (ITP file) for GROMACS.

**What it does:**
- Parses CG GRO, CG NDX, all-atom ITP, and a bead type mapping JSON
- Maps each CG bead to a Martini 3 bead type (e.g., BB → SN4a, SC1 → P3)
- Computes bead mass and charge either by:
  - Summing constituent atom properties (default), or
  - Using default Martini 3 values (72 amu, 0 charge)
- Outputs a GROMACS ITP file with `[ moleculetype ]` and `[ atoms ]` sections
- Maintains sequential residue numbering

---

## Script 5: `5-bonds_general.py`
**Purpose:** Detects internal and external bonds between CG beads.

**What it does:**
- Reads a CG GRO file with bead coordinates
- **Internal bonds:** Connects beads within the same residue based on proximity (max 4 bonds per bead)
- **External bonds:** Connects beads between different residues (glycosidic bonds) with rules:
  - Branching residues (4 beads): only C1 and C4 can bond externally
  - Non-branching residues (3 beads): C1, C3, C4 can bond externally
  - C1 can make at most one external bond
  - Residues 2 and 5 can make two external bonds (special cases)
  - Distance cutoff with tolerance
- Outputs a bonds NDX file for GROMACS

---

## Script 6: `6-angles.py`
**Purpose:** Calculates internal and external angles between CG beads.

**What it does:**
- Reads CG GRO file and bonds NDX file
- Builds adjacency lists from bond connections
- **Internal angles:** Finds i-j-k triplets where all three beads are in the same residue
- **External angles:** Finds angles spanning residue boundaries (e.g., neighbor - connecting bond - other residue)
- Calculates angles in degrees using the dot product formula
- Outputs an angles NDX file with annotations and angle values

---

## Script 7: `7-dihedrals.py`
**Purpose:** Calculates dihedral angles from CG bead connections.

**What it does:**
- Reads CG GRO file and bonds NDX file
- Separates bonds into internal (same residue) and external (different residues)
- Finds all possible i-j-k-l quadruplets where consecutive pairs are bonded
- **Internal dihedrals:** All four beads within the same residue
- **External dihedrals:** Spans two residues connected by an external bond (i from residue A, j-k external bond, l from residue B)
- Removes duplicate dihedrals (considering reverse direction as identical)
- Calculates dihedral angles using vector cross product method (with sign)
- Outputs a dihedrals NDX file for GROMACS

---

## Script 8: `8-itp_final.py`
**Purpose:** Generates the final complete ITP file with bonds, angles, and dihedrals.

**What it does:**
- Reads bond, angle, and dihedral definitions from NDX files
- Merges them with the incomplete ITP topology
- Applies force constants (fixed or random) to bonds and angles
- Configures exclusions for non-bonded interactions
- Outputs a complete, ready-to-use GROMACS ITP file

---

**What it does:**
- Reads CG GRO file and bonds NDX file
- Separates bonds into internal (same residue) and external (different residues)
- Finds all possible i-j-k-l quadruplets where consecutive pairs are bonded
- **Internal dihedrals:** All four beads within the same residue
- **External dihedrals:** Spans two residues connected by an external bond (i from residue A, j-k external bond, l from residue B)
- Removes duplicate dihedrals (considering reverse direction as identical)
- Calculates dihedral angles using vector cross product method (with sign)
- Outputs a dihedrals NDX file for GROMACS

---

## Workflow Summary (Typical Usage)

**Script 1** → Convert AA GRO to JSON bead definitions
**Script 2** → Convert JSON bead definitions to CG GRO
**Script 3** → Generate NDX and MAP files
**Script 4** → Build CG topology (ITP)
**Script 5** → Detect bonds between beads
**Script 6** → Calculate angles from bonds
**Script 7** → Calculate dihedrals from bonds
**Script 8** → Calculate final topology


## Output Files Detailed Description

### JSON Directory
- **`beads_definition_com.json`**: CG bead definitions with calculated positions

### GRO Directory
- **`carb.gro`**: Copy of original all-atom structure
- **`cg.gro`**: Coarse-grained bead coordinates

### NDX Directory
- **`cg.ndx`**: GROMACS index file mapping AA atoms to CG beads
- **`cg.map`**: Simple mapping file (AA index → CG bead name)
- **`bonds.ndx`**: Bond definitions for GROMACS
- **`angles.ndx`**: Angle definitions for GROMACS
- **`dihedrals.ndx`**: Dihedral definitions for GROMACS

### ITP Directory
- **`carb_cg.itp`**: Complete CG topology for GROMACS

## Bond Detection Logic

The bond detection algorithm (`5-bonds_general.py`) uses:

1. **Distance-based bonding**: Beads within `max_distance ± tolerance` are considered bonded
2. **Connectivity analysis**: Identifies bonded networks
3. **Report generation**: Lists all detected bonds with distances

**Example output:**
```
Bond: 1BB (1) - 1SC1 (2)  distance: 3.45 Å
Bond: 1SC1 (2) - BG1 (3)   distance: 4.12 Å
```

## Angle Generation Details

### Internal Angles (Within Residues)

Internal angles are defined when a single residue contains 3 or more beads that form connected bonds. The angle is calculated using three consecutive indices (i, j, k) where:

j is the central/common index (vertex)

i and k are two different neighbors of j connected by bonds

All three indices belong to the same residue number

Example from residue 2 (indices 3,4,5,6 with bonds: 3-6, 4-5, 5-3, 6-5):

Angle 3-6-5: j=6 (common), i=3 and k=5 are both bonded to 6

Angle 5-3-6: j=3 (common), i=5 and k=6 bonded to 3

Output format:

```
3     6     5     ;  2_C1--2_C4--2_C3  angle = 105.4 degree
```

### External Angles (between different residues)

External angles are defined when two different residues are linked by a bond. The angle is calculated in two possible ways:

Type 1: i (residue A) - j (residue A/B bond point) - k (residue B)

j is the bond index connecting residues

i is another bead from the same residue as j (connected to j by a bond)

k is a bead from the other residue (connected to j by the existing bond)

Type 2: i (residue A) - j (residue A) - k (residue B)

i and j are from same residue (connected by bond)

k is from different residue (connected to j by external bond)

Example: Bond 2-3 connects residue 1 (index 2) and residue 2 (index 3)

Angle 2-3-6: j=3, i=2 (from residue 1), k=6 (neighbor of 3 in residue 2)

Angle 1-2-3: j=2, i=1 (neighbor of 2 in residue 1), k=3 (from residue 2)

Output format:

```
2     3     6     ;  1_CB--2_C1--2_C4  angle = 95.2 degree
```

- Annotations showing {residue numbers}_{atom representation of the bead}

- Calculated angle values in degrees

## Dihedral Generation Details

The algorithm:

For each external bond (a-b connecting residues A and B)

Find all internal connections to 'a' (from residue A)

Find all internal connections to 'b' (from residue B)

Combine them to form dihedrals: i - j - k - l

```
        i
         \
          j --- k
               /
              l
```

Automatic detection: No manual dihedral definition needed

Handles branching: Works with complex carbohydrate topologies

Duplicate removal: Prevents redundant calculations

Sign-aware: Preserves chirality information (positive/negative angles)

- Each dihedral is labeled with residue numbers and bead names

## Complete Command Reference

### Main Pipeline Script

Args to cg-martini3:

### Main Pipeline Script

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--input_gro` | Yes | - | Input all-atom GRO file |
| `--input_beads_definitions` | Yes | - | JSON configuration file for bead definitions |
| `--output_dir` | Yes | - | Output directory for all generated files |
| `--beads_position` | Yes | - | Position calculation method: 'com' (center of mass) or 'geom' (geometric center) |
| `--aa_itp` | Yes | - | All-atom ITP file with [ atoms ] section |
| `--definition_martini3_bead_type_to_itp` | No | Auto-detected | Path to Martini3 mapping file (definitions_atoms_ff_martini3.json) |
| `--max_distance` | No | 5.0 | Maximum distance for external bonds in Å |
| `--tolerance` | No | 0.5 | Tolerance for external bonds in Å |
| `--keep_intermediate` | No | False | Keep all intermediate files |
| `--exclusion` | No | "all" | Exclusions: 'all', 'none', or list like '1-5,7-9,20-26' |
| `--force_application` | No | "fix=[1250;25]" | Force constants: 'fix=[bond_k;angle_k]' or 'random=[mean_bond,sd_bond;mean_angle,sd_angle]' |
| `--python_dir` | No | (not used) | Python directory (kept for compatibility, not used) |

### Viewing CG Structure

```bash
# Convert GRO to PDB for visualization
gmx editconf -f cg.gro -o cg.pdb

# View in VMD or PyMOL
vmd cg.pdb
```

## References

1. **Martini 3 Force Field**: Souza, P. C. T., et al. (2021). "Martini 3: a general purpose force field for coarse-grained molecular dynamics." Nature Methods, 18(4), 382-388.

## License

This pipeline is distributed under the MIT License. See LICENSE file for details.

## Author

Anacleto Silva de Souza

## Version History

- **v1.0** (2026-03-31): Initial release with complete pipeline
- **v1.1** (2026-04-04): Added bond tolerance parameter, simplified input structure

## Citation

If you use this pipeline in your research, please cite:

```
Souza, A. S. Marrink, S.J. (2026). CG Topology Generation Pipeline for Glycans. 
```

---

**For questions, issues, or contributions, please contact the author or open an issue on the repository.**
```

