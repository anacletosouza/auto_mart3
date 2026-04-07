#!/usr/bin/env python3
"""
General script to generate a Martini-style .itp file from statistics and bead data.
"""

import os
import json
import argparse
import pandas as pd
from datetime import datetime

# ============================================
# Helper functions
# ============================================

def read_tsv(file_path):
    """Read TSV file into pandas DataFrame."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} not found.")
    df = pd.read_csv(file_path, sep="\t")
    return df

def read_ndx(file_path):
    """Read .ndx index file and return list of tuples."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} not found.")
    
    indices = []
    with open(file_path, 'r') as f:
        current_group = None
        for line in f:
            line = line.strip()
            if line.startswith('['):
                current_group = line
            elif line and not line.startswith(';'):
                parts = line.split()
                if len(parts) == 2:
                    indices.append((int(parts[0]), int(parts[1])))
                elif len(parts) == 3:
                    indices.append(tuple(int(x) for x in parts))
                elif len(parts) == 4:
                    indices.append(tuple(int(x) for x in parts))
    return indices

def format_bond(i,j,length,force):
    return f"  {i:2d}  {j:2d}   1    {length:.3f}    {force:d}   ;"

def format_angle(i,j,k,angle,force):
    return f"  {i:2d}  {j:2d}  {k:2d}   2    {angle:.1f}    {force:d}   ;"

def format_dihedral(i,j,k,l,angle,force):
    return f"  {i:2d}  {j:2d}  {k:2d}  {l:2d}   2    {angle:.1f}    {force:d}   ;"

def get_bead_desc(beads, idx):
    """Return bead description string from ATOMS_JSON."""
    bead = next((b for b in beads if b['nr']==idx), None)
    if bead:
        return f"{bead['residue']} {bead['atom']}"
    else:
        return f"Bead {idx}"

# ============================================
# Argument parsing
# ============================================

parser = argparse.ArgumentParser(description="Generate a Martini .itp file from stats and ndx files.")
parser.add_argument("--bonds", required=True, help="TSV file with bond statistics")
parser.add_argument("--angles", required=True, help="TSV file with angle statistics")
parser.add_argument("--dihedrals", required=True, help="TSV file with dihedral statistics")
parser.add_argument("--atoms_json", required=True, help="JSON file with atoms/beads (type, mass, charge)")
parser.add_argument("--ndx_bounds", required=True, help="Ndx file for bonds")
parser.add_argument("--ndx_angles", required=True, help="Ndx file for angles")
parser.add_argument("--ndx_dihedrals", required=True, help="Ndx file for dihedrals")
parser.add_argument("--molecule_name", default="FA2", help="Molecule name")
parser.add_argument("--dihedrals_target", action='store_true', help="Use mean angle for dihedrals; default False uses 0")
parser.add_argument("--exclusion", action='store_true', help="Enable custom exclusions (requires --exclusion_json)")
parser.add_argument("--exclusion_json", help="JSON file with exclusion atom indices (required if --exclusion)")

args = parser.parse_args()

if args.exclusion and not args.exclusion_json:
    parser.error("--exclusion requires --exclusion_json")

# ============================================
# Load files
# ============================================

bonds_df = read_tsv(args.bonds)
angles_df = read_tsv(args.angles)
dihedrals_df = read_tsv(args.dihedrals)

with open(args.atoms_json, 'r') as f:
    atoms_data = json.load(f)

bond_pairs = read_ndx(args.ndx_bounds)
angle_triplets = read_ndx(args.ndx_angles)
dihedral_quartets = read_ndx(args.ndx_dihedrals)

exclusion_list = []
if args.exclusion:
    with open(args.exclusion_json, 'r') as f:
        exclusion_list = json.load(f)

# ============================================
# Force constants (default Martini 3)
# ============================================
FORCE_BOND = 1250
FORCE_ANGLE = 25
FORCE_DIHEDRAL = 25

# ============================================
# Generate ITP
# ============================================

output_file = f"{args.molecule_name}_final.itp"
timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

with open(output_file, 'w') as f:
    f.write(f";;;;;; {args.molecule_name} - Final topology\n")
    f.write(f"; Generated on {timestamp}\n\n")

    # Molecule type
    f.write("[ moleculetype ]\n")
    f.write(f"; molname    nrexcl\n  {args.molecule_name}     1\n\n")

    # Atoms
    f.write("[ atoms ]\n")
    f.write("; nr type resnr residue atom cgnr charge mass\n")
    for bead in atoms_data:
        nr = bead['nr']
        typ = bead['type']
        resnr = bead['resnr']
        residue = bead['residue']
        atom = bead['atom']
        cgnr = bead['cgnr']
        charge = bead['charge']
        mass = bead['mass']
        desc = f"{residue} {atom}"
        f.write(f"  {nr:2d}  {typ:4s}  {resnr:2d}  {residue:4s}  {atom:3s}  {cgnr:2d}  {charge:6.2f}  {mass:6.2f}   ; {desc}\n")
    f.write("\n")

    # Bonds
    f.write("[ bonds ]\n")
    f.write("; i  j  funct length force.k\n")
    for idx, row in bonds_df.iterrows():
        if idx < len(bond_pairs):
            i,j = bond_pairs[idx]
            f.write(format_bond(i,j,row['mean'],FORCE_BOND) + f" ; σ={row['sd']:.4f} {get_bead_desc(atoms_data,i)}-{get_bead_desc(atoms_data,j)}\n")

    f.write("\n")

    # Angles
    f.write("[ angles ]\n")
    f.write("; i  j  k  funct angle force.k\n")
    for idx, row in angles_df.iterrows():
        if idx < len(angle_triplets):
            i,j,k = angle_triplets[idx]
            f.write(format_angle(i,j,k,row['mean'],FORCE_ANGLE) + f" ; σ={row['sd']:.4f} {get_bead_desc(atoms_data,i)}-{get_bead_desc(atoms_data,j)}-{get_bead_desc(atoms_data,k)}\n")

    f.write("\n")

    # Dihedrals
    f.write("[ dihedrals ]\n")
    f.write("; i  j  k  l  funct angle force.k\n")
    for idx, row in dihedrals_df.iterrows():
        if idx < len(dihedral_quartets):
            i,j,k,l = dihedral_quartets[idx]
            angle_val = row['mean'] if args.dihedrals_target else 0
            f.write(format_dihedral(i,j,k,l,angle_val,FORCE_DIHEDRAL) + f" ; original mean={row['mean']:.1f} σ={row['sd']:.1f}\n")

    # Exclusions
    f.write("\n[ exclusions ]\n")
    if args.exclusion:
        f.write("; custom exclusions provided\n")
        for line in exclusion_list:
            f.write(" ".join(str(x) for x in line) + "\n")
    else:
        f.write("; full exclusions for 1-4 and 1-5\n")
        all_indices = [str(bead['nr']) for bead in atoms_data]
        f.write(" ".join(all_indices) + "\n")

print(f"\nITP file generated: {output_file}")
