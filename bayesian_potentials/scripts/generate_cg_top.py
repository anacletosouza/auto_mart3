#!/usr/bin/env python3
"""
Description: Generate a GROMACS topology (.top) file with specified force fields, ligands, molecules.

Usage:
python ../2-obtaining_cg_top.py \
    --path_ff ../ff_files/ \
    --ff martini_v3.0.0.itp \
    --ions martini_v3.0.0_ions_v1.itp \
    --solvent martini_v3.0.0_solvents_v1.itp \
    --itp_ligand ../setup/cg.itp \
    --name_molecule molecule \
    --number_molecule 2 \
    --title_comments "Topology system in Martini 3" \
    --title_system "molecule in aqueous solution" \
    --output_topol topol_cg.top
"""

import argparse
import os

def clean_path(base_path, filename):
    """Join base path and filename, handling trailing slashes properly."""
    base_path = base_path.rstrip('/')
    filename = filename.lstrip('/')
    return f"{base_path}/{filename}"

def main():
    parser = argparse.ArgumentParser(
        description="Generate a GROMACS topology file (.top)."
    )

    # Required arguments
    parser.add_argument("--path_ff", type=str, required=True,
                        help="Directory containing force field ITP files.")
    
    parser.add_argument("--ff", type=str, default="martini_v3.0.0.itp",
                        help="Force field ITP filename (default: martini_v3.0.0.itp)")
    
    parser.add_argument("--ions", type=str, default="martini_v3.0.0_ions_v1.itp",
                        help="Ions ITP filename (default: martini_v3.0.0_ions_v1.itp)")
    
    parser.add_argument("--solvent", type=str, default="martini_v3.0.0_solvents_v1.itp",
                        help="Solvent ITP filename (default: martini_v3.0.0_solvents_v1.itp)")
    
    parser.add_argument("--itp_ligand", type=str, required=True,
                        help="Ligand ITP file to include (use None to skip)")
    
    parser.add_argument("--name_molecule", type=str, required=True,
                        help="Molecule name for the ligand")
    
    parser.add_argument("--number_molecule", type=int, required=True,
                        help="Number of molecules for the ligand")

    # Optional arguments
    parser.add_argument("--title_comments", type=str, default="",
                        help="Optional comment title at the top of the file.")
    parser.add_argument("--title_system", type=str, default="",
                        help="Optional system title.")
    parser.add_argument("--output_topol", type=str, required=True,
                        help="Output topology filename.")

    args = parser.parse_args()

    # Open output file
    with open(args.output_topol, 'w') as f:
        # Write comments if provided
        if args.title_comments:
            f.write(f";;; {args.title_comments}\n\n")
        
        # Write #include statements for force field files
        f.write(f'#include "{clean_path(args.path_ff, args.ff)}"\n')
        f.write(f'#include "{clean_path(args.path_ff, args.ions)}"\n')
        f.write(f'#include "{clean_path(args.path_ff, args.solvent)}"\n')
        
        # Write #include statement for ligand if not None
        if args.itp_ligand.lower() != "none":
            f.write(f'#include "{args.itp_ligand}"\n')
        
        f.write("\n[ system ]\n")
        f.write("; Name\n")
        if args.title_system:
            f.write(f"{args.title_system}\n")
        else:
            f.write("Unnamed system\n")
        
        f.write("\n[ molecules ]\n")
        f.write("; Compound        #mols\n")
        f.write(f"{args.name_molecule}               {args.number_molecule}\n")

    print(f"Topology file '{args.output_topol}' successfully generated.")

if __name__ == "__main__":
    main()
