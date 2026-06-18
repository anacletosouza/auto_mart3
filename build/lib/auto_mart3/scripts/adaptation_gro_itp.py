#!/usr/bin/env python3

"""
Update atom names in a GROMACS ITP file to match those from a reference GRO file.

This script ensures consistency between topology (.itp) and coordinate (.gro)
files by replacing atom names in the [ atoms ] section of the ITP file with
the corresponding names from the GRO file. This is useful when atom naming
differs between parametrization and structural files.

The mapping is performed sequentially, assuming that atoms are ordered
identically in both files.

Usage:
    python 5-adaptation_gro_itp.py \
        --input_itp input.itp \
        --input_gro_ref structure.gro \
        --output_itp_adapted output.itp

Arguments:
    --input_itp             Input ITP file to be modified
    --input_gro_ref        Reference GRO file containing correct atom names
    --output_itp_adapted   Output ITP file with updated atom names
    --verbose              Print detailed debug information

Output:
    A new ITP file where atom names in the [ atoms ] section are updated
    to match the GRO reference.

Notes:
    - The number and order of atoms must be identical in both ITP and GRO files.
    - Only the atom name field (5th column in [ atoms ]) is modified.
    - The script preserves formatting and ignores comments or malformed lines.
    - An error is raised if the number of atoms does not match.
"""

import argparse
import sys


def read_gro_atom_names(gro_file):
    names = []
    with open(gro_file) as f:
        lines = f.readlines()

    if len(lines) < 3:
        raise ValueError("Invalid GRO file")

    atom_lines = lines[2:-1]

    for line in atom_lines:
        if len(line) < 15:
            continue
        names.append(line[10:15].strip())

    return names


def fix_itp_atom_names(itp_file, gro_names, output_file, verbose=False):
    with open(itp_file) as f:
        lines = f.readlines()

    new_lines = []
    atom_index = 0
    in_atoms = False
    total_atoms_itp = 0

    for line in lines:
        stripped = line.strip()

        # Detect start of [ atoms ]
        if stripped.startswith("[ atoms ]"):
            in_atoms = True
            new_lines.append(line)
            continue

        # Detect end of section
        if in_atoms and stripped.startswith("[") and not stripped.startswith("[ atoms ]"):
            in_atoms = False

        if in_atoms:
            if stripped.startswith(";") or stripped == "":
                new_lines.append(line)
                continue

            parts = line.split()

            if len(parts) >= 8:
                if atom_index >= len(gro_names):
                    raise ValueError(
                        f"More atoms in ITP than in GRO (index {atom_index})"
                    )

                old_name = parts[4]
                new_name = gro_names[atom_index]

                parts[4] = new_name
                atom_index += 1
                total_atoms_itp += 1

                if verbose and old_name != new_name:
                    print(f"[DEBUG] Atom {parts[0]}: {old_name} → {new_name}")

                new_line = "{:>5} {:<6} {:>5} {:<6} {:<6} {:>5} {:>10} {:>10}\n".format(
                    parts[0], parts[1], parts[2], parts[3],
                    parts[4], parts[5], parts[6], parts[7]
                )

                new_lines.append(new_line)
            else:
                new_lines.append(line)
        else:
            new_lines.append(line)

    if atom_index != len(gro_names):
        raise ValueError(
            f"Mismatch: GRO has {len(gro_names)} atoms but ITP has {atom_index}"
        )

    with open(output_file, "w") as f:
        f.writelines(new_lines)

    return total_atoms_itp


def main():
    parser = argparse.ArgumentParser(
        description="Adapt ITP atom names to match GRO reference"
    )

    parser.add_argument("--input_itp", required=True, help="Input ITP file")
    parser.add_argument("--input_gro_ref", required=True, help="Reference GRO file")
    parser.add_argument("--output_itp_adapted", required=True, help="Output ITP file")
    parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args()

    try:
        gro_names = read_gro_atom_names(args.input_gro_ref)

        if args.verbose:
            print(f"[INFO] Read {len(gro_names)} atoms from GRO")

        n_atoms = fix_itp_atom_names(
            args.input_itp,
            gro_names,
            args.output_itp_adapted,
            verbose=args.verbose
        )

        print(f"✓ ITP successfully adapted")
        print(f"  Atoms processed: {n_atoms}")
        print(f"  Output: {args.output_itp_adapted}")

    except Exception as e:
        print(f"✗ Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
