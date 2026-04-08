#!/usr/bin/env python3

import re
import argparse

# -------------------------
# command line arguments
# -------------------------

parser = argparse.ArgumentParser(
    description="Filter and renumber atoms in a GROMACS .itp topology using atom index intervals."
)

parser.add_argument(
    "--input",
    required=True,
    help="Input .itp file"
)

parser.add_argument(
    "--output",
    required=True,
    help="Output filtered .itp file"
)

parser.add_argument(
    "--interval",
    required=True,
    help='Atom index intervals, e.g. "[465:477]U[675:866]"'
)

args = parser.parse_args()

input_file = args.input
output_file = args.output
interval_string = args.interval


# -------------------------
# parse intervals
# -------------------------

intervals = []
for start, end in re.findall(r"\[(\d+):(\d+)\]", interval_string):
    intervals.append((int(start), int(end)))


def in_intervals(n):
    for a, b in intervals:
        if a >= 0 and a <= n <= b:
            return True
    return False


# -------------------------
# first pass: build atom map
# -------------------------

atom_map = {}
new_index = 1

inside_atoms = False

with open(input_file) as f:

    for line in f:

        stripped = line.strip()

        if stripped.startswith("[ atoms ]"):
            inside_atoms = True
            continue

        if inside_atoms and stripped.startswith("["):
            inside_atoms = False

        if inside_atoms:

            if stripped.startswith(";") or stripped == "":
                continue

            data = line.split(";")[0].strip()
            parts = data.split()

            if not parts:
                continue

            try:
                atom_id = int(parts[0])
            except:
                continue

            if in_intervals(atom_id):
                atom_map[atom_id] = new_index
                new_index += 1


# -------------------------
# second pass: write output
# -------------------------

section = None

with open(input_file) as fin, open(output_file, "w") as fout:

    for line in fin:

        stripped = line.strip()

        # detect sections

        if stripped.startswith("["):

            if "[ atoms ]" in stripped:
                section = "atoms"

            elif "[ bonds ]" in stripped:
                section = "bonds"

            elif "[ pairs ]" in stripped:
                section = "pairs"

            elif "[ angles ]" in stripped:
                section = "angles"

            elif "[ dihedrals ]" in stripped:
                section = "dihedrals"

            elif "[ cmap ]" in stripped:
                section = "cmap"

            elif "[ position_restraints ]" in stripped:
                section = "posres"

            elif "[ dihedral_restraints ]" in stripped:
                section = "dihres"

            else:
                section = None

            fout.write(line)
            continue


        # preserve comments / preprocessor

        if stripped.startswith(";") or stripped == "" or stripped.startswith("#"):
            fout.write(line)
            continue


        data = line.split(";")[0].strip()
        parts = data.split()

        if not parts:
            continue

        try:

            # -------------------------
            # atoms
            # -------------------------

            if section == "atoms":

                old = int(parts[0])

                if old in atom_map:

                    parts[0] = str(atom_map[old])
                    parts[5] = str(atom_map[old])  # update cgnr as well

                    fout.write("{:<6s} {:<8s} {:<6s} {:<8s} {:<6s} {:<6s} {:<10s} {:<10s}\n".format(*parts[:8]))


            # -------------------------
            # bonds / pairs
            # -------------------------

            elif section in ["bonds", "pairs"]:

                a1 = int(parts[0])
                a2 = int(parts[1])

                if a1 in atom_map and a2 in atom_map:

                    parts[0] = str(atom_map[a1])
                    parts[1] = str(atom_map[a2])

                    fout.write(" ".join(parts) + "\n")


            # -------------------------
            # angles
            # -------------------------

            elif section == "angles":

                a1, a2, a3 = map(int, parts[:3])

                if a1 in atom_map and a2 in atom_map and a3 in atom_map:

                    parts[0] = str(atom_map[a1])
                    parts[1] = str(atom_map[a2])
                    parts[2] = str(atom_map[a3])

                    fout.write(" ".join(parts) + "\n")


            # -------------------------
            # dihedrals
            # -------------------------

            elif section == "dihedrals":

                a1, a2, a3, a4 = map(int, parts[:4])

                if (a1 in atom_map and
                    a2 in atom_map and
                    a3 in atom_map and
                    a4 in atom_map):

                    parts[0] = str(atom_map[a1])
                    parts[1] = str(atom_map[a2])
                    parts[2] = str(atom_map[a3])
                    parts[3] = str(atom_map[a4])

                    fout.write(" ".join(parts) + "\n")


            # -------------------------
            # cmap
            # -------------------------

            elif section == "cmap":

                a1, a2, a3, a4, a5 = map(int, parts[:5])

                if (a1 in atom_map and
                    a2 in atom_map and
                    a3 in atom_map and
                    a4 in atom_map and
                    a5 in atom_map):

                    parts[0] = str(atom_map[a1])
                    parts[1] = str(atom_map[a2])
                    parts[2] = str(atom_map[a3])
                    parts[3] = str(atom_map[a4])
                    parts[4] = str(atom_map[a5])

                    fout.write(" ".join(parts) + "\n")


            # -------------------------
            # position restraints
            # -------------------------

            elif section == "posres":

                a1 = int(parts[0])

                if a1 in atom_map:

                    parts[0] = str(atom_map[a1])
                    fout.write(" ".join(parts) + "\n")


            # -------------------------
            # dihedral restraints
            # -------------------------

            elif section == "dihres":

                a1, a2, a3, a4 = map(int, parts[:4])

                if (a1 in atom_map and
                    a2 in atom_map and
                    a3 in atom_map and
                    a4 in atom_map):

                    parts[0] = str(atom_map[a1])
                    parts[1] = str(atom_map[a2])
                    parts[2] = str(atom_map[a3])
                    parts[3] = str(atom_map[a4])

                    fout.write(" ".join(parts) + "\n")


            else:
                fout.write(line)

        except:
            continue


print(f"Filtered topology written to: {output_file}")
print(f"Total atoms kept: {len(atom_map)}")
