#!/usr/bin/env python3
"""
Map atomistic trajectory to coarse-grained (CG) resolution using GROMACS.
Automatically generate CG .gro and .xtc mapped from AA MD simulation.

It is necessary to have:

- md.xtc (from AA model, with removing the pbc)
- cg.ndx (bead definitions obtained from cg-martini3 or cgbuilder)
- md.tpr (from AA model simulation)

(1) Considering to correct md.xtc (with pbc) to obtain only trajectory mapped from cg.ndx (with flag --pbc) 
python ../1-map_aa_traj_to_cg.py --index_cg cg.ndx --aa_tpr md.tpr --aa_xtc md.xtc --output_mapped mapped.xtc --remove_pbc --output_cg_gro molecule.gro

(2) Considering md_no_pbc.xtc (with pbc removed) to obtain only trajectory mapped from cg.ndx (with flag --no-pbc)
python ../1-map_aa_traj_to_cg.py --index_cg cg.ndx --aa_tpr md.tpr --aa_xtc md.xtc --output_mapped mapped.xtc --corrected_pbc --output_cg_gro molecule.gro
"""

import argparse
import subprocess
import os
import sys
import tempfile
import shutil
import re


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Map atomistic trajectory to coarse-grained (CG) resolution",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic mapping only
  %(prog)s -i cg.ndx -t md.tpr -x md.xtc
  
  # Generate CG .gro automatically after mapping
  %(prog)s -i cg.ndx -t md.tpr -x md.xtc --output_cg_gro molecule.gro
  
  # Skip PBC correction (use when trajectory already whole)
  %(prog)s -i cg.ndx -t md.tpr -x md.xtc --corrected_pbc --output_cg_gro molecule.gro
        """
    )
    
    # Required arguments
    parser.add_argument("--index_cg", "-i", required=True, help="CG index file (.ndx)")
    parser.add_argument("--aa_tpr", "-t", required=True, help="Atomistic .tpr file")
    parser.add_argument("--aa_xtc", "-x", required=True, help="Atomistic trajectory (.xtc)")
    
    # Optional arguments
    parser.add_argument("--output_mapped", "-o", default="mapped.xtc", help="Output CG trajectory (default: mapped.xtc)")
    parser.add_argument("--output_cg_gro", default="molecule.gro", help="Output CG .gro file (first frame) - auto-generated after mapping")
    
    # PBC options
    parser.add_argument("--remove_pbc", action="store_true", default=True, 
                        help="Remove PBC with gmx trjconv -pbc whole (default: True)")
    parser.add_argument("--corrected_pbc", action="store_false", dest="remove_pbc",
                        help="Use when trajectory already has corrected PBC (skip trjconv)")
    
    # Other options
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--dry-run", action="store_true", help="Dry run (print commands only)")
    parser.add_argument("--keep-temp", action="store_true", help="Keep temporary files")
    parser.add_argument("--gmx_cmd", default="gmx", help="GROMACS command (default: gmx)")
    
    return parser.parse_args()


def validate_args(args):
    """Validate arguments and auto-generate missing ones."""
    # Check required files exist
    for f, name in [(args.index_cg, "index"), (args.aa_tpr, "AA .tpr"), (args.aa_xtc, "AA .xtc")]:
        if not os.path.exists(f):
            print(f"ERROR: {name} file not found: {f}")
            sys.exit(1)
    
    return args


def count_beads(index_file):
    """Count number of bead groups in index file."""
    with open(index_file, 'r') as f:
        return sum(1 for line in f if line.strip().startswith('[') and line.strip().endswith(']'))


def run_command_with_sequential_input(cmd, num_inputs, verbose=False, dry_run=False, description=None):
    """
    Run a command that requires sequential numeric inputs (0 to num_inputs-1).
    This mimics: seq 0 $((num_inputs-1)) | gmx command
    """
    if description:
        print(f"\n>>> {description}")
    
    # Create the input string with newlines (simulating seq)
    input_data = '\n'.join(str(i) for i in range(num_inputs))
    
    # Use printf instead of echo for better newline handling
    cmd_str = f"printf '{input_data}\\n' | {cmd}"
    
    if verbose or dry_run:
        print(f"Command: {cmd_str}")
    
    if dry_run:
        return True
    
    try:
        result = subprocess.run(cmd_str, shell=True, check=True, capture_output=True, text=True)
        if verbose and result.stdout:
            print(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        print(f"ERROR: {e.stderr}")
        return False


def map_trajectory_to_cg(args, no_of_beads, traj_file):
    """Map trajectory to CG using gmx traj."""
    cmd = f"{args.gmx_cmd} traj -f {traj_file} -s {args.aa_tpr} -oxt {args.output_mapped} -n {args.index_cg} -ng {no_of_beads} -com"
    
    return run_command_with_sequential_input(cmd, no_of_beads, args.verbose, args.dry_run, 
                                               f"Mapping {no_of_beads} beads to CG")


def extract_cg_gro_after_mapping(args, no_of_beads):
    """
    Extract first frame as CG .gro using ORIGINAL AA trajectory.
    This avoids mismatch between mapped.xtc and md.tpr.
    """

    print(f"\n>>> Extracting first frame as CG .gro (using AA trajectory)")

    cmd = (
        f"{args.gmx_cmd} traj "
        f"-f {args.aa_xtc} "
        f"-s {args.aa_tpr} "
        f"-oxt {args.output_cg_gro} "
        f"-n {args.index_cg} "
        f"-ng {no_of_beads} "
        f"-com -e 0"
    )

    if run_command_with_sequential_input(
        cmd,
        no_of_beads,
        args.verbose,
        args.dry_run,
        "Extracting first frame (AA → CG COM)"
    ):
        if os.path.exists(args.output_cg_gro):
            with open(args.output_cg_gro, 'r') as f:
                lines = f.readlines()
                if len(lines) >= 2:
                    print(f"✓ CG .gro saved to: {args.output_cg_gro} ({lines[1].strip()} atoms)")
                else:
                    print(f"✓ CG .gro saved to: {args.output_cg_gro}")
            return True
        else:
            print(f"✗ Failed to generate CG .gro")
            return False
    else:
        print(f"✗ Failed to generate CG .gro")
        return False

def remove_pbc(args):
    """Remove PBC from trajectory using gmx trjconv."""
    temp_dir = tempfile.mkdtemp(prefix="gmx_map_")
    whole_xtc = os.path.join(temp_dir, "whole.xtc")
    
    cmd = f"{args.gmx_cmd} trjconv -f {args.aa_xtc} -s {args.aa_tpr} -o {whole_xtc} -pbc whole"
    
    # For trjconv, we need to select "0" (system)
    input_data = "0"
    cmd_str = f"echo '{input_data}' | {cmd}"
    
    if args.verbose or args.dry_run:
        print(f"\n>>> Removing PBC (gmx trjconv -pbc whole)")
        print(f"Command: {cmd_str}")
    
    if args.dry_run:
        return temp_dir, whole_xtc
    
    try:
        result = subprocess.run(cmd_str, shell=True, check=True, capture_output=True, text=True)
        if args.verbose and result.stdout:
            print(result.stdout)
        return temp_dir, whole_xtc
    except subprocess.CalledProcessError as e:
        print(f"WARNING: PBC correction failed: {e.stderr}")
        print("Using original trajectory")
        shutil.rmtree(temp_dir)
        return None, args.aa_xtc


def main():
    args = parse_arguments()
    args = validate_args(args)
    
    # Count beads
    no_of_beads = count_beads(args.index_cg)
    print(f"Found {no_of_beads} bead groups in {args.index_cg}")
    
    temp_dir = None
    traj_for_mapping = args.aa_xtc
    
    # Step 1: Apply PBC correction if requested
    if args.remove_pbc:
        temp_dir, traj_for_mapping = remove_pbc(args)
    else:
        print("Skipping PBC correction (using trajectory as is)")
    
    # Step 2: Map trajectory to CG
    if os.path.exists(args.output_mapped) and not args.dry_run:
        print(f"Using existing CG trajectory: {args.output_mapped}")
    else:
        if map_trajectory_to_cg(args, no_of_beads, traj_for_mapping):
            print(f"✓ CG trajectory saved to: {args.output_mapped}")
        else:
            print("ERROR: Trajectory mapping failed")
            if temp_dir and not args.keep_temp:
                shutil.rmtree(temp_dir)
            sys.exit(1)
    
    # Step 3: Automatically extract CG .gro (first frame) after mapping
    if args.output_cg_gro:
        if extract_cg_gro_after_mapping(args, no_of_beads):
            print(f"\n✓ CG .gro file generated automatically after mapping: {args.output_cg_gro}")
        else:
            print(f"✗ Failed to generate CG .gro automatically")
            if temp_dir and not args.keep_temp:
                shutil.rmtree(temp_dir)
            sys.exit(1)
    
    # Cleanup
    if temp_dir and not args.keep_temp and not args.dry_run:
        shutil.rmtree(temp_dir)
        if args.verbose:
            print(f"\nRemoved temporary directory: {temp_dir}")
    elif temp_dir and args.keep_temp:
        print(f"\nTemporary files kept in: {temp_dir}")
    
    # Summary
    print("\n" + "="*50)
    print("SUMMARY")
    print("="*50)
    print(f"Input files:")
    print(f"  - Index file:      {args.index_cg}")
    print(f"  - AA .tpr:         {args.aa_tpr}")
    print(f"  - AA .xtc:         {args.aa_xtc}")
    print(f"\nParameters:")
    print(f"  - Number of beads: {no_of_beads}")
    print(f"  - PBC correction:  {'removed' if args.remove_pbc else 'skipped'}")
    print(f"\nOutput files:")
    
    for name, path in [("CG trajectory", args.output_mapped),
                       ("CG .gro", args.output_cg_gro)]:
        if path and os.path.exists(path):
            size = os.path.getsize(path) / (1024*1024)
            print(f"  - {name}: {path} ({size:.2f} MB)")
        elif path:
            print(f"  - {name}: {path} [NOT CREATED]")
    
    if args.dry_run:
        print("\n*** DRY RUN - No files were actually created ***")


if __name__ == "__main__":
    main()
