#!/usr/bin/env python3

"""
Calculate bonds, angles, and dihedrals from molecular dynamics trajectories.

Usage:
    python3 generate_bonds_angles_dihedrals.py \\
        --bonds_ndx      bonds.ndx          \\
        --angles_ndx     angles.ndx         \\
        --dihedrals_ndx  dihedrals.ndx      \\
        --xtc_file       md.xtc             \\
        --tpr_file       md.tpr             \\
                                           \\
        # Optional (required if --remove_pbc is activated): \\
        # Removes PBC and fits the molecule to eliminate rotation/translation \\
        # Recommended if your system is not aligned \\
                                           \\
        --index          None (default) or index.ndx   \\
        --remove_pbc                         # Flag to remove PBC \\
        --group_1        "name group"                   \\
        --group_2        "System"                       \\
        --keep_intermediate                   # Keeps all intermediate files
"""

import argparse
import subprocess
import os
import sys
import re
from pathlib import Path

def setup_directories(analysis_type):
    """Setup directories for each analysis type without removing existing ones"""
    dir_name = analysis_type  # bonds, angles, dihedrals
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    return dir_name

def extract_indices(input_file, pattern):
    """Extract indices from file based on pattern"""
    with open(input_file, 'r') as f:
        content = f.read()
    
    # Find all lines matching the pattern
    matches = re.findall(pattern, content, re.MULTILINE)
    indices_list = []
    for match in matches:
        # Convert to list of integers
        indices = [int(x) for x in match.strip().split()]
        indices_list.append(indices)
    
    return indices_list

def remove_pbc_and_align(xtc_file, tpr_file, index_file, group_1, group_2, keep_intermediate):
    """
    Remove PBC and align the trajectory by fitting to remove rotation and translation.
    Returns the path to the processed trajectory file.
    """
    print("\n=== Removing PBC and Aligning Trajectory ===")
    
    # Define output files
    base_name = os.path.splitext(xtc_file)[0]
    whole_xtc = f"{base_name}_whole.xtc"
    center_fit_xtc = f"{base_name}_center_fit.xtc"
    reference_pdb = f"{base_name}_t-0ns.pdb"
    
    # Step 1: Remove PBC to make molecules whole
    print("Step 1: Removing PBC (making molecules whole)...")
    if index_file and os.path.exists(index_file):
        cmd1 = f"echo {group_2} | gmx trjconv -s {tpr_file} -f {xtc_file} -o {whole_xtc} -pbc whole -n {index_file}"
    else:
        cmd1 = f"echo {group_2} | gmx trjconv -s {tpr_file} -f {xtc_file} -o {whole_xtc} -pbc whole"
    
    result1 = subprocess.run(cmd1, shell=True, capture_output=True, text=True)
    if result1.returncode != 0:
        print("ERROR: Failed to remove PBC")
        print(result1.stderr)
        sys.exit(1)
    print("PBC removal completed successfully")
    
    # Step 2: Align trajectory (fit rotation + translation)
    print("Step 2: Aligning trajectory (removing rotation and translation)...")
    if index_file and os.path.exists(index_file):
        cmd2 = f"echo {group_1} {group_2} | gmx trjconv -s {tpr_file} -f {whole_xtc} -o {center_fit_xtc} -fit rot+trans -n {index_file}"
    else:
        cmd2 = f"echo {group_1} {group_2} | gmx trjconv -s {tpr_file} -f {whole_xtc} -o {center_fit_xtc} -fit rot+trans"
    
    result2 = subprocess.run(cmd2, shell=True, capture_output=True, text=True)
    if result2.returncode != 0:
        print("ERROR: Failed to align trajectory")
        print(result2.stderr)
        sys.exit(1)
    print("Alignment completed successfully")
    
    # Step 3: Extract reference structure at t=0 ns (optional, for visualization)
    print("Step 3: Extracting reference structure at t=0 ns...")
    if index_file and os.path.exists(index_file):
        cmd3 = f"echo {group_2} | gmx trjconv -s {tpr_file} -f {center_fit_xtc} -o {reference_pdb} -dump 0 -n {index_file}"
    else:
        cmd3 = f"echo {group_2} | gmx trjconv -s {tpr_file} -f {center_fit_xtc} -o {reference_pdb} -dump 0"
    
    subprocess.run(cmd3, shell=True, capture_output=True, text=True)
    print(f"Reference structure saved as {reference_pdb}")
    
    # Clean up intermediate files if not keeping them
    if not keep_intermediate:
        print("Cleaning up intermediate files...")
        if os.path.exists(whole_xtc) and whole_xtc != center_fit_xtc:
            os.remove(whole_xtc)
    else:
        print(f"Intermediate files kept: {whole_xtc}, {center_fit_xtc}, {reference_pdb}")
    
    return center_fit_xtc

def process_bonds(bonds_ndx, xtc_file, tpr_file):
    """Process bond distances"""
    print("\n=== Processing Bonds ===")
    
    dir_name = setup_directories("bonds")
    report_path = os.path.join(dir_name, "report_bonds.txt")
    data_path = os.path.join(dir_name, "data_bonds.txt")
    errors_path = os.path.join(dir_name, "errors.log")
    
    # Extract bonds (lines with exactly 2 numbers)
    bonds = extract_indices(bonds_ndx, r'^\s*(\d+\s+\d+)\s*$')
    
    print(f"Found {len(bonds)} bonds")
    
    # Create empty files if they don't exist (append mode later)
    if not os.path.exists(report_path):
        Path(report_path).touch()
    if not os.path.exists(data_path):
        Path(data_path).touch()
    
    for i, bond in enumerate(bonds):
        bead_pair = f"{bond[0]} {bond[1]}"
        print(f"Processing bond {i}: {bead_pair}")
        
        # Create temporary index file
        temp_ndx = f"temp_bond_{i}.ndx"
        with open(temp_ndx, 'w') as f:
            f.write(f"[ bond_{i} ]\n")
            f.write(f"{bead_pair}\n")
        
        # Run gmx distance
        cmd = f"echo 0 | gmx distance -f {xtc_file} -n {temp_ndx} -s {tpr_file} -oall {dir_name}/bond_{i}.xvg -xvg none"
        log_path = os.path.join(dir_name, f"bond_{i}.log")
        
        with open(log_path, 'w') as log_file:
            result = subprocess.run(cmd, shell=True, stdout=log_file, stderr=subprocess.STDOUT)
        
        bond_xvg = os.path.join(dir_name, f"bond_{i}.xvg")
        
        if result.returncode != 0 or not os.path.exists(bond_xvg):
            with open(errors_path, 'a') as errors_file:
                errors_file.write(f"ERROR on bond {i}\n")
                with open(log_path, 'r') as log_file:
                    errors_file.write(log_file.read())
        else:
            # Extract average and std deviation from log
            with open(log_path, 'r') as log_file:
                log_content = log_file.read()
                avg_match = re.search(r'Average distance\s+([\d.]+)', log_content)
                std_match = re.search(r'Standard deviation\s+([\d.]+)', log_content)
                
                if avg_match and std_match:
                    with open(data_path, 'a') as data_file:
                        data_file.write(f"---- bond {i} ({bead_pair}) ----\n")
                        data_file.write(f"{avg_match.group(1)}\n")
                        data_file.write(f"{std_match.group(1)}\n")
            
            # Run gmx analyze for distribution
            dist_cmd = f"gmx analyze -f {dir_name}/bond_{i}.xvg -dist {dir_name}/distr_bond_{i}.xvg -xvg none -bw 0.001"
            subprocess.run(dist_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            
            # Write to report
            with open(report_path, 'a') as report_file:
                report_file.write(f"{i}: {bead_pair}\n")
        
        # Clean up temp file
        os.remove(temp_ndx)
    
    print(f"Bonds processing completed. Results in {dir_name}/")

def process_angles(angles_ndx, xtc_file):
    """Process angles"""
    print("\n=== Processing Angles ===")
    
    dir_name = setup_directories("angles")
    report_path = os.path.join(dir_name, "report_angles.txt")
    data_path = os.path.join(dir_name, "data_angles.txt")
    errors_path = os.path.join(dir_name, "errors.log")
    
    # Extract angles (lines with exactly 3 numbers)
    angles = extract_indices(angles_ndx, r'^\s*(\d+\s+\d+\s+\d+)\s*$')
    
    print(f"Found {len(angles)} angles")
    
    if not os.path.exists(report_path):
        Path(report_path).touch()
    if not os.path.exists(data_path):
        Path(data_path).touch()
    
    for i, angle in enumerate(angles):
        bead_trio = f"{angle[0]} {angle[1]} {angle[2]}"
        print(f"Processing angle {i}: {bead_trio}")
        
        # Create temporary index file
        temp_ndx = f"temp_angle_{i}.ndx"
        with open(temp_ndx, 'w') as f:
            f.write(f"[ angle_{i} ]\n")
            f.write(f"{bead_trio}\n")
        
        # Run gmx angle
        cmd = f"echo 0 | gmx angle -f {xtc_file} -n {temp_ndx} -ov {dir_name}/ang_{i}.xvg"
        log_path = os.path.join(dir_name, f"ang_{i}.log")
        
        with open(log_path, 'w') as log_file:
            result = subprocess.run(cmd, shell=True, stdout=log_file, stderr=subprocess.STDOUT)
        
        ang_xvg = os.path.join(dir_name, f"ang_{i}.xvg")
        
        if result.returncode != 0 or not os.path.exists(ang_xvg):
            with open(errors_path, 'a') as errors_file:
                errors_file.write(f"ERROR on angle {i}\n")
                with open(log_path, 'r') as log_file:
                    errors_file.write(log_file.read())
        else:
            # Extract average and std deviation from log
            with open(log_path, 'r') as log_file:
                log_content = log_file.read()
                avg_match = re.search(r'< angle >\s+([\d.]+)', log_content)
                std_match = re.search(r'Std\. Dev\.\s+([\d.]+)', log_content)
                
                if avg_match and std_match:
                    with open(data_path, 'a') as data_file:
                        data_file.write(f"---- ang {i} ({bead_trio}) ----\n")
                        data_file.write(f"{avg_match.group(1)}\n")
                        data_file.write(f"{std_match.group(1)}\n")
            
            # Run gmx analyze for distribution
            dist_cmd = f"gmx analyze -f {dir_name}/ang_{i}.xvg -dist {dir_name}/distr_ang_{i}.xvg -xvg none -bw 1.0"
            subprocess.run(dist_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            
            # Write to report
            with open(report_path, 'a') as report_file:
                report_file.write(f"{i}: {bead_trio}\n")
        
        # Clean up temp file
        os.remove(temp_ndx)
    
    print(f"Angles processing completed. Results in {dir_name}/")

def process_dihedrals(dihedrals_ndx, xtc_file):
    """Process dihedrals"""
    print("\n=== Processing Dihedrals ===")
    
    dir_name = setup_directories("dihedrals")
    report_path = os.path.join(dir_name, "report_dihedrals.txt")
    data_path = os.path.join(dir_name, "data_dihedrals.txt")
    errors_path = os.path.join(dir_name, "errors.log")
    
    # Extract dihedrals (lines with exactly 4 numbers)
    dihedrals = extract_indices(dihedrals_ndx, r'^\s*(\d+\s+\d+\s+\d+\s+\d+)\s*$')
    
    print(f"Found {len(dihedrals)} dihedrals")
    
    if not os.path.exists(report_path):
        Path(report_path).touch()
    if not os.path.exists(data_path):
        Path(data_path).touch()
    
    for i, dihedral in enumerate(dihedrals):
        bead_quartet = f"{dihedral[0]} {dihedral[1]} {dihedral[2]} {dihedral[3]}"
        print(f"Processing dihedral {i}: {bead_quartet}")
        
        # Validate bead numbers (1-30 as in original)
        valid = True
        for num in dihedral:
            if num < 1 or num > 30:
                valid = False
                with open(errors_path, 'a') as errors_file:
                    errors_file.write(f"ERROR: bead {num} out of range (1-30)\n")
        
        if valid:
            # Create temporary index file
            temp_ndx = f"temp_dih_{i}.ndx"
            with open(temp_ndx, 'w') as f:
                f.write(f"[ dihedral_{i} ]\n")
                f.write(f"{bead_quartet}\n")
            
            # Run gmx angle for dihedrals
            cmd = f"echo 0 | gmx angle -type dihedral -f {xtc_file} -n {temp_ndx} -ov {dir_name}/dih_{i}.xvg"
            log_path = os.path.join(dir_name, f"dih_{i}.log")
            
            with open(log_path, 'w') as log_file:
                result = subprocess.run(cmd, shell=True, stdout=log_file, stderr=subprocess.STDOUT)
            
            dih_xvg = os.path.join(dir_name, f"dih_{i}.xvg")
            
            if result.returncode != 0 or not os.path.exists(dih_xvg):
                with open(errors_path, 'a') as errors_file:
                    errors_file.write(f"ERROR on dihedral {i}\n")
                    with open(log_path, 'r') as log_file:
                        errors_file.write(log_file.read())
            else:
                # Extract average and std deviation from log
                with open(log_path, 'r') as log_file:
                    log_content = log_file.read()
                    avg_match = re.search(r'< angle >\s+([\d.]+)', log_content)
                    std_match = re.search(r'Std\. Dev\.\s+([\d.]+)', log_content)
                    
                    if avg_match and std_match:
                        with open(data_path, 'a') as data_file:
                            data_file.write(f"---- dih {i} ({bead_quartet}) ----\n")
                            data_file.write(f"{avg_match.group(1)}\n")
                            data_file.write(f"{std_match.group(1)}\n")
                
                # Run gmx analyze for distribution
                dist_cmd = f"gmx analyze -f {dir_name}/dih_{i}.xvg -dist {dir_name}/distr_dih_{i}.xvg -xvg none -bw 1.0"
                subprocess.run(dist_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                
                # Write to report
                with open(report_path, 'a') as report_file:
                    report_file.write(f"{i}: {bead_quartet}\n")
            
            # Clean up temp file
            os.remove(temp_ndx)
    
    print(f"Dihedrals processing completed. Results in {dir_name}/")

def main():
    parser = argparse.ArgumentParser(
        description='Calculate bonds, angles, and dihedrals from molecular dynamics trajectories.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage
    python3 generate_bonds_angles_dihedrals.py \\
        --bonds_ndx bonds.ndx \\
        --angles_ndx angles.ndx \\
        --dihedrals_ndx dihedrals.ndx \\
        --xtc_file md.xtc \\
        --tpr_file md.tpr
    
    # With PBC removal and alignment
    python3 generate_bonds_angles_dihedrals.py \\
        --bonds_ndx bonds.ndx \\
        --angles_ndx angles.ndx \\
        --dihedrals_ndx dihedrals.ndx \\
        --xtc_file md.xtc \\
        --tpr_file md.tpr \\
        --index index.ndx \\
        --remove_pbc \\
        --group_1 "Backbone" \\
        --group_2 "System" \\
        --keep_intermediate
        """
    )
    
    # Required arguments
    parser.add_argument('--bonds_ndx', required=True, help='Bonds index file')
    parser.add_argument('--angles_ndx', required=True, help='Angles index file')
    parser.add_argument('--dihedrals_ndx', required=True, help='Dihedrals index file')
    parser.add_argument('--xtc_file', required=True, help='XTC trajectory file')
    parser.add_argument('--tpr_file', required=True, help='TPR topology file')
    
    # Optional arguments for PBC removal and alignment
    parser.add_argument('--index', default=None, help='Index file for selecting groups (required if --remove_pbc is used)')
    parser.add_argument('--remove_pbc', action='store_true', help='Remove PBC and align trajectory (recommended for non-aligned systems)')
    parser.add_argument('--group_1', default=None, help='Group for fitting (e.g., "Backbone", "C-alpha", or a custom name)')
    parser.add_argument('--group_2', default="System", help='Group for output (default: "System")')
    parser.add_argument('--keep_intermediate', action='store_true', help='Keep all intermediate files (whole.xtc, center_fit.xtc, reference.pdb)')
    
    args = parser.parse_args()
    
    # Check if required files exist
    required_files = [args.bonds_ndx, args.angles_ndx, args.dihedrals_ndx, 
                      args.xtc_file, args.tpr_file]
    for file_path in required_files:
        if not os.path.exists(file_path):
            print(f"Error: File not found - {file_path}")
            sys.exit(1)
    
    # Process trajectory if --remove_pbc is activated
    xtc_file_to_use = args.xtc_file
    
    if args.remove_pbc:
        # Check required arguments for PBC removal
        if args.group_1 is None:
            print("ERROR: --group_1 is required when --remove_pbc is activated")
            print("Please specify a group for fitting (e.g., --group_1 'Backbone')")
            sys.exit(1)
        
        if args.index is None:
            print("WARNING: No index file provided. Continuing without index selection.")
            print("This may cause issues if your system has multiple molecules.")
        
        # Remove PBC and align trajectory
        xtc_file_to_use = remove_pbc_and_align(
            args.xtc_file, 
            args.tpr_file, 
            args.index, 
            args.group_1, 
            args.group_2, 
            args.keep_intermediate
        )
    else:
        print("\n=== Using original trajectory without PBC removal ===")
        print("Note: If your system has PBC artifacts or is not aligned, consider using --remove_pbc")
    
    # Process all analyses with the (potentially processed) trajectory
    process_bonds(args.bonds_ndx, xtc_file_to_use, args.tpr_file)
    process_angles(args.angles_ndx, xtc_file_to_use)
    process_dihedrals(args.dihedrals_ndx, xtc_file_to_use)
    
    print("\n" + "="*50)
    print("Processing completed successfully!")
    print("Results are located in the following directories:")
    print("  - bonds/")
    print("  - angles/")
    print("  - dihedrals/")
    print("\nCheck errors.log in each directory if any warnings or errors occurred.")
    print("="*50)

if __name__ == "__main__":
    main()
