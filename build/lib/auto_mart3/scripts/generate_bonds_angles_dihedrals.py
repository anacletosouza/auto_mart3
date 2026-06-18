#!/usr/bin/env python3

"""
Compute bond distances, angles, and dihedral distributions from molecular
dynamics trajectories using GROMACS tools.

This script reads index files defining bonds (2 atoms), angles (3 atoms),
and dihedrals (4 atoms), and extracts structural properties from a trajectory
(.xtc) using the corresponding topology (.tpr). Results include averages,
standard deviations, and distributions for each interaction.

Optionally, the trajectory can be preprocessed to remove periodic boundary
conditions (PBC) and perform structural alignment (rotation/translation fitting),
which is recommended for non-centered or non-aligned systems.

IMPORTANT:
    The trajectory (.xtc) and topology (.tpr) must correspond to the same system
    and contain the same number of atoms.

Usage:
    python3 3-generate_bonds_angles_dihedrals.py \\
        --bonds_ndx bonds.ndx \\
        --angles_ndx angles.ndx \\
        --dihedrals_ndx dihedrals.ndx \\
        --xtc_file md.xtc \\
        --tpr_file md.tpr \\
        --output_all_files "OUTPUT"

Optional preprocessing (PBC removal and alignment):
    --index index.ndx        Index file for group selection (recommended)
    --remove_pbc             Enable PBC removal and alignment
    --group_1 "GroupName"    Group used for fitting (required if --remove_pbc)
    --group_2 "GroupName"    Group for output (default: "System")
    --keep_intermediate      Keep intermediate trajectory files

Arguments:
    --bonds_ndx         Index file defining bonds (pairs of atoms)
    --angles_ndx        Index file defining angles (triplets of atoms)
    --dihedrals_ndx     Index file defining dihedrals (quartets of atoms)
    --xtc_file          Input trajectory file (.xtc)
    --tpr_file          Input topology file (.tpr)
    --output_all_files  Base output directory path (creates bonds/, angles/, dihedrals/ subdirectories)

Output:
    The output directory specified by --output_all_files will contain three subdirectories:
        bonds/       Bond distance data and distributions
        angles/      Angle data and distributions
        dihedrals/   Dihedral data and distributions

    Each directory contains:
        - report_*.txt   Mapping of indices
        - data_*.txt     Average and standard deviation values
        - distr_*.xvg    Distribution files
        - errors.log     Execution errors (if any)

Notes:
    - Temporary files are automatically cleaned unless --keep_intermediate is used.
    - GROMACS (gmx) must be available in the system PATH.
    - For accurate structural analysis, preprocessing with --remove_pbc is recommended
      when dealing with periodic systems.
    - The output directory structure: OUTPUT/bonds/, OUTPUT/angles/, OUTPUT/dihedrals/
"""
import argparse
import subprocess
import os
import sys
import re
import tempfile
import shutil
import glob
from pathlib import Path

def cleanup_backup_files():
    """Remove GROMACS backup files (# files) from current directory"""
    backup_files = glob.glob('#*') + glob.glob('*.#*') + glob.glob('*.1#') + glob.glob('*.2#') + glob.glob('*.3#') + glob.glob('*.4#') + glob.glob('*.5#')
    for file in backup_files:
        try:
            os.remove(file)
            print(f"Removed old backup file: {file}")
        except OSError:
            pass

def setup_directories(analysis_type, output_base_dir):
    """Setup directories for each analysis type without removing existing ones"""
    # Create the full path: output_base_dir/analysis_type
    dir_path = os.path.join(output_base_dir, analysis_type)
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    return dir_path

def clean_index_file(input_file):
    """
    Clean index file by removing headers, comments, and empty lines.
    Returns path to temporary cleaned file.
    """
    # Create a temporary file
    temp_fd, temp_path = tempfile.mkstemp(suffix='.ndx', prefix='cleaned_', text=True)
    
    with open(temp_fd, 'w') as temp_file:
        with open(input_file, 'r') as original_file:
            for line in original_file:
                # Skip lines that start with '[' (headers) or are empty
                if line.strip().startswith('[') or not line.strip():
                    continue
                # Remove comments (anything after ';')
                line_without_comments = line.split(';')[0]
                # Skip if line is empty after removing comments
                if line_without_comments.strip():
                    temp_file.write(line_without_comments + '\n')
    
    return temp_path

def extract_indices(cleaned_file, pattern):
    """Extract indices from cleaned file based on pattern"""
    with open(cleaned_file, 'r') as f:
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
    
    # Create a temporary directory for PBC processing
    temp_dir = tempfile.mkdtemp(prefix='gmx_pbc_')
    print(f"Using temporary directory: {temp_dir}")
    
    # Define output files in temporary directory
    base_name = os.path.splitext(os.path.basename(xtc_file))[0]
    whole_xtc = os.path.join(temp_dir, f"{base_name}_whole.xtc")
    center_fit_xtc = os.path.join(temp_dir, f"{base_name}_center_fit.xtc")
    reference_pdb = os.path.join(temp_dir, f"{base_name}_t-0ns.pdb")
    
    # Final output path (outside temp dir if keeping intermediate)
    if keep_intermediate:
        final_xtc = f"{base_name}_center_fit.xtc"
        final_pdb = f"{base_name}_t-0ns.pdb"
        final_whole = f"{base_name}_whole.xtc"
    else:
        final_xtc = center_fit_xtc
    
    # Step 1: Remove PBC to make molecules whole
    print("Step 1: Removing PBC (making molecules whole)...")
    if index_file and os.path.exists(index_file):
        cmd1 = f"echo {group_1} | gmx trjconv -s {tpr_file} -f {xtc_file} -o {whole_xtc} -pbc whole -n {index_file}"
    else:
        cmd1 = f"echo {group_1} | gmx trjconv -s {tpr_file} -f {xtc_file} -o {whole_xtc} -pbc whole"
    
    result1 = subprocess.run(cmd1, shell=True, capture_output=True, text=True)
    if result1.returncode != 0:
        print("ERROR: Failed to remove PBC")
        print(result1.stderr)
        shutil.rmtree(temp_dir)
        sys.exit(1)
    print("PBC removal completed successfully")
    
    # Step 2: Align trajectory (fit rotation + translation)
    print("Step 2: Aligning trajectory (removing rotation and translation)...")
    if index_file and os.path.exists(index_file):
        cmd2 = f"echo {group_2} {group_1} | gmx trjconv -s {tpr_file} -f {whole_xtc} -o {center_fit_xtc} -fit rot+trans -n {index_file}"
    else:
        cmd2 = f"echo {group_2} {group_1} | gmx trjconv -s {tpr_file} -f {whole_xtc} -o {center_fit_xtc} -fit rot+trans"
    
    result2 = subprocess.run(cmd2, shell=True, capture_output=True, text=True)
    if result2.returncode != 0:
        print("ERROR: Failed to align trajectory")
        print(result2.stderr)
        shutil.rmtree(temp_dir)
        sys.exit(1)
    print("Alignment completed successfully")
    
    # Step 3: Extract reference structure at t=0 ns (optional, for visualization)
    print("Step 3: Extracting reference structure at t=0 ns...")
    if index_file and os.path.exists(index_file):
        cmd3 = f"echo {group_2} | gmx trjconv -s {tpr_file} -f {center_fit_xtc} -o {reference_pdb} -dump 0 -n {index_file}"
    else:
        cmd3 = f"echo {group_2} | gmx trjconv -s {tpr_file} -f {center_fit_xtc} -o {reference_pdb} -dump 0"
    
    subprocess.run(cmd3, shell=True, capture_output=True, text=True)
    print(f"Reference structure extracted")
    
    # Copy files to final location if keeping intermediate
    if keep_intermediate:
        print("Copying intermediate files to current directory...")
        shutil.copy2(center_fit_xtc, final_xtc)
        shutil.copy2(reference_pdb, final_pdb)
        shutil.copy2(whole_xtc, final_whole)
        print(f"Intermediate files kept: {final_whole}, {final_xtc}, {final_pdb}")
        result_path = final_xtc
    else:
        result_path = center_fit_xtc
        print(f"Processed trajectory: {result_path}")
    
    # Clean up temporary directory
    shutil.rmtree(temp_dir)
    print("Temporary files cleaned up")
    
    return result_path

def process_bonds(bonds_ndx, xtc_file, tpr_file, output_base_dir):
    """Process bond distances"""
    print("\n=== Processing Bonds ===")
    
    # Convert to absolute paths
    bonds_ndx_abs = os.path.abspath(bonds_ndx)
    xtc_file_abs = os.path.abspath(xtc_file)
    tpr_file_abs = os.path.abspath(tpr_file)
    
    # Create a temporary directory for the entire bond processing
    temp_work_dir = tempfile.mkdtemp(prefix='bonds_work_')
    original_dir = os.getcwd()
    
    # Change to temp directory
    os.chdir(temp_work_dir)
    
    # Get the output directory path for bonds
    bonds_output_dir = setup_directories("bonds", output_base_dir)
    
    # Create bonds directory in temp
    os.makedirs("bonds", exist_ok=True)
    
    # Clean the bonds index file
    cleaned_bonds = clean_index_file(bonds_ndx_abs)
    
    # Extract bonds (lines with exactly 2 numbers)
    bonds = extract_indices(cleaned_bonds, r'^\s*(\d+\s+\d+)\s*$')
    
    print(f"Found {len(bonds)} bonds")
    
    report_path = os.path.join("bonds", "report_bonds.txt")
    data_path = os.path.join("bonds", "data_bonds.txt")
    errors_path = os.path.join("bonds", "errors.log")
    
    # Create empty files if they don't exist
    Path(report_path).touch()
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
        cmd = f"echo 0 | gmx distance -f {xtc_file_abs} -n {temp_ndx} -s {tpr_file_abs} -oall bonds/bond_{i}.xvg -xvg none"
        log_path = os.path.join("bonds", f"bond_{i}.log")
        
        with open(log_path, 'w') as log_file:
            result = subprocess.run(cmd, shell=True, stdout=log_file, stderr=subprocess.STDOUT)
        
        bond_xvg = os.path.join("bonds", f"bond_{i}.xvg")
        
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
            
            # Run gmx analyze for distribution (suppress output)
            dist_cmd = f"gmx analyze -f bonds/bond_{i}.xvg -dist bonds/distr_bond_{i}.xvg -xvg none -bw 0.001"
            subprocess.run(dist_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            
            # Write to report
            with open(report_path, 'a') as report_file:
                report_file.write(f"{i}: {bead_pair}\n")
        
        # Clean up temp file
        os.remove(temp_ndx)
    
    # Clean up cleaned index file
    os.unlink(cleaned_bonds)
    
    # Copy results back to original output directory
    if os.path.exists("bonds"):
        # Ensure the output directory exists
        os.makedirs(bonds_output_dir, exist_ok=True)
        # Copy each file individually
        for item in os.listdir("bonds"):
            src = os.path.join("bonds", item)
            dst = os.path.join(bonds_output_dir, item)
            if os.path.isfile(src):
                shutil.copy2(src, dst)
            elif os.path.isdir(src):
                dest_subdir = os.path.join(bonds_output_dir, item)
                if os.path.exists(dest_subdir):
                    shutil.rmtree(dest_subdir)
                shutil.copytree(src, dest_subdir)
    
    # Return to original directory and clean up
    os.chdir(original_dir)
    shutil.rmtree(temp_work_dir)
    
    print(f"Bonds processing completed. Results in {bonds_output_dir}")

def process_angles(angles_ndx, xtc_file, output_base_dir):
    """Process angles"""
    print("\n=== Processing Angles ===")
    
    # Convert to absolute paths
    angles_ndx_abs = os.path.abspath(angles_ndx)
    xtc_file_abs = os.path.abspath(xtc_file)
    
    # Create a temporary directory for the entire angle processing
    temp_work_dir = tempfile.mkdtemp(prefix='angles_work_')
    original_dir = os.getcwd()
    
    # Change to temp directory
    os.chdir(temp_work_dir)
    
    # Get the output directory path for angles
    angles_output_dir = setup_directories("angles", output_base_dir)
    
    # Create angles directory in temp
    os.makedirs("angles", exist_ok=True)
    
    # Clean the angles index file
    cleaned_angles = clean_index_file(angles_ndx_abs)
    
    # Extract angles (lines with exactly 3 numbers)
    angles = extract_indices(cleaned_angles, r'^\s*(\d+\s+\d+\s+\d+)\s*$')
    
    print(f"Found {len(angles)} angles")
    
    report_path = os.path.join("angles", "report_angles.txt")
    data_path = os.path.join("angles", "data_angles.txt")
    errors_path = os.path.join("angles", "errors.log")
    
    Path(report_path).touch()
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
        cmd = f"echo 0 | gmx angle -f {xtc_file_abs} -n {temp_ndx} -ov angles/ang_{i}.xvg"
        log_path = os.path.join("angles", f"ang_{i}.log")
        
        with open(log_path, 'w') as log_file:
            result = subprocess.run(cmd, shell=True, stdout=log_file, stderr=subprocess.STDOUT)
        
        ang_xvg = os.path.join("angles", f"ang_{i}.xvg")
        
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
            
            # Run gmx analyze for distribution (suppress output)
            dist_cmd = f"gmx analyze -f angles/ang_{i}.xvg -dist angles/distr_ang_{i}.xvg -xvg none -bw 1.0"
            subprocess.run(dist_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            
            # Write to report
            with open(report_path, 'a') as report_file:
                report_file.write(f"{i}: {bead_trio}\n")
        
        # Clean up temp file
        os.remove(temp_ndx)
    
    # Clean up cleaned index file
    os.unlink(cleaned_angles)
    
    # Copy results back to original output directory
    if os.path.exists("angles"):
        # Ensure the output directory exists
        os.makedirs(angles_output_dir, exist_ok=True)
        # Copy each file individually
        for item in os.listdir("angles"):
            src = os.path.join("angles", item)
            dst = os.path.join(angles_output_dir, item)
            if os.path.isfile(src):
                shutil.copy2(src, dst)
            elif os.path.isdir(src):
                dest_subdir = os.path.join(angles_output_dir, item)
                if os.path.exists(dest_subdir):
                    shutil.rmtree(dest_subdir)
                shutil.copytree(src, dest_subdir)
    
    # Return to original directory and clean up
    os.chdir(original_dir)
    shutil.rmtree(temp_work_dir)
    
    print(f"Angles processing completed. Results in {angles_output_dir}")

def process_dihedrals(dihedrals_ndx, xtc_file, output_base_dir):
    """Process dihedrals"""
    print("\n=== Processing Dihedrals ===")
    
    # Convert to absolute paths
    dihedrals_ndx_abs = os.path.abspath(dihedrals_ndx)
    xtc_file_abs = os.path.abspath(xtc_file)
    
    # Create a temporary directory for the entire dihedral processing
    temp_work_dir = tempfile.mkdtemp(prefix='dihedrals_work_')
    original_dir = os.getcwd()
    
    # Change to temp directory
    os.chdir(temp_work_dir)
    
    # Get the output directory path for dihedrals
    dihedrals_output_dir = setup_directories("dihedrals", output_base_dir)
    
    # Create dihedrals directory in temp
    os.makedirs("dihedrals", exist_ok=True)
    
    # Clean the dihedrals index file
    cleaned_dihedrals = clean_index_file(dihedrals_ndx_abs)
    
    # Extract dihedrals (lines with exactly 4 numbers)
    dihedrals = extract_indices(cleaned_dihedrals, r'^\s*(\d+\s+\d+\s+\d+\s+\d+)\s*$')
    
    print(f"Found {len(dihedrals)} dihedrals")
    
    report_path = os.path.join("dihedrals", "report_dihedrals.txt")
    data_path = os.path.join("dihedrals", "data_dihedrals.txt")
    errors_path = os.path.join("dihedrals", "errors.log")
    
    Path(report_path).touch()
    Path(data_path).touch()
    
    for i, dihedral in enumerate(dihedrals):
        bead_quartet = f"{dihedral[0]} {dihedral[1]} {dihedral[2]} {dihedral[3]}"
        print(f"Processing dihedral {i}: {bead_quartet}")
        
        # Create temporary index file
        temp_ndx = f"temp_dih_{i}.ndx"
        with open(temp_ndx, 'w') as f:
            f.write(f"[ dihedral_{i} ]\n")
            f.write(f"{bead_quartet}\n")
        
        # Run gmx angle for dihedrals
        cmd = f"echo 0 | gmx angle -type dihedral -f {xtc_file_abs} -n {temp_ndx} -ov dihedrals/dih_{i}.xvg"
        log_path = os.path.join("dihedrals", f"dih_{i}.log")
        
        with open(log_path, 'w') as log_file:
            result = subprocess.run(cmd, shell=True, stdout=log_file, stderr=subprocess.STDOUT)
        
        dih_xvg = os.path.join("dihedrals", f"dih_{i}.xvg")
        
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
            
            # Run gmx analyze for distribution (suppress output)
            dist_cmd = f"gmx analyze -f dihedrals/dih_{i}.xvg -dist dihedrals/distr_dih_{i}.xvg -xvg none -bw 1.0"
            subprocess.run(dist_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            
            # Write to report
            with open(report_path, 'a') as report_file:
                report_file.write(f"{i}: {bead_quartet}\n")
        
        # Clean up temp file
        os.remove(temp_ndx)
    
    # Clean up cleaned index file
    os.unlink(cleaned_dihedrals)
    
    # Copy results back to original output directory
    if os.path.exists("dihedrals"):
        # Ensure the output directory exists
        os.makedirs(dihedrals_output_dir, exist_ok=True)
        # Copy each file individually
        for item in os.listdir("dihedrals"):
            src = os.path.join("dihedrals", item)
            dst = os.path.join(dihedrals_output_dir, item)
            if os.path.isfile(src):
                shutil.copy2(src, dst)
            elif os.path.isdir(src):
                dest_subdir = os.path.join(dihedrals_output_dir, item)
                if os.path.exists(dest_subdir):
                    shutil.rmtree(dest_subdir)
                shutil.copytree(src, dest_subdir)
    
    # Return to original directory and clean up
    os.chdir(original_dir)
    shutil.rmtree(temp_work_dir)
    
    print(f"Dihedrals processing completed. Results in {dihedrals_output_dir}")

def main():
    # Clean up old backup files before starting
    cleanup_backup_files()
    
    parser = argparse.ArgumentParser(
        description='Calculate bonds, angles, and dihedrals from molecular dynamics trajectories.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Compute bond distances, angles, and dihedral distributions from molecular
dynamics trajectories using GROMACS tools.

This script reads index files defining bonds (2 atoms), angles (3 atoms),
and dihedrals (4 atoms), and extracts structural properties from a trajectory
(.xtc) using the corresponding topology (.tpr). Results include averages,
standard deviations, and distributions for each interaction.

Optionally, the trajectory can be preprocessed to remove periodic boundary
conditions (PBC) and perform structural alignment (rotation/translation fitting),
which is recommended for non-centered or non-aligned systems.

IMPORTANT:
    The trajectory (.xtc) and topology (.tpr) must correspond to the same system
    and contain the same number of atoms.

Usage:
    python3 generate_bonds_angles_dihedrals.py \\
        --bonds_ndx bonds.ndx \\
        --angles_ndx angles.ndx \\
        --dihedrals_ndx dihedrals.ndx \\
        --xtc_file md.xtc \\
        --tpr_file md.tpr \\
        --output_all_files "OUTPUT"

Optional preprocessing (PBC removal and alignment):
    --index index.ndx        Index file for group selection (recommended)
    --remove_pbc             Enable PBC removal and alignment
    --group_1 "GroupName"    Group used for fitting (required if --remove_pbc)
    --group_2 "GroupName"    Group for output (default: "System")
    --keep_intermediate      Keep intermediate trajectory files

Arguments:
    --bonds_ndx         Index file defining bonds (pairs of atoms)
    --angles_ndx        Index file defining angles (triplets of atoms)
    --dihedrals_ndx     Index file defining dihedrals (quartets of atoms)
    --xtc_file          Input trajectory file (.xtc)
    --tpr_file          Input topology file (.tpr)
    --output_all_files  Base output directory path (creates bonds/, angles/, dihedrals/ subdirectories)

Output:
    The output directory specified by --output_all_files will contain three subdirectories:
        bonds/       Bond distance data and distributions
        angles/      Angle data and distributions
        dihedrals/   Dihedral data and distributions

    Each directory contains:
        - report_*.txt   Mapping of indices
        - data_*.txt     Average and standard deviation values
        - distr_*.xvg    Distribution files
        - errors.log     Execution errors (if any)

Notes:
    - Temporary files are automatically cleaned unless --keep_intermediate is used.
    - GROMACS (gmx) must be available in the system PATH.
    - For accurate structural analysis, preprocessing with --remove_pbc is recommended
      when dealing with periodic systems.
    - The output directory structure: OUTPUT/bonds/, OUTPUT/angles/, OUTPUT/dihedrals/
"""
    )
    
    # Required arguments
    parser.add_argument('--bonds_ndx', required=True, help='Bonds index file')
    parser.add_argument('--angles_ndx', required=True, help='Angles index file')
    parser.add_argument('--dihedrals_ndx', required=True, help='Dihedrals index file')
    parser.add_argument('--xtc_file', required=True, help='XTC trajectory file')
    parser.add_argument('--tpr_file', required=True, help='TPR topology file')
    parser.add_argument('--output_all_files', required=True, help='Base output directory path (creates bonds/, angles/, dihedrals/ subdirectories)')
    
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
    
    # Create the base output directory if it doesn't exist
    os.makedirs(args.output_all_files, exist_ok=True)
    
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
    process_bonds(args.bonds_ndx, xtc_file_to_use, args.tpr_file, args.output_all_files)
    process_angles(args.angles_ndx, xtc_file_to_use, args.output_all_files)
    process_dihedrals(args.dihedrals_ndx, xtc_file_to_use, args.output_all_files)
    
    # Final cleanup of any remaining backup files
    cleanup_backup_files()
    
    print("\n" + "="*50)
    print("Processing completed successfully!")
    print(f"Results are located in the following directories under {args.output_all_files}:")
    print(f"  - {os.path.join(args.output_all_files, 'bonds/')}")
    print(f"  - {os.path.join(args.output_all_files, 'angles/')}")
    print(f"  - {os.path.join(args.output_all_files, 'dihedrals/')}")
    print("\nCheck errors.log in each directory if any warnings or errors occurred.")
    print("="*50)

if __name__ == "__main__":
    main()
