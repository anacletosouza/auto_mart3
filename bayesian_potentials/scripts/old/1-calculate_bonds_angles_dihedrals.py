#!/usr/bin/env python3
"""
Script to analyze bonds, angles, and dihedrals from GROMACS trajectories.
Equivalent to the original bash script but with -- arguments.
"""

import argparse
import subprocess
import os
import sys
import re
from pathlib import Path

def parse_arguments():
    """Parse command line arguments with -- format."""
    parser = argparse.ArgumentParser(
        description='Analyze bonds, angles, and dihedrals from GROMACS trajectory'
    )
    parser.add_argument('--bonds_ndx', required=True, 
                       help='Index file with bond definitions')
    parser.add_argument('--angles_ndx', required=True,
                       help='Index file with angle definitions')
    parser.add_argument('--dihedrals_ndx', required=True,
                       help='Index file with dihedral definitions')
    parser.add_argument('--xtc_file', required=True,
                       help='Input trajectory file (XTC format)')
    parser.add_argument('--tpr_file', required=True,
                       help='Input run input file (TPR format)')
    parser.add_argument('--out_prefix', required=True,
                       help='Prefix for output directories')
    parser.add_argument('--index', default=None,
                       help='Optional custom index file (index.ndx)')
    
    return parser.parse_args()

def run_gmx_command(cmd, log_file, error_log, description):
    """Run a GROMACS command and handle errors."""
    try:
        # For gmx commands that expect input, provide '0' (for selecting group)
        if cmd[0] == 'gmx' and cmd[1] in ['distance', 'angle']:
            process = subprocess.Popen(cmd, stdin=subprocess.PIPE, 
                                      stdout=subprocess.PIPE, 
                                      stderr=subprocess.PIPE,
                                      universal_newlines=True)
            stdout, stderr = process.communicate(input='0\n')
            returncode = process.returncode
        else:
            result = subprocess.run(cmd, capture_output=True, text=True)
            stdout, stderr = result.stdout, result.stderr
            returncode = result.returncode
        
        # Write output to log file
        with open(log_file, 'w') as f:
            f.write(stdout)
            f.write(stderr)
        
        if returncode != 0:
            with open(error_log, 'a') as f:
                f.write(f"ERROR on {description}\n")
                f.write(f"Command: {' '.join(cmd)}\n")
                f.write(stderr)
                f.write("\n")
            return False
        return True
        
    except Exception as e:
        with open(error_log, 'a') as f:
            f.write(f"EXCEPTION on {description}: {str(e)}\n")
        return False

def process_bonds(args):
    """Process all bonds."""
    print("\n=== Processing Bonds ===")
    
    dir_name = f"{args.out_prefix}_bonds"
    if os.path.exists(dir_name):
        subprocess.run(['rm', '-rf', dir_name])
    os.makedirs(dir_name)
    
    # Extract bonds from index file
    bonds_list = []
    with open(args.bonds_ndx, 'r') as f:
        for line in f:
            if re.match(r'^\s*\d+\s+\d+', line):
                bonds_list.append(line.strip())
    
    nbonds = len(bonds_list)
    print(f"Found {nbonds} bonds")
    
    report_file = os.path.join(dir_name, "report_bonds.txt")
    with open(report_file, 'w') as f:
        f.write("")
    
    data_file = os.path.join(dir_name, "data_bonds.txt")
    error_log = os.path.join(dir_name, "errors.log")
    
    for ibond, bead_pair in enumerate(bonds_list):
        print(f"Processing bond {ibond}: {bead_pair}")
        
        # Create temporary index file
        temp_ndx = f"temp_bond_{ibond}.ndx"
        with open(temp_ndx, 'w') as f:
            f.write(f"[ bond_{ibond} ]\n")
            f.write(f"{bead_pair}\n")
        
        # Run gmx distance
        xvg_file = os.path.join(dir_name, f"bond_{ibond}.xvg")
        log_file = os.path.join(dir_name, f"bond_{ibond}.log")
        
        cmd = ['gmx', 'distance', '-f', args.xtc_file, 
               '-n', temp_ndx, '-s', args.tpr_file,
               '-oall', xvg_file, '-xvg', 'none']
        
        if run_gmx_command(cmd, log_file, error_log, f"bond {ibond}"):
            # Extract average and standard deviation
            with open(log_file, 'r') as f:
                content = f.read()
                avg_match = re.search(r'Average distance\s+([\d.]+)', content)
                std_match = re.search(r'Standard deviation\s+([\d.]+)', content)
                
                if avg_match and std_match:
                    with open(data_file, 'a') as f:
                        f.write(f"---- bond {ibond} ({bead_pair}) ----\n")
                        f.write(f"{avg_match.group(1)}\n")
                        f.write(f"{std_match.group(1)}\n")
            
            # Run gmx analyze for distribution
            distr_file = os.path.join(dir_name, f"distr_bond_{ibond}.xvg")
            analyze_cmd = ['gmx', 'analyze', '-f', xvg_file, 
                          '-dist', distr_file, '-xvg', 'none', 
                          '-bw', '0.001']
            subprocess.run(analyze_cmd, capture_output=True, text=True)
            
            # Update report
            with open(report_file, 'a') as f:
                f.write(f"{ibond}: {bead_pair}\n")
        
        # Cleanup
        os.remove(temp_ndx)
    
    return True

def process_angles(args):
    """Process all angles."""
    print("\n=== Processing Angles ===")
    
    dir_name = f"{args.out_prefix}_angles"
    if os.path.exists(dir_name):
        subprocess.run(['rm', '-rf', dir_name])
    os.makedirs(dir_name)
    
    # Extract angles from index file
    angles_list = []
    with open(args.angles_ndx, 'r') as f:
        for line in f:
            if re.match(r'^\s*\d+\s+\d+\s+\d+', line):
                angles_list.append(line.strip())
    
    nangles = len(angles_list)
    print(f"Found {nangles} angles")
    
    report_file = os.path.join(dir_name, "report_angles.txt")
    with open(report_file, 'w') as f:
        f.write("")
    
    data_file = os.path.join(dir_name, "data_angles.txt")
    error_log = os.path.join(dir_name, "errors.log")
    
    for iang, bead_trio in enumerate(angles_list):
        print(f"Processing angle {iang}: {bead_trio}")
        
        # Create temporary index file
        temp_ndx = f"temp_angle_{iang}.ndx"
        with open(temp_ndx, 'w') as f:
            f.write(f"[ angle_{iang} ]\n")
            f.write(f"{bead_trio}\n")
        
        # Run gmx angle
        xvg_file = os.path.join(dir_name, f"ang_{iang}.xvg")
        log_file = os.path.join(dir_name, f"ang_{iang}.log")
        
        cmd = ['gmx', 'angle', '-f', args.xtc_file, 
               '-n', temp_ndx, '-ov', xvg_file]
        
        if run_gmx_command(cmd, log_file, error_log, f"angle {iang}"):
            # Extract angle and standard deviation
            with open(log_file, 'r') as f:
                content = f.read()
                angle_match = re.search(r'<\s*angle\s*>\s+([\d.]+)', content)
                std_match = re.search(r'Std\. Dev\.\s+([\d.]+)', content)
                
                if angle_match and std_match:
                    with open(data_file, 'a') as f:
                        f.write(f"---- ang {iang} ({bead_trio}) ----\n")
                        f.write(f"{angle_match.group(1)}\n")
                        f.write(f"{std_match.group(1)}\n")
            
            # Run gmx analyze for distribution
            distr_file = os.path.join(dir_name, f"distr_ang_{iang}.xvg")
            analyze_cmd = ['gmx', 'analyze', '-f', xvg_file, 
                          '-dist', distr_file, '-xvg', 'none', 
                          '-bw', '1.0']
            subprocess.run(analyze_cmd, capture_output=True, text=True)
            
            # Update report
            with open(report_file, 'a') as f:
                f.write(f"{iang}: {bead_trio}\n")
        
        # Cleanup
        os.remove(temp_ndx)
    
    return True

def process_dihedrals(args):
    """Process all dihedrals."""
    print("\n=== Processing Dihedrals ===")
    
    dir_name = f"{args.out_prefix}_dihedrals"
    if os.path.exists(dir_name):
        subprocess.run(['rm', '-rf', dir_name])
    os.makedirs(dir_name)
    
    # Extract dihedrals from index file
    dihedrals_list = []
    with open(args.dihedrals_ndx, 'r') as f:
        for line in f:
            if re.match(r'^\s*\d+\s+\d+\s+\d+\s+\d+', line):
                dihedrals_list.append(line.strip())
    
    ndihedrals = len(dihedrals_list)
    print(f"Found {ndihedrals} dihedrals")
    
    report_file = os.path.join(dir_name, "report_dihedrals.txt")
    with open(report_file, 'w') as f:
        f.write("")
    
    data_file = os.path.join(dir_name, "data_dihedrals.txt")
    error_log = os.path.join(dir_name, "errors.log")
    
    for idih, bead_quartet in enumerate(dihedrals_list):
        print(f"Processing dihedral {idih}: {bead_quartet}")
        
        # Validate bead numbers (1-30 as in original)
        valid = True
        for num in bead_quartet.split():
            if int(num) < 1 or int(num) > 30:
                valid = False
                with open(error_log, 'a') as f:
                    f.write(f"ERROR: bead {num} out of range (1-30)\n")
        
        if valid:
            # Create temporary index file
            temp_ndx = f"temp_dih_{idih}.ndx"
            with open(temp_ndx, 'w') as f:
                f.write(f"[ dihedral_{idih} ]\n")
                f.write(f"{bead_quartet}\n")
            
            # Run gmx angle for dihedral
            xvg_file = os.path.join(dir_name, f"dih_{idih}.xvg")
            log_file = os.path.join(dir_name, f"dih_{idih}.log")
            
            cmd = ['gmx', 'angle', '-type', 'dihedral', 
                   '-f', args.xtc_file, '-n', temp_ndx, 
                   '-ov', xvg_file]
            
            if run_gmx_command(cmd, log_file, error_log, f"dihedral {idih}"):
                # Extract angle and standard deviation
                with open(log_file, 'r') as f:
                    content = f.read()
                    angle_match = re.search(r'<\s*angle\s*>\s+([\d.]+)', content)
                    std_match = re.search(r'Std\. Dev\.\s+([\d.]+)', content)
                    
                    if angle_match and std_match:
                        with open(data_file, 'a') as f:
                            f.write(f"---- dih {idih} ({bead_quartet}) ----\n")
                            f.write(f"{angle_match.group(1)}\n")
                            f.write(f"{std_match.group(1)}\n")
                
                # Run gmx analyze for distribution
                distr_file = os.path.join(dir_name, f"distr_dih_{idih}.xvg")
                analyze_cmd = ['gmx', 'analyze', '-f', xvg_file, 
                              '-dist', distr_file, '-xvg', 'none', 
                              '-bw', '1.0']
                subprocess.run(analyze_cmd, capture_output=True, text=True)
                
                # Update report
                with open(report_file, 'a') as f:
                    f.write(f"{idih}: {bead_quartet}\n")
            
            # Cleanup
            os.remove(temp_ndx)
    
    return True

def main():
    """Main execution function."""
    args = parse_arguments()
    
    print("=" * 60)
    print("Starting GROMACS analysis")
    print("=" * 60)
    
    # Process bonds, angles, and dihedrals
    process_bonds(args)
    process_angles(args)
    process_dihedrals(args)
    
    print("\n" + "=" * 60)
    print("Processing completed!")
    print(f"Results in: {args.out_prefix}_bonds/, {args.out_prefix}_angles/, {args.out_prefix}_dihedrals/")
    print("Check errors.log in each directory if needed.")
    print("=" * 60)

if __name__ == "__main__":
    main()
