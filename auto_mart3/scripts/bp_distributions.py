"""
Compute statistical properties (mean and standard deviation) of bond,
angle, and dihedral distributions from GROMACS XVG output files.

This script processes .xvg files generated from structural analyses
(e.g., gmx distance, gmx angle, gmx analyze) and extracts statistics
either from time series data or from precomputed distributions.

For each interaction type (bonds, angles, dihedrals), the script:
    - Reads corresponding XVG files
    - Detects whether the data is a time series or a distribution
    - Computes mean and standard deviation
    - Saves results as tab-separated values (.tsv)

Usage:
    python bp_distributions.py \
        --bonds_dir bonds/ \
        --angles_dir angles/ \
        --dihedrals_dir dihedrals/ \
        --dir_to_output results/

Arguments:
    --bonds_dir        Directory containing bond XVG files (default: bonds/)
    --angles_dir       Directory containing angle XVG files (default: angles/)
    --dihedrals_dir    Directory containing dihedral XVG files (default: dihedrals/)
    
    --dir_to_output    Output directory for TSV files (default: TSV_statistics/)
    --bond_out         Output filename for bond statistics (default: bond_statistics.tsv)
    --angle_out        Output filename for angle statistics (default: angle_statistics.tsv)
    --dihedral_out     Output filename for dihedral statistics (default: dihedral_statistics.tsv)

Output:
    Three TSV files are generated:
        - bond_statistics.tsv
        - angle_statistics.tsv
        - dihedral_statistics.tsv

    Each file contains:
        index    Identifier of the interaction
        mean     Average value
        sd       Standard deviation

Notes:
    - The script supports both time series data (value vs time) and
      distribution data (histograms).
    - Distribution files (e.g., distr_*.xvg) are used when available.
    - XVG comment lines (starting with '#' or '@') are ignored.
    - Requires NumPy and pandas.
"""

import os
import glob
import numpy as np
import pandas as pd
from scipy import stats
import re
import argparse

# ============================================
# Helper functions
# ============================================

def extract_index(filename):
    return int(re.search(r'\d+', os.path.basename(filename)).group())


def read_xvg(filename):
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(('#', '@')):
                continue
            if line.strip():
                try:
                    values = line.strip().split()
                    if len(values) == 2:
                        data.append([float(values[0]), float(values[1])])
                    elif len(values) == 1:
                        data.append([float(values[0])])
                except ValueError:
                    continue
    return np.array(data)


def calculate_statistics(data):
    """Calculate mean and standard deviation from data"""
    if data.shape[1] >= 2:
        # Time series data
        values = data[:, 1]
        mean_val = np.mean(values)
        std_val = np.std(values)
    else:
        # Distribution data
        bins = data[:, 0]
        density = data[:, 1] if data.shape[1] > 1 else None
        
        if density is not None:
            mean_val = np.sum(bins * density) / np.sum(density) if np.sum(density) > 0 else 0
            var_val = np.sum(((bins - mean_val)**2) * density) / np.sum(density) if np.sum(density) > 0 else 0
            std_val = np.sqrt(var_val)
        else:
            mean_val = np.mean(bins)
            std_val = np.std(bins)
    
    return mean_val, std_val


# ============================================
# Main function
# ============================================

def main():
    # ============================================
    # Argument parser
    # ============================================
    
    parser = argparse.ArgumentParser(description="Process XVG files (bonds, angles, dihedrals)")
    
    # Inputs
    parser.add_argument("--bonds_dir", default="bonds")
    parser.add_argument("--angles_dir", default="angles")
    parser.add_argument("--dihedrals_dir", default="dihedrals")
    
    # Outputs (tables)
    parser.add_argument("--dir_to_output", default="TSV_statistics", 
                        help="Directory to save *_statistics.tsv files (optional)")
    parser.add_argument("--bond_out", default="bond_statistics.tsv")
    parser.add_argument("--angle_out", default="angle_statistics.tsv")
    parser.add_argument("--dihedral_out", default="dihedral_statistics.tsv")
    
    args = parser.parse_args()
    
    
    # ============================================
    # Create output directory if it doesn't exist
    # ============================================
    
    if args.dir_to_output:
        os.makedirs(args.dir_to_output, exist_ok=True)
    
    
    # ============================================
    # Process Bonds
    # ============================================
    
    print("Processing Bonds...")
    bond_stats = []
    
    bond_files = sorted(
        glob.glob(os.path.join(args.bonds_dir, "bond_*.xvg")),
        key=extract_index
    )
    
    for bond_file in bond_files:
        idx = int(os.path.basename(bond_file).replace('bond_', '').replace('.xvg', ''))
        data = read_xvg(bond_file)
    
        if len(data) > 0:
            if data.shape[1] == 2 and np.std(data[:, 0]) > 0.1:
                # Time series data
                mean_val, std_val = calculate_statistics(data)
                
                bond_stats.append({
                    'index': idx,
                    'mean': mean_val,
                    'sd': std_val,
                    'type': 'bond'
                })
                
                print(f"  Bond {idx}: mean={mean_val:.4f}, sd={std_val:.4f}")
                
            else:
                # Distribution data
                dist_file = os.path.join(args.bonds_dir, f"distr_bond_{idx}.xvg")
                
                if os.path.exists(dist_file):
                    dist_data = read_xvg(dist_file)
                    
                    if len(dist_data) > 0:
                        mean_val, std_val = calculate_statistics(dist_data)
                        
                        bond_stats.append({
                            'index': idx,
                            'mean': mean_val,
                            'sd': std_val,
                            'type': 'bond'
                        })
                        
                        print(f"  Bond {idx} (dist): mean={mean_val:.4f}, sd={std_val:.4f}")
    
    if bond_stats:
        df_bonds = pd.DataFrame(bond_stats)
        bond_output_path = os.path.join(args.dir_to_output, args.bond_out)
        df_bonds.to_csv(bond_output_path, sep='\t', index=False)
        print(f"Saved {bond_output_path} with {len(bond_stats)} bonds")
    
    
    # ============================================
    # Process Angles
    # ============================================
    
    print("\nProcessing Angles...")
    angle_stats = []
    
    angle_files = sorted(
        glob.glob(os.path.join(args.angles_dir, "ang_*.xvg")),
        key=extract_index
    )
    
    for angle_file in angle_files:
        idx = int(os.path.basename(angle_file).replace('ang_', '').replace('.xvg', ''))
        data = read_xvg(angle_file)
    
        if len(data) > 0:
            if data.shape[1] == 2 and np.std(data[:, 0]) > 0.1:
                # Time series data
                mean_val, std_val = calculate_statistics(data)
                
                angle_stats.append({
                    'index': idx,
                    'mean': mean_val,
                    'sd': std_val,
                    'type': 'angle'
                })
                
                print(f"  Angle {idx}: mean={mean_val:.4f}, sd={std_val:.4f}")
                
            else:
                # Distribution data
                dist_file = os.path.join(args.angles_dir, f"distr_ang_{idx}.xvg")
                
                if os.path.exists(dist_file):
                    dist_data = read_xvg(dist_file)
                    
                    if len(dist_data) > 0:
                        mean_val, std_val = calculate_statistics(dist_data)
                        
                        angle_stats.append({
                            'index': idx,
                            'mean': mean_val,
                            'sd': std_val,
                            'type': 'angle'
                        })
                        
                        print(f"  Angle {idx} (dist): mean={mean_val:.4f}, sd={std_val:.4f}")
    
    if angle_stats:
        df_angles = pd.DataFrame(angle_stats)
        angle_output_path = os.path.join(args.dir_to_output, args.angle_out)
        df_angles.to_csv(angle_output_path, sep='\t', index=False)
        print(f"Saved {angle_output_path} with {len(angle_stats)} angles")
    
    
    # ============================================
    # Process Dihedrals
    # ============================================
    
    print("\nProcessing Dihedrals...")
    dihedral_stats = []
    
    dihedral_files = sorted(
        glob.glob(os.path.join(args.dihedrals_dir, "dih_*.xvg")),
        key=extract_index
    )
    
    for dihedral_file in dihedral_files:
        idx = int(os.path.basename(dihedral_file).replace('dih_', '').replace('.xvg', ''))
        data = read_xvg(dihedral_file)
    
        if len(data) > 0:
            if data.shape[1] == 2 and np.std(data[:, 0]) > 0.1:
                # Time series data
                mean_val, std_val = calculate_statistics(data)
                
                dihedral_stats.append({
                    'index': idx,
                    'mean': mean_val,
                    'sd': std_val
                })
                
                print(f"  Dihedral {idx}: mean={mean_val:.4f}, sd={std_val:.4f}")
                
            else:
                # Distribution data
                dist_file = os.path.join(args.dihedrals_dir, f"distr_dih_{idx}.xvg")
                
                if os.path.exists(dist_file):
                    dist_data = read_xvg(dist_file)
                    
                    if len(dist_data) > 0:
                        mean_val, std_val = calculate_statistics(dist_data)
                        
                        dihedral_stats.append({
                            'index': idx,
                            'mean': mean_val,
                            'sd': std_val
                        })
                        
                        print(f"  Dihedral {idx} (dist): mean={mean_val:.4f}, sd={std_val:.4f}")
    
    if dihedral_stats:
        df_dihedrals = pd.DataFrame(dihedral_stats)
        dihedral_output_path = os.path.join(args.dir_to_output, args.dihedral_out)
        df_dihedrals.to_csv(
            dihedral_output_path,
            sep='\t',
            index=False,
            columns=['index', 'mean', 'sd']
        )
        
        print(f"Saved {dihedral_output_path} with {len(dihedral_stats)} dihedrals")
    
    
    # ============================================
    # Final summary
    # ============================================
    
    print("\nProcessing completed!")
    
    print("Generated files:")
    print(f"  - {os.path.join(args.dir_to_output, args.bond_out)}")
    print(f"  - {os.path.join(args.dir_to_output, args.angle_out)}")
    print(f"  - {os.path.join(args.dir_to_output, args.dihedral_out)}")


# ============================================
# Script entry point
# ============================================

if __name__ == "__main__":
    main()
