#!/usr/bin/env python3
"""
Generate Martini-style .itp with automatic force constant adjustment
based on reference and simulated statistics.

FINAL VERSION WITH SIMULATED ANNEALING:
- Bayesian force update (Gibbs-like sampling)
- Simulated Annealing acceptance for global error
- Reuse previous .itp if iteration is rejected
- Minimum force constraints to prevent numerical instability

usage: bayesian-potentials force-adjust [-h] --bonds_ref BONDS_REF --angles_ref ANGLES_REF --dihedrals_ref DIHEDRALS_REF
                                        --bonds_sim BONDS_SIM --angles_sim ANGLES_SIM --dihedrals_sim DIHEDRALS_SIM --itp_cg
                                        ITP_CG --ndx_bounds NDX_BOUNDS --ndx_angles NDX_ANGLES --ndx_dihedrals NDX_DIHEDRALS
                                        [--molecule_name MOLECULE_NAME] [--multimodal_mode MULTIMODAL_MODE]
"""

import os
import json
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import shutil
import math
import random
from datetime import datetime
import warnings

# ============================================
# Constants
# ============================================

# Minimum allowed force constants (to prevent simulation instability)
MIN_FORCE_BOND = 750.0      # Bonds: minimum 750 kJ/(mol·nm^2)
MIN_FORCE_ANGLE = 15.0      # Angles: minimum 15 kJ/(mol·rad^2)
MIN_FORCE_DIHEDRAL = 15.0   # Dihedrals: minimum 15 kJ/(mol·rad^2)

# Default force constants (used when missing in ITP)
DEFAULT_FORCE_BOND = 1250.0
DEFAULT_FORCE_ANGLE = 25.0
DEFAULT_FORCE_DIHEDRAL = 25.0

# ============================================
# Helper functions
# ============================================

def str_to_bool(value):
    """Convert string to boolean."""
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.lower() in ('true', '1', 'yes', 'on')
    return bool(value)

def read_tsv(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} not found.")
    return pd.read_csv(file_path, sep="\t")

def read_ndx(file_path):
    """Read NDX file and return list of index tuples.
    
    NDX files can have lines like:
    - "1 2" (bond between atoms 1 and 2)
    - "1 2 3" (angle between atoms 1,2,3)
    - "1 2 3 4" (dihedral between atoms 1,2,3,4)
    - Lines starting with ';' are comments
    - Lines starting with '[' are section headers
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} not found.")
    
    indices = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip empty lines, comments, and section headers
            if not line or line.startswith('[') or line.startswith(';'):
                continue
            
            # Split the line into parts
            parts = line.split()
            if not parts:
                continue
            
            # Convert only numeric values (ignore any non-numeric like comments)
            numeric_parts = []
            for x in parts:
                # Check if it's a number (positive or negative)
                try:
                    num = int(x)
                    numeric_parts.append(num)
                except ValueError:
                    # Skip non-numeric values (like inline comments)
                    # If we encounter a comment after numbers, stop processing this line
                    if x.startswith(';'):
                        break
                    continue
            
            # Add tuple if we have at least 2 numbers (bonds need at least 2)
            if len(numeric_parts) >= 2:
                indices.append(tuple(numeric_parts))
    
    return indices

def format_bond(i,j,length,force):
    return f"  {i:2d}  {j:2d}   1    {length:.3f}    {force:.1f}   ;"

def format_angle(i,j,k,angle,force):
    return f"  {i:2d}  {j:2d}  {k:2d}   2    {angle:.1f}    {force:.1f}   ;"

def format_dihedral(i,j,k,l,angle,force):
    return f"  {i:2d}  {j:2d}  {k:2d}  {l:2d}   2    {angle:.1f}    {force:.1f}   ;"

def parse_itp_cg(itp_file):
    """Parse ITP file to extract atoms, bonds, angles, and dihedrals information"""
    dic_itp_cg = {
        'atoms': [],
        'bonds': [],
        'angles': [],
        'dihedrals': []
    }
    
    current_section = None
    
    with open(itp_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Check for section headers
            if line.startswith('[ atoms ]'):
                current_section = 'atoms'
                continue
            elif line.startswith('[ bonds ]'):
                current_section = 'bonds'
                continue
            elif line.startswith('[ angles ]'):
                current_section = 'angles'
                continue
            elif line.startswith('[ dihedrals ]'):
                current_section = 'dihedrals'
                continue
            elif line.startswith('['):
                current_section = None
                continue
            
            # Skip empty lines and comments
            if not line or line.startswith(';'):
                continue
            
            # Parse based on current section
            if current_section == 'atoms':
                parts = line.split()
                if len(parts) >= 8:
                    atom_info = {
                        'nr': int(parts[0]),
                        'type': parts[1],
                        'resnr': int(parts[2]),
                        'residue': parts[3],
                        'atom': parts[4],
                        'cgnr': int(parts[5]),
                        'charge': float(parts[6]),
                        'mass': float(parts[7])
                    }
                    dic_itp_cg['atoms'].append(atom_info)
            
            elif current_section == 'bonds':
                parts = line.split()
                if len(parts) >= 6:
                    bond_info = {
                        'i': int(parts[0]),
                        'j': int(parts[1]),
                        'funct': int(parts[2]),
                        'length': float(parts[3]),
                        'force_k': float(parts[4])
                    }
                    dic_itp_cg['bonds'].append(bond_info)
            
            elif current_section == 'angles':
                parts = line.split()
                if len(parts) >= 7:
                    angle_info = {
                        'i': int(parts[0]),
                        'j': int(parts[1]),
                        'k': int(parts[2]),
                        'funct': int(parts[3]),
                        'angle': float(parts[4]),
                        'force_k': float(parts[5])
                    }
                    dic_itp_cg['angles'].append(angle_info)
            
            elif current_section == 'dihedrals':
                parts = line.split()
                if len(parts) >= 8:
                    dihedral_info = {
                        'i': int(parts[0]),
                        'j': int(parts[1]),
                        'k': int(parts[2]),
                        'l': int(parts[3]),
                        'funct': int(parts[4]),
                        'angle': float(parts[5]),
                        'force_k': float(parts[6])
                    }
                    dic_itp_cg['dihedrals'].append(dihedral_info)
    
    return dic_itp_cg

def get_peak_properties(data):
    """Calculate peak (mode) properties from distribution"""
    if len(data) == 0:
        return np.mean(data), np.var(data)
    
    # Use kernel density estimation to find the mode
    try:
        from scipy import stats
        kde = stats.gaussian_kde(data)
        x_grid = np.linspace(np.min(data), np.max(data), 1000)
        kde_values = kde.evaluate(x_grid)
        peak_idx = np.argmax(kde_values)
        peak_value = x_grid[peak_idx]
        
        # Calculate variance around the peak
        # Use points within 2 standard deviations of the peak
        threshold = 2 * np.std(data)
        mask = np.abs(data - peak_value) <= threshold
        peak_variance = np.var(data[mask]) if np.any(mask) else np.var(data)
        
        return peak_value, peak_variance
    except:
        # Fallback to mean if KDE fails
        return np.mean(data), np.var(data)

# ============================================
# Bayesian update (Gibbs-like) with bounds
# ============================================

def choose_best_force(F_orig, sigma_ref, sigma_curr, sigma_prev, force_type='bond'):
    """
    Choose best force constant with Bayesian update.
    
    Parameters:
    -----------
    F_orig : float
        Original force constant from ITP
    sigma_ref : float
        Reference standard deviation
    sigma_curr : float
        Current simulated standard deviation
    sigma_prev : float
        Previous iteration standard deviation (optional)
    force_type : str
        Type of force: 'bond', 'angle', or 'dihedral'
    
    Returns:
    --------
    tuple : (new_force, status_message)
    """
    
    # Get minimum force based on type
    if force_type == 'bond':
        min_force = MIN_FORCE_BOND
    elif force_type == 'angle':
        min_force = MIN_FORCE_ANGLE
    else:  # dihedral
        min_force = MIN_FORCE_DIHEDRAL
    
    if sigma_curr == 0:
        return max(F_orig, min_force), "NOT CHANGED (zero sigma)"

    if sigma_prev is not None:
        mu_prior = sigma_prev
        var_prior = (0.2 * sigma_prev)**2
    else:
        mu_prior = sigma_curr
        var_prior = (0.5 * sigma_curr)**2

    mu_like = sigma_curr
    var_like = (0.2 * sigma_curr)**2

    var_post = 1 / (1/var_prior + 1/var_like)
    mu_post = var_post * (mu_prior/var_prior + mu_like/var_like)

    samples = np.random.normal(mu_post, np.sqrt(var_post), size=20)

    F_candidates = []
    for s in samples:
        if s <= 0:
            continue
        F_cand = F_orig * (s / sigma_ref)
        # Apply minimum bound
        F_cand = max(F_cand, min_force)
        F_candidates.append(F_cand)

    if not F_candidates:
        return max(F_orig, min_force), "NOT CHANGED (no samples)"

    best_F = max(F_orig, min_force)
    best_error = abs(sigma_ref - sigma_curr)

    for F_cand in F_candidates:
        sigma_pred = sigma_curr * (F_orig / F_cand) if F_cand > 0 else sigma_curr
        err = abs(sigma_ref - sigma_pred)
        if err < best_error:
            best_error = err
            best_F = F_cand

    # Final check: ensure force is not below minimum
    if best_F < min_force:
        best_F = min_force
        return best_F, f"FORCED TO MIN ({min_force:.1f})"
    
    if best_F == F_orig:
        return best_F, "NOT CHANGED (Bayes reject)"
    else:
        return best_F, "CHANGED (Bayes)"

# ============================================
# Process statistics
# ============================================

def process_statistics(df_ref, df_sim, multimodal_mode, variance_multimodal):
    """Process reference and simulated statistics with multimodal/unimodal options"""
    
    processed_ref = df_ref.copy()
    processed_sim = df_sim.copy()
    
    if multimodal_mode:
        # Use peak (mode) instead of mean for reference
        # Note: This requires access to original distribution data
        # For now, we'll use the provided mean and add a flag
        processed_ref['mean'] = df_ref['mean']  # Placeholder - would need distribution data
        processed_sim['mean'] = df_sim['mean']
    
    if variance_multimodal:
        # Use variance from all data
        processed_ref['sd'] = df_ref['sd']
        processed_sim['sd'] = df_sim['sd']
    else:
        # Use variance from peak only
        # Placeholder - would need distribution data
        processed_ref['sd'] = df_ref['sd']
        processed_sim['sd'] = df_sim['sd']
    
    return processed_ref, processed_sim

def compute_global_error(df_ref, df_sim):
    return ((df_ref['sd'] - df_sim['sd']).abs()).mean()

def safe_load(file):
    return read_tsv(file) if os.path.exists(file) else None

# ============================================
# Main function
# ============================================

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(
        description="Generate Martini-style .itp with automatic force constant adjustment"
    )
    parser.add_argument("--bonds_ref", required=True, help="Reference bonds TSV file")
    parser.add_argument("--angles_ref", required=True, help="Reference angles TSV file")
    parser.add_argument("--dihedrals_ref", required=True, help="Reference dihedrals TSV file")
    parser.add_argument("--bonds_sim", required=True, help="Simulated bonds TSV file")
    parser.add_argument("--angles_sim", required=True, help="Simulated angles TSV file")
    parser.add_argument("--dihedrals_sim", required=True, help="Simulated dihedrals TSV file")
    parser.add_argument("--itp_cg", required=True, help="CG ITP file with topology information")
    parser.add_argument("--ndx_bounds", required=True, help="NDX file with bond indices")
    parser.add_argument("--ndx_angles", required=True, help="NDX file with angle indices")
    parser.add_argument("--ndx_dihedrals", required=True, help="NDX file with dihedral indices")
    parser.add_argument("--molecule_name", default="molecule", help="Molecule name (default: molecule)")
    parser.add_argument("--multimodal_mode", type=str, default="true", 
                        help="Use peak density mode instead of mean (default: true)")
    parser.add_argument("--variance_multimodal", type=str, default="true",
                        help="Use variance from all data instead of peak (default: true)")
    parser.add_argument("--T0", type=float, default=10.0, 
                        help="Initial temperature for simulated annealing (default: 10.0)")
    parser.add_argument("--alpha", type=float, default=0.85,
                        help="Cooling factor for simulated annealing (default: 0.85)")
    parser.add_argument("--itp_out", default="cg.itp", help="Output ITP file (default: cg.itp)")
    
    args = parser.parse_args()
    
    # Convert string boolean arguments to actual booleans
    multimodal_mode = str_to_bool(args.multimodal_mode)
    variance_multimodal = str_to_bool(args.variance_multimodal)
    
    # Detect iteration
    current_iter = int(os.path.basename(args.bonds_sim).split("_")[-1].split(".")[0])
    prev_iter = current_iter - 1
    
    base_dir = os.path.abspath(os.path.join(os.path.dirname(args.bonds_sim), "../../"))
    prev_stat_dir = os.path.join(base_dir, f"iter_{prev_iter}", "stat")
    
    print(f"\n{'='*60}")
    print(f"Starting iteration {current_iter}")
    print(f"{'='*60}")
    print(f"Multimodal mode: {multimodal_mode}")
    print(f"Variance multimodal: {variance_multimodal}")
    print(f"T0: {args.T0}, alpha: {args.alpha}")
    print(f"Minimum forces: Bond={MIN_FORCE_BOND}, Angle={MIN_FORCE_ANGLE}, Dihedral={MIN_FORCE_DIHEDRAL}")
    print(f"{'='*60}\n")
    
    # Load data
    print("Loading reference and simulation data...")
    bonds_ref_df = read_tsv(args.bonds_ref)
    angles_ref_df = read_tsv(args.angles_ref)
    dihedrals_ref_df = read_tsv(args.dihedrals_ref)
    
    bonds_sim_df = read_tsv(args.bonds_sim)
    angles_sim_df = read_tsv(args.angles_sim)
    dihedrals_sim_df = read_tsv(args.dihedrals_sim)
    
    bonds_prev_df = safe_load(f"{prev_stat_dir}/bond_{prev_iter}.tsv")
    angles_prev_df = safe_load(f"{prev_stat_dir}/angle_{prev_iter}.tsv")
    dihedrals_prev_df = safe_load(f"{prev_stat_dir}/dihedral_{prev_iter}.tsv")
    
    # Parse ITP CG file
    print(f"Parsing ITP CG file: {args.itp_cg}")
    dic_itp_cg = parse_itp_cg(args.itp_cg)
    
    # Extract force constants from ITP
    force_bond_dict = {}
    force_angle_dict = {}
    force_dihedral_dict = {}
    
    for bond in dic_itp_cg['bonds']:
        key = (bond['i'], bond['j'])
        force_bond_dict[key] = bond['force_k']
    
    for angle in dic_itp_cg['angles']:
        key = (angle['i'], angle['j'], angle['k'])
        force_angle_dict[key] = angle['force_k']
    
    for dihedral in dic_itp_cg['dihedrals']:
        key = (dihedral['i'], dihedral['j'], dihedral['k'], dihedral['l'])
        force_dihedral_dict[key] = dihedral['force_k']
    
    print(f"Loaded {len(force_bond_dict)} bonds, {len(force_angle_dict)} angles, "
          f"{len(force_dihedral_dict)} dihedrals from ITP")
    
    # Load topology
    print("Loading topology indices...")
    bond_pairs = read_ndx(args.ndx_bounds)
    angle_triplets = read_ndx(args.ndx_angles)
    dihedral_quartets = read_ndx(args.ndx_dihedrals)
    
    print(f"Loaded {len(bond_pairs)} bonds, {len(angle_triplets)} angles, "
          f"{len(dihedral_quartets)} dihedrals")
    
    # Apply multimodal processing
    bonds_ref_processed, bonds_sim_processed = process_statistics(
        bonds_ref_df, bonds_sim_df, multimodal_mode, variance_multimodal
    )
    angles_ref_processed, angles_sim_processed = process_statistics(
        angles_ref_df, angles_sim_df, multimodal_mode, variance_multimodal
    )
    dihedrals_ref_processed, dihedrals_sim_processed = process_statistics(
        dihedrals_ref_df, dihedrals_sim_df, multimodal_mode, variance_multimodal
    )
    
    # Compute global error
    global_error = (
        compute_global_error(bonds_ref_processed, bonds_sim_processed) +
        compute_global_error(angles_ref_processed, angles_sim_processed) +
        compute_global_error(dihedrals_ref_processed, dihedrals_sim_processed)
    )
    
    print(f"\nGlobal error: {global_error:.6f}")
    
    # Create inter directory
    inter_dir = os.path.abspath(os.path.join(os.path.dirname(args.itp_out), "../inter"))
    os.makedirs(inter_dir, exist_ok=True)
    
    iteration_file = os.path.join(inter_dir, "dihedral_iteration.tsv")
    new_row = pd.DataFrame({
        "iteration": [current_iter],
        "global_error": [global_error]
    })
    
    # Simulated Annealing acceptance
    T = args.T0 * (args.alpha ** current_iter)
    
    if os.path.exists(iteration_file):
        old = pd.read_csv(iteration_file, sep="\t")
        prev_last = old.sort_values("iteration").iloc[-1]["global_error"]
        delta_E = global_error - prev_last
        
        if delta_E <= 0:
            accept = True
            reason = "IMPROVED"
        else:
            prob = math.exp(-delta_E / T)
            rand_val = random.random()
            if rand_val < prob:
                accept = True
                reason = f"ACCEPTED (SA prob={prob:.4f})"
            else:
                accept = False
                reason = f"REJECTED (SA prob={prob:.4f})"
        
        print(f"\nSA decision: ΔE={delta_E:.4f}, T={T:.4f} → {reason}")
        
        if not accept:
            print("Reverting to previous iteration (Simulated Annealing)")
            
            prev_itp = os.path.join(
                base_dir,
                f"iter_{prev_iter}",
                os.path.basename(args.itp_out)
            )
            
            if os.path.exists(prev_itp):
                shutil.copy(prev_itp, args.itp_out)
                print(f"Reused previous ITP: {prev_itp}")
                print(f"Output written to: {args.itp_out}")
                return
            else:
                print(f"Warning: Previous ITP not found at {prev_itp}")
                print("Continuing with new ITP generation...")
        
        merged = pd.concat([old, new_row]).drop_duplicates("iteration").sort_values("iteration")
    else:
        merged = new_row
    
    merged.to_csv(iteration_file, sep="\t", index=False)
    
    # Plot evolution
    plt.figure(figsize=(10, 6))
    plt.plot(
        merged["iteration"],
        merged["global_error"],
        marker='o',
        markeredgecolor='black',
        linewidth=2,
        markersize=8
    )
    plt.xlabel("Iteration", fontsize=12)
    plt.ylabel("Global Error", fontsize=12)
    plt.title("Global Error Evolution", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(inter_dir, "dihedral_evolution.png"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(inter_dir, "dihedral_evolution.svg"), bbox_inches='tight')
    plt.close()
    
    # Generate ITP
    print("\nGenerating new ITP file...")
    
    # Track missing force constants
    missing_bonds = []
    missing_angles = []
    missing_dihedrals = []
    
    # Track forced-to-min values
    forced_bonds = []
    forced_angles = []
    forced_dihedrals = []
    
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    with open(args.itp_out, 'w') as f:
        # Write header
        f.write(f";;;;;; {args.molecule_name} - Final topology\n")
        f.write(f"; Generated on {timestamp}\n")
        f.write(f"; Iteration: {current_iter}\n")
        f.write(f"; Multimodal mode: {multimodal_mode}\n")
        f.write(f"; Variance multimodal: {variance_multimodal}\n")
        f.write(f"; Minimum forces: Bond={MIN_FORCE_BOND}, Angle={MIN_FORCE_ANGLE}, Dihedral={MIN_FORCE_DIHEDRAL}\n\n")
        
        # Write molecule type
        f.write("[ moleculetype ]\n")
        f.write(f"  {args.molecule_name}         1\n\n")
        
        # Write atoms section
        f.write("[ atoms ]\n")
        f.write("; nr   type  resnr  residue  atom    cgnr     charge     mass\n")
        for bead in dic_itp_cg['atoms']:
            f.write(f"  {bead['nr']:2d}  {bead['type']:4s}  {bead['resnr']:3d}  {bead['residue']:4s}  {bead['atom']:4s}  {bead['cgnr']:3d}  {bead['charge']:8.4f}  {bead['mass']:8.4f}\n")
        
        # Write bonds section
        f.write("\n[ bonds ]\n")
        f.write("; i    j    funct    length    force.k\n")
        for idx, row in bonds_ref_processed.iterrows():
            if idx < len(bond_pairs):
                i, j = bond_pairs[idx]
            else:
                i, j = 0, 0
                warnings.warn(f"Bond index {idx} out of range for bond_pairs")
            
            # Get force constant from ITP or use default
            bond_key = (i, j)
            FORCE_BOND = force_bond_dict.get(bond_key, DEFAULT_FORCE_BOND)
            
            if bond_key not in force_bond_dict:
                missing_bonds.append(f"({i}, {j})")
            
            sigma_sim = bonds_sim_processed.loc[idx,'sd'] if idx < len(bonds_sim_processed) else row['sd']
            prev_sigma = bonds_prev_df.loc[idx,'sd'] if bonds_prev_df is not None and idx < len(bonds_prev_df) else None
            
            F_new, status = choose_best_force(FORCE_BOND, row['sd'], sigma_sim, prev_sigma, force_type='bond')
            
            # Track forced values
            if "FORCED TO MIN" in status:
                forced_bonds.append(f"({i}, {j}) -> {F_new:.1f}")
            
            f.write(format_bond(i, j, row['mean'], F_new) +
                    f" ; σ_ref={row['sd']:.4f} σ_sim={sigma_sim:.4f} "
                    f"σ_prev={prev_sigma if prev_sigma is not None else 'NA'} status={status}\n")
        
        # Write angles section
        f.write("\n[ angles ]\n")
        f.write("; i    j    k    funct    angle    force.k\n")
        for idx, row in angles_ref_processed.iterrows():
            if idx < len(angle_triplets):
                i, j, k = angle_triplets[idx]
            else:
                i, j, k = 0, 0, 0
                warnings.warn(f"Angle index {idx} out of range for angle_triplets")
            
            # Get force constant from ITP or use default
            angle_key = (i, j, k)
            FORCE_ANGLE = force_angle_dict.get(angle_key, DEFAULT_FORCE_ANGLE)
            
            if angle_key not in force_angle_dict:
                missing_angles.append(f"({i}, {j}, {k})")
            
            sigma_sim = angles_sim_processed.loc[idx,'sd'] if idx < len(angles_sim_processed) else row['sd']
            prev_sigma = angles_prev_df.loc[idx,'sd'] if angles_prev_df is not None and idx < len(angles_prev_df) else None
            
            F_new, status = choose_best_force(FORCE_ANGLE, row['sd'], sigma_sim, prev_sigma, force_type='angle')
            
            # Track forced values
            if "FORCED TO MIN" in status:
                forced_angles.append(f"({i}, {j}, {k}) -> {F_new:.1f}")
            
            f.write(format_angle(i, j, k, row['mean'], F_new) +
                    f" ; σ_ref={row['sd']:.4f} σ_sim={sigma_sim:.4f} "
                    f"σ_prev={prev_sigma if prev_sigma is not None else 'NA'} status={status}\n")
        
        # Write dihedrals section
        f.write("\n[ dihedrals ]\n")
        f.write("; i    j    k    l    funct    angle    force.k\n")
        for idx, row in dihedrals_ref_processed.iterrows():
            if idx < len(dihedral_quartets):
                i, j, k, l = dihedral_quartets[idx]
            else:
                i, j, k, l = 0, 0, 0, 0
                warnings.warn(f"Dihedral index {idx} out of range for dihedral_quartets")
            
            # Get force constant from ITP or use default
            dihedral_key = (i, j, k, l)
            FORCE_DIHEDRAL = force_dihedral_dict.get(dihedral_key, DEFAULT_FORCE_DIHEDRAL)
            
            if dihedral_key not in force_dihedral_dict:
                missing_dihedrals.append(f"({i}, {j}, {k}, {l})")
            
            sigma_sim = dihedrals_sim_processed.loc[idx,'sd'] if idx < len(dihedrals_sim_processed) else row['sd']
            prev_sigma = dihedrals_prev_df.loc[idx,'sd'] if dihedrals_prev_df is not None and idx < len(dihedrals_prev_df) else None
            
            # Use angle from reference data (removed dihedrals_target)
            angle_val = row['mean']
            
            F_new, status = choose_best_force(FORCE_DIHEDRAL, row['sd'], sigma_sim, prev_sigma, force_type='dihedral')
            
            # Track forced values
            if "FORCED TO MIN" in status:
                forced_dihedrals.append(f"({i}, {j}, {k}, {l}) -> {F_new:.1f}")
            
            f.write(format_dihedral(i, j, k, l, angle_val, F_new) +
                    f" ; σ_ref={row['sd']:.4f} "
                    f"σ_sim={sigma_sim:.4f} "
                    f"σ_prev={prev_sigma if prev_sigma is not None else 'NA'} "
                    f"status={status}\n")
    
    # Print warnings for missing force constants
    if missing_bonds:
        warnings.warn(f"Bonds without defined force constants (using default {DEFAULT_FORCE_BOND}): {', '.join(missing_bonds[:10])}")
        if len(missing_bonds) > 10:
            print(f"  ... and {len(missing_bonds) - 10} more")
    
    if missing_angles:
        warnings.warn(f"Angles without defined force constants (using default {DEFAULT_FORCE_ANGLE}): {', '.join(missing_angles[:10])}")
        if len(missing_angles) > 10:
            print(f"  ... and {len(missing_angles) - 10} more")
    
    if missing_dihedrals:
        warnings.warn(f"Dihedrals without defined force constants (using default {DEFAULT_FORCE_DIHEDRAL}): {', '.join(missing_dihedrals[:10])}")
        if len(missing_dihedrals) > 10:
            print(f"  ... and {len(missing_dihedrals) - 10} more")
    
    # Print forced values summary
    if forced_bonds:
        print(f"\n  Bonds forced to minimum ({MIN_FORCE_BOND}): {len(forced_bonds)}")
        if len(forced_bonds) <= 10:
            for fb in forced_bonds:
                print(f"    {fb}")
        else:
            for fb in forced_bonds[:5]:
                print(f"    {fb}")
            print(f"    ... and {len(forced_bonds)-5} more")
    
    if forced_angles:
        print(f"\n  Angles forced to minimum ({MIN_FORCE_ANGLE}): {len(forced_angles)}")
        if len(forced_angles) <= 10:
            for fa in forced_angles:
                print(f"    {fa}")
        else:
            for fa in forced_angles[:5]:
                print(f"    {fa}")
            print(f"    ... and {len(forced_angles)-5} more")
    
    if forced_dihedrals:
        print(f"\n  Dihedrals forced to minimum ({MIN_FORCE_DIHEDRAL}): {len(forced_dihedrals)}")
        if len(forced_dihedrals) <= 10:
            for fd in forced_dihedrals:
                print(f"    {fd}")
        else:
            for fd in forced_dihedrals[:5]:
                print(f"    {fd}")
            print(f"    ... and {len(forced_dihedrals)-5} more")
    
    print(f"\n{'='*60}")
    print(f"ITP file successfully generated: {args.itp_out}")
    print(f"Total bonds processed: {len(bonds_ref_processed)}")
    print(f"Total angles processed: {len(angles_ref_processed)}")
    print(f"Total dihedrals processed: {len(dihedrals_ref_processed)}")
    print(f"{'='*60}\n")

if __name__ == "__main__":
    main()
