#!/usr/bin/env python3
"""
Description:

Automated optimization script for Martini coarse-grained force field parameters (.itp files),
using reference and simulation data derived from XVG time series.

The script iteratively adjusts bond, angle, and dihedral force constants to improve statistical
agreement between reference and simulated distributions, while preserving the original ITP
structure (only specific interaction sections are modified).

Core Features:

- XVG data processing:
  - Reads time-series data from XVG files for bonds, angles, and dihedrals.
  - Builds normalized probability distributions using KDE or histogram fallback.

- Statistical comparison:
  - Computes R² (coefficient of determination) between reference and simulation distributions.
  - Extracts mean and standard deviation for each interaction type.

- Force constant optimization:
  - Bayesian update of force constants based on uncertainty in distributions.
  - R²-guided correction factor to improve agreement with reference data.
  - Enforces physical bounds (min/max force constants per interaction type).
  - Stochastic sampling to explore multiple candidate solutions.

- Iterative optimization control:
  - Simulated Annealing acceptance criterion based on global error improvement (1 - mean R²).
  - Temperature decay across iterations to reduce randomness over time.

- Structural preservation:
  - Preserves exact ordering and formatting of the original ITP file.
  - Only replaces bonds, angles, and dihedrals sections.
  - Keeps all other topology information unchanged.

- Analysis and output:
  - Tracks evolution of R² and error per iteration.
  - Generates plots for global, bond, angle, and dihedral convergence.
  - Stores iteration history in TSV format.
  - Outputs updated optimized ITP file.

------------------------------------------------------------
Usage Example:

1. Basic execution:

python 7-potential_adjustment.py \
    --bonds_ref_xvg_dir ref/bonds \
    --angles_ref_xvg_dir ref/angles \
    --dihedrals_ref_xvg_dir ref/dihedrals \
    --bonds_sim_xvg_dir sim/iter_1/bonds \
    --angles_sim_xvg_dir sim/iter_1/angles \
    --dihedrals_sim_xvg_dir sim/iter_1/dihedrals \
    --itp_cg input/cg.itp \
    --ndx_bounds NDX/bonds.ndx \
    --ndx_angles NDX/angles.ndx \
    --ndx_dihedrals NDX/dihedrals.ndx \
    --itp_out output/cg_optimized.itp

------------------------------------------------------------

2. Example with custom parameters:

python 7-potential_adjustment.py \
    --bonds_ref_xvg_dir ref/bonds \
    --angles_ref_xvg_dir ref/angles \
    --dihedrals_ref_xvg_dir ref/dihedrals \
    --bonds_sim_xvg_dir sim/iter_5/bonds \
    --angles_sim_xvg_dir sim/iter_5/angles \
    --dihedrals_sim_xvg_dir sim/iter_5/dihedrals \
    --itp_cg input/cg.itp \
    --ndx_bounds NDX/bonds.ndx \
    --ndx_angles NDX/angles.ndx \
    --ndx_dihedrals NDX/dihedrals.ndx \
    --itp_out output/cg_optimized.itp \
    --distribution_points 100 \
    --T0 15.0 \
    --alpha 0.9

------------------------------------------------------------

3. Typical workflow:

- Run molecular simulation → generate XVG files
- Compare reference vs simulation distributions
- Update force field parameters iteratively
- Monitor convergence via R² evolution plots
- Repeat until convergence is achieved

------------------------------------------------------------

Output:
- Optimized ITP file (cg.itp)
- Iteration log (force_iteration.tsv)
- Convergence plots:
  - global_r2.png
  - bonds_r2.png
  - angles_r2.png
  - dihedrals_r2.png
"""
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import shutil
import math
import random
import re
from datetime import datetime
import warnings
from scipy import stats

# ============================================
# Constants
# ============================================

MIN_FORCE_BOND = 750.0
MIN_FORCE_ANGLE = 15.0
MIN_FORCE_DIHEDRAL = 15.0
MAX_FORCE_BOND = 10000.0
MAX_FORCE_ANGLE = 150.0
MAX_FORCE_DIHEDRAL = 150.0
DEFAULT_FORCE_BOND = 1250.0
DEFAULT_FORCE_ANGLE = 25.0
DEFAULT_FORCE_DIHEDRAL = 25.0
DEFAULT_DISTRIBUTION_POINTS = 50
DEFAULT_PREFIX_BOND = "bond_"
DEFAULT_PREFIX_ANGLE = "ang_"
DEFAULT_PREFIX_DIHEDRAL = "dih_"

# ============================================
# Helper functions
# ============================================

def read_ndx(file_path):
    """Read NDX file and return list of index tuples."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} not found.")
    
    indices = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('[') or line.startswith(';'):
                continue
            
            parts = line.split()
            if not parts:
                continue
            
            numeric_parts = []
            for x in parts:
                try:
                    num = int(x)
                    numeric_parts.append(num)
                except ValueError:
                    if x.startswith(';'):
                        break
                    continue
            
            if len(numeric_parts) >= 2:
                indices.append(tuple(numeric_parts))
    
    return indices

def format_bond(i, j, length, force):
    return f"  {i:2d}  {j:2d}   1    {length:.3f}    {force:.1f}   ;"

def format_angle(i, j, k, angle, force):
    return f"  {i:2d}  {j:2d}  {k:2d}   2    {angle:.1f}    {force:.1f}   ;"

def format_dihedral(i, j, k, l, angle, force):
    return f"  {i:2d}  {j:2d}  {k:2d}  {l:2d}   2    {angle:.1f}    {force:.1f}   ;"

def parse_xvg(xvg_file):
    """Parse XVG file and extract time series data."""
    if not os.path.exists(xvg_file):
        raise FileNotFoundError(f"XVG file {xvg_file} not found.")
    
    times = []
    values = []
    
    with open(xvg_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or line.startswith('@'):
                continue
            if not line:
                continue
            
            parts = line.split()
            if len(parts) >= 2:
                try:
                    time = float(parts[0])
                    value = float(parts[1])
                    times.append(time)
                    values.append(value)
                except ValueError:
                    continue
    
    return np.array(times), np.array(values)

def load_xvg_series(xvg_dir, prefix, num_interactions):
    """Load all XVG files for a specific interaction type."""
    if not os.path.exists(xvg_dir):
        raise FileNotFoundError(f"XVG directory {xvg_dir} not found.")
    
    series_dict = {}
    
    for idx in range(num_interactions):
        xvg_file = os.path.join(xvg_dir, f"{prefix}{idx}.xvg")
        
        if not os.path.exists(xvg_file):
            warnings.warn(f"XVG file not found: {xvg_file}")
            continue
        
        try:
            times, values = parse_xvg(xvg_file)
            series_dict[idx] = (times, values)
        except Exception as e:
            warnings.warn(f"Error parsing {xvg_file}: {e}")
            continue
    
    return series_dict

def compute_probability_density(data, n_points=DEFAULT_DISTRIBUTION_POINTS, 
                                x_min=None, x_max=None):
    """Compute normalized probability density function from data using KDE."""
    if len(data) == 0:
        return np.array([]), np.array([])
    
    if x_min is None:
        x_min = np.min(data)
    if x_max is None:
        x_max = np.max(data)
    
    padding = 0.05 * (x_max - x_min) if x_max > x_min else 0.1
    x_min_padded = x_min - padding
    x_max_padded = x_max + padding
    
    x_grid = np.linspace(x_min_padded, x_max_padded, n_points)
    
    try:
        kde = stats.gaussian_kde(data)
        pdf = kde.evaluate(x_grid)
        
        try:
            integral = np.trapezoid(pdf, x_grid)
        except AttributeError:
            integral = np.trapz(pdf, x_grid)
        
        if integral > 0:
            pdf_normalized = pdf / integral
        else:
            pdf_normalized = pdf
            
        return x_grid, pdf_normalized
        
    except Exception as e:
        warnings.warn(f"KDE failed: {e}. Using histogram instead.")
        hist, bin_edges = np.histogram(data, bins=min(n_points//2, 20), density=True)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        return bin_centers, hist

def calculate_r2_from_xvg(ref_series, sim_series, n_points=DEFAULT_DISTRIBUTION_POINTS):
    """Calculate R² between reference and simulated distributions using XVG data."""
    ref_times, ref_values = ref_series
    sim_times, sim_values = sim_series
    
    if len(ref_values) == 0 or len(sim_values) == 0:
        return 0.0, np.nan, np.nan, np.nan, np.nan
    
    ref_mean = np.mean(ref_values)
    ref_std = np.std(ref_values)
    sim_mean = np.mean(sim_values)
    sim_std = np.std(sim_values)
    
    x_min = min(np.min(ref_values), np.min(sim_values))
    x_max = max(np.max(ref_values), np.max(sim_values))
    
    try:
        x_grid, ref_pdf = compute_probability_density(ref_values, n_points, x_min, x_max)
        _, sim_pdf = compute_probability_density(sim_values, n_points, x_min, x_max)
        
        if len(x_grid) != len(ref_pdf) or len(x_grid) != len(sim_pdf):
            return 0.0, ref_mean, ref_std, sim_mean, sim_std
        
        ss_res = np.sum((ref_pdf - sim_pdf)**2)
        ss_tot = np.sum((ref_pdf - np.mean(ref_pdf))**2)
        
        r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
        r2 = max(0.0, min(1.0, r2))
        
        return r2, ref_mean, ref_std, sim_mean, sim_std
        
    except Exception as e:
        warnings.warn(f"R² calculation failed: {e}. Using fallback method.")
        try:
            hist_ref, bins = np.histogram(ref_values, bins=min(n_points//2, 20), density=True)
            hist_sim, _ = np.histogram(sim_values, bins=bins, density=True)
            
            min_len = min(len(hist_ref), len(hist_sim))
            hist_ref = hist_ref[:min_len]
            hist_sim = hist_sim[:min_len]
            
            ss_res = np.sum((hist_ref - hist_sim)**2)
            ss_tot = np.sum((hist_ref - np.mean(hist_ref))**2)
            r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
            r2 = max(0.0, min(1.0, r2))
            
            return r2, ref_mean, ref_std, sim_mean, sim_std
        except Exception as e2:
            warnings.warn(f"Fallback method also failed: {e2}. Returning 0.0")
            return 0.0, ref_mean, ref_std, sim_mean, sim_std

def parse_itp_cg_preserve_order(itp_file):
    """
    Parse ITP file preserving EXACT order of all sections and lines.
    Returns a list of sections in order, each section is (section_name, list_of_lines).
    """
    if not os.path.exists(itp_file):
        raise FileNotFoundError(f"ITP file {itp_file} not found.")
    
    sections = []  # List of (section_name, lines_list)
    current_section = None
    current_lines = []
    
    with open(itp_file, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].rstrip('\n')
        line_stripped = line.strip()
        
        # Check for section header
        if line_stripped.startswith('[') and line_stripped.endswith(']'):
            # Save previous section
            if current_section is not None:
                sections.append((current_section, current_lines))
            
            # Start new section
            current_section = line_stripped[1:-1].strip()
            current_lines = [line]
            i += 1
            continue
        
        # Add line to current section or as header
        if current_section is not None:
            current_lines.append(line)
        else:
            # Lines before first section (header/comments)
            if len(sections) == 0 or sections[0][0] != '__header__':
                sections.insert(0, ('__header__', []))
            sections[0][1].append(line)
        
        i += 1
    
    # Save last section
    if current_section is not None:
        sections.append((current_section, current_lines))
    
    return sections

def parse_bonds_from_section(lines):
    """Parse bond information from a bonds section."""
    bonds = []
    for line in lines:
        line_stripped = line.strip()
        if not line_stripped or line_stripped.startswith(';') or line_stripped.startswith('['):
            continue
        parts = line.split()
        if len(parts) >= 6:
            try:
                bond = {
                    'i': int(parts[0]),
                    'j': int(parts[1]),
                    'funct': int(parts[2]),
                    'length': float(parts[3]),
                    'force_k': float(parts[4]),
                    'line': line
                }
                bonds.append(bond)
            except ValueError:
                continue
    return bonds

def parse_angles_from_section(lines):
    """Parse angle information from an angles section."""
    angles = []
    for line in lines:
        line_stripped = line.strip()
        if not line_stripped or line_stripped.startswith(';') or line_stripped.startswith('['):
            continue
        parts = line.split()
        if len(parts) >= 7:
            try:
                angle = {
                    'i': int(parts[0]),
                    'j': int(parts[1]),
                    'k': int(parts[2]),
                    'funct': int(parts[3]),
                    'angle': float(parts[4]),
                    'force_k': float(parts[5]),
                    'line': line
                }
                angles.append(angle)
            except ValueError:
                continue
    return angles

def parse_dihedrals_from_section(lines):
    """Parse dihedral information from a dihedrals section."""
    dihedrals = []
    for line in lines:
        line_stripped = line.strip()
        if not line_stripped or line_stripped.startswith(';') or line_stripped.startswith('['):
            continue
        parts = line.split()
        if len(parts) >= 8:
            try:
                dihedral = {
                    'i': int(parts[0]),
                    'j': int(parts[1]),
                    'k': int(parts[2]),
                    'l': int(parts[3]),
                    'funct': int(parts[4]),
                    'angle': float(parts[5]),
                    'force_k': float(parts[6]),
                    'line': line
                }
                dihedrals.append(dihedral)
            except ValueError:
                continue
    return dihedrals

def choose_best_force(F_orig, sigma_ref, sigma_curr, sigma_prev, r2, force_type='bond'):
    """Choose best force constant with Bayesian update and R² correction."""
    if force_type == 'bond':
        min_force = MIN_FORCE_BOND
        max_force = MAX_FORCE_BOND
    elif force_type == 'angle':
        min_force = MIN_FORCE_ANGLE
        max_force = MAX_FORCE_ANGLE
    else:
        min_force = MIN_FORCE_DIHEDRAL
        max_force = MAX_FORCE_DIHEDRAL
    
    if sigma_curr == 0:
        return max(min_force, min(max_force, F_orig)), "NOT CHANGED (zero sigma)"
    
    r2_correction = 1.0 + (1.0 - r2)
    
    if sigma_prev is not None and not np.isnan(sigma_prev):
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
        F_cand = F_orig * (s / sigma_ref) * r2_correction
        F_cand = max(min_force, min(max_force, F_cand))
        F_candidates.append(F_cand)
    
    if not F_candidates:
        return max(min_force, min(max_force, F_orig)), "NOT CHANGED (no samples)"
    
    best_F = max(min_force, min(max_force, F_orig))
    best_error = abs(sigma_ref - sigma_curr)
    
    for F_cand in F_candidates:
        sigma_pred = sigma_curr * (F_orig / F_cand) if F_cand > 0 else sigma_curr
        err = abs(sigma_ref - sigma_pred)
        if err < best_error:
            best_error = err
            best_F = F_cand
    
    if best_F < min_force:
        best_F = min_force
        return best_F, f"FORCED TO MIN ({min_force:.1f})"
    elif best_F > max_force:
        best_F = max_force
        return best_F, f"FORCED TO MAX ({max_force:.1f})"
    
    if best_F == F_orig:
        return best_F, "NOT CHANGED (Bayes reject)"
    else:
        return best_F, f"CHANGED (Bayes, R²={r2:.3f})"

def compute_global_error_from_r2(r2_values):
    if not r2_values:
        return 1.0
    return 1.0 - np.mean(r2_values)

def plot_evolution(impact_dir, iteration_data, title_prefix, ylabel_error, ylabel_r2, color_error='red', color_r2='blue'):
    """
    Generic function to plot evolution with integer x-axis spacing.
    """
    # Determine column names based on title_prefix
    error_col = f"{title_prefix.lower()}_error"
    
    # Special case for 'global' - use 'global_r2' instead of 'global_r2' (same name)
    r2_col = f"{title_prefix.lower()}_r2"
    
    # Check if columns exist
    if error_col not in iteration_data.columns:
        print(f"Warning: Column '{error_col}' not found in iteration data")
        print(f"Available columns: {list(iteration_data.columns)}")
        return
    
    if r2_col not in iteration_data.columns:
        print(f"Warning: Column '{r2_col}' not found in iteration data")
        print(f"Available columns: {list(iteration_data.columns)}")
        return
    
    # Plot error evolution
    plt.figure(figsize=(10, 5))
    plt.plot(
        iteration_data["iteration"],
        iteration_data[error_col],
        marker='o',
        markeredgecolor='black',
        linewidth=2,
        markersize=8,
        color=color_error
    )
    plt.xlabel("Iteration", fontsize=12)
    plt.ylabel(ylabel_error, fontsize=12)
    plt.title(f"{title_prefix} Error Evolution", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    plt.tight_layout()
    plt.savefig(os.path.join(impact_dir, f"{title_prefix.lower()}_error.png"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(impact_dir, f"{title_prefix.lower()}_error.svg"), bbox_inches='tight')
    plt.close()
    
    # Plot R² evolution
    plt.figure(figsize=(10, 5))
    plt.plot(
        iteration_data["iteration"],
        iteration_data[r2_col],
        marker='s',
        markeredgecolor='black',
        linewidth=2,
        markersize=8,
        color=color_r2
    )
    plt.xlabel("Iteration", fontsize=12)
    plt.ylabel(ylabel_r2, fontsize=12)
    plt.title(f"{title_prefix} R² Evolution", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    plt.tight_layout()
    plt.savefig(os.path.join(impact_dir, f"{title_prefix.lower()}_r2.png"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(impact_dir, f"{title_prefix.lower()}_r2.svg"), bbox_inches='tight')
    plt.close()

# ============================================
# Main function
# ============================================

def main():
    parser = argparse.ArgumentParser(description="Generate Martini-style .itp with automatic force constant adjustment")
    parser.add_argument("--bonds_ref_xvg_dir", required=True)
    parser.add_argument("--angles_ref_xvg_dir", required=True)
    parser.add_argument("--dihedrals_ref_xvg_dir", required=True)
    parser.add_argument("--bonds_sim_xvg_dir", required=True)
    parser.add_argument("--angles_sim_xvg_dir", required=True)
    parser.add_argument("--dihedrals_sim_xvg_dir", required=True)
    parser.add_argument("--itp_cg", required=True)
    parser.add_argument("--ndx_bounds", required=True)
    parser.add_argument("--ndx_angles", required=True)
    parser.add_argument("--ndx_dihedrals", required=True)
    parser.add_argument("--molecule_name", default="molecule")
    parser.add_argument("--distribution_points", type=int, default=DEFAULT_DISTRIBUTION_POINTS)
    parser.add_argument("--T0", type=float, default=10.0)
    parser.add_argument("--alpha", type=float, default=0.85)
    parser.add_argument("--itp_out", default="cg.itp")
    parser.add_argument("--prefix_xvg_bond_ref", default=DEFAULT_PREFIX_BOND)
    parser.add_argument("--prefix_xvg_angle_ref", default=DEFAULT_PREFIX_ANGLE)
    parser.add_argument("--prefix_xvg_dihedral_ref", default=DEFAULT_PREFIX_DIHEDRAL)
    parser.add_argument("--prefix_xvg_bond_sim", default=DEFAULT_PREFIX_BOND)
    parser.add_argument("--prefix_xvg_angle_sim", default=DEFAULT_PREFIX_ANGLE)
    parser.add_argument("--prefix_xvg_dihedral_sim", default=DEFAULT_PREFIX_DIHEDRAL)
    
    args = parser.parse_args()
    
    # Detect iteration
    iter_match = re.search(r'iter_(\d+)', args.bonds_sim_xvg_dir)
    current_iter = int(iter_match.group(1)) if iter_match else 1
    prev_iter = current_iter - 1
    
    base_dir = os.path.abspath(os.path.join(args.bonds_sim_xvg_dir, "../../../"))
    prev_sim_dir = os.path.join(base_dir, f"iter_{prev_iter}", "sim", "XVG") if prev_iter >= 0 else None
    
    print(f"\n{'='*60}")
    print(f"Starting iteration {current_iter}")
    print(f"{'='*60}")
    print(f"Distribution points: {args.distribution_points}")
    print(f"T0: {args.T0}, alpha: {args.alpha}")
    print(f"Minimum forces: Bond={MIN_FORCE_BOND}, Angle={MIN_FORCE_ANGLE}, Dihedral={MIN_FORCE_DIHEDRAL}")
    print(f"Maximum forces: Bond={MAX_FORCE_BOND}, Angle={MAX_FORCE_ANGLE}, Dihedral={MAX_FORCE_DIHEDRAL}")
    print(f"{'='*60}\n")
    
    # Load topology indices
    print("Loading topology indices...")
    bond_pairs = read_ndx(args.ndx_bounds)
    angle_triplets = read_ndx(args.ndx_angles)
    dihedral_quartets = read_ndx(args.ndx_dihedrals)
    
    num_bonds = len(bond_pairs)
    num_angles = len(angle_triplets)
    num_dihedrals = len(dihedral_quartets)
    
    print(f"Loaded {num_bonds} bonds, {num_angles} angles, {num_dihedrals} dihedrals from NDX files")
    
    # Load reference XVG data
    print("\nLoading reference XVG data...")
    ref_bonds_dict = load_xvg_series(args.bonds_ref_xvg_dir, args.prefix_xvg_bond_ref, num_bonds)
    ref_angles_dict = load_xvg_series(args.angles_ref_xvg_dir, args.prefix_xvg_angle_ref, num_angles)
    ref_dihedrals_dict = load_xvg_series(args.dihedrals_ref_xvg_dir, args.prefix_xvg_dihedral_ref, num_dihedrals)
    
    print(f"  Reference bonds loaded: {len(ref_bonds_dict)}/{num_bonds}")
    print(f"  Reference angles loaded: {len(ref_angles_dict)}/{num_angles}")
    print(f"  Reference dihedrals loaded: {len(ref_dihedrals_dict)}/{num_dihedrals}")
    
    # Load simulated XVG data
    print("\nLoading simulated XVG data...")
    sim_bonds_dict = load_xvg_series(args.bonds_sim_xvg_dir, args.prefix_xvg_bond_sim, num_bonds)
    sim_angles_dict = load_xvg_series(args.angles_sim_xvg_dir, args.prefix_xvg_angle_sim, num_angles)
    sim_dihedrals_dict = load_xvg_series(args.dihedrals_sim_xvg_dir, args.prefix_xvg_dihedral_sim, num_dihedrals)
    
    print(f"  Simulated bonds loaded: {len(sim_bonds_dict)}/{num_bonds}")
    print(f"  Simulated angles loaded: {len(sim_angles_dict)}/{num_angles}")
    print(f"  Simulated dihedrals loaded: {len(sim_dihedrals_dict)}/{num_dihedrals}")
    
    # Load previous data
    prev_bonds_dict = {}
    prev_angles_dict = {}
    prev_dihedrals_dict = {}
    if prev_sim_dir and os.path.exists(prev_sim_dir):
        print(f"\nLoading previous iteration (iter_{prev_iter}) simulated data...")
        prev_bonds_dict = load_xvg_series(os.path.join(prev_sim_dir, "bonds"), args.prefix_xvg_bond_sim, num_bonds)
        prev_angles_dict = load_xvg_series(os.path.join(prev_sim_dir, "angles"), args.prefix_xvg_angle_sim, num_angles)
        prev_dihedrals_dict = load_xvg_series(os.path.join(prev_sim_dir, "dihedrals"), args.prefix_xvg_dihedral_sim, num_dihedrals)
    
    # Calculate R² and statistics
    print("\nCalculating R² values and statistics...")
    
    bonds_data = []
    for idx in range(num_bonds):
        if idx in ref_bonds_dict and idx in sim_bonds_dict:
            r2, ref_mean, ref_std, sim_mean, sim_std = calculate_r2_from_xvg(
                ref_bonds_dict[idx], sim_bonds_dict[idx], args.distribution_points)
            prev_std = None
            if idx in prev_bonds_dict:
                _, _, _, _, prev_std = calculate_r2_from_xvg(
                    ref_bonds_dict[idx], prev_bonds_dict[idx], args.distribution_points)
            bonds_data.append((r2, ref_mean, ref_std, sim_mean, sim_std, prev_std))
        else:
            bonds_data.append((0.0, np.nan, np.nan, np.nan, np.nan, None))
    
    angles_data = []
    for idx in range(num_angles):
        if idx in ref_angles_dict and idx in sim_angles_dict:
            r2, ref_mean, ref_std, sim_mean, sim_std = calculate_r2_from_xvg(
                ref_angles_dict[idx], sim_angles_dict[idx], args.distribution_points)
            prev_std = None
            if idx in prev_angles_dict:
                _, _, _, _, prev_std = calculate_r2_from_xvg(
                    ref_angles_dict[idx], prev_angles_dict[idx], args.distribution_points)
            angles_data.append((r2, ref_mean, ref_std, sim_mean, sim_std, prev_std))
        else:
            angles_data.append((0.0, np.nan, np.nan, np.nan, np.nan, None))
    
    dihedrals_data = []
    for idx in range(num_dihedrals):
        if idx in ref_dihedrals_dict and idx in sim_dihedrals_dict:
            r2, ref_mean, ref_std, sim_mean, sim_std = calculate_r2_from_xvg(
                ref_dihedrals_dict[idx], sim_dihedrals_dict[idx], args.distribution_points)
            prev_std = None
            if idx in prev_dihedrals_dict:
                _, _, _, _, prev_std = calculate_r2_from_xvg(
                    ref_dihedrals_dict[idx], prev_dihedrals_dict[idx], args.distribution_points)
            dihedrals_data.append((r2, ref_mean, ref_std, sim_mean, sim_std, prev_std))
        else:
            dihedrals_data.append((0.0, np.nan, np.nan, np.nan, np.nan, None))
    
    bonds_r2 = [d[0] for d in bonds_data]
    bonds_ref_mean = [d[1] for d in bonds_data]
    bonds_ref_std = [d[2] for d in bonds_data]
    bonds_sim_std = [d[4] for d in bonds_data]
    bonds_prev_std = [d[5] for d in bonds_data]
    
    angles_r2 = [d[0] for d in angles_data]
    angles_ref_mean = [d[1] for d in angles_data]
    angles_ref_std = [d[2] for d in angles_data]
    angles_sim_std = [d[4] for d in angles_data]
    angles_prev_std = [d[5] for d in angles_data]
    
    dihedrals_r2 = [d[0] for d in dihedrals_data]
    dihedrals_ref_mean = [d[1] for d in dihedrals_data]
    dihedrals_ref_std = [d[2] for d in dihedrals_data]
    dihedrals_sim_std = [d[4] for d in dihedrals_data]
    dihedrals_prev_std = [d[5] for d in dihedrals_data]
    
    all_r2 = bonds_r2 + angles_r2 + dihedrals_r2
    global_error = compute_global_error_from_r2(all_r2)
    
    print(f"\nGlobal error (1 - mean R²): {global_error:.6f}")
    print(f"Mean R²: {1.0 - global_error:.6f}")
    print(f"  Bonds mean R²: {np.mean(bonds_r2) if bonds_r2 else 0.0:.6f}")
    print(f"  Angles mean R²: {np.mean(angles_r2) if angles_r2 else 0.0:.6f}")
    print(f"  Dihedrals mean R²: {np.mean(dihedrals_r2) if dihedrals_r2 else 0.0:.6f}")
    
    # Create iteration tracking
    impact_dir = os.path.abspath(os.path.join(os.path.dirname(args.itp_out), "../impact_of_potentials"))
    os.makedirs(impact_dir, exist_ok=True)
    
    iteration_file = os.path.join(impact_dir, "force_iteration.tsv")
    
    bonds_error = compute_global_error_from_r2(bonds_r2) if bonds_r2 else 1.0
    angles_error = compute_global_error_from_r2(angles_r2) if angles_r2 else 1.0
    dihedrals_error = compute_global_error_from_r2(dihedrals_r2) if dihedrals_r2 else 1.0
    bonds_mean_r2 = np.mean(bonds_r2) if bonds_r2 else 0.0
    angles_mean_r2 = np.mean(angles_r2) if angles_r2 else 0.0
    dihedrals_mean_r2 = np.mean(dihedrals_r2) if dihedrals_r2 else 0.0
    
    new_row = pd.DataFrame({
        "iteration": [current_iter],
        "global_error": [global_error],
        "global_r2": [1.0 - global_error],
        "bonds_error": [bonds_error],
        "bonds_r2": [bonds_mean_r2],
        "angles_error": [angles_error],
        "angles_r2": [angles_mean_r2],
        "dihedrals_error": [dihedrals_error],
        "dihedrals_r2": [dihedrals_mean_r2]
    })
    
    T = args.T0 * (args.alpha ** current_iter)
    
    if os.path.exists(iteration_file):
        old = pd.read_csv(iteration_file, sep="\t")
        prev_last = old.sort_values("iteration").iloc[-1]["global_error"]
        delta_E = global_error - prev_last
        
        if delta_E <= 0:
            accept = True
            reason = "IMPROVED (R² increased)"
        else:
            prob = math.exp(-delta_E / T) if T > 0 else 0
            rand_val = random.random()
            if rand_val < prob:
                accept = True
                reason = f"ACCEPTED (SA prob={prob:.4f})"
            else:
                accept = False
                reason = f"REJECTED (SA prob={prob:.4f})"
        
        print(f"\nSA decision: ΔError={delta_E:.4f}, T={T:.4f} → {reason}")
        
        if not accept:
            print("Reverting to previous iteration (Simulated Annealing)")
            prev_itp = os.path.join(base_dir, f"iter_{prev_iter}", os.path.basename(args.itp_out))
            if os.path.exists(prev_itp):
                shutil.copy(prev_itp, args.itp_out)
                print(f"Reused previous ITP: {prev_itp}")
                print(f"Output written to: {args.itp_out}")
                return
        
        merged = pd.concat([old, new_row]).drop_duplicates("iteration").sort_values("iteration")
    else:
        merged = new_row
    
    merged.to_csv(iteration_file, sep="\t", index=False)
    
    print("\nGenerating evolution plots...")
    plot_evolution(impact_dir, merged, "global", "Global Error (1 - R²)", "Global R²", 'red', 'blue')
    plot_evolution(impact_dir, merged, "bonds", "Bonds Error (1 - R²)", "Bonds Mean R²", 'red', 'blue')
    plot_evolution(impact_dir, merged, "angles", "Angles Error (1 - R²)", "Angles Mean R²", 'red', 'blue')
    plot_evolution(impact_dir, merged, "dihedrals", "Dihedrals Error (1 - R²)", "Dihedrals Mean R²", 'red', 'blue')
    
    # Parse original ITP preserving exact order
    print(f"\nParsing ITP CG file: {args.itp_cg}")
    sections = parse_itp_cg_preserve_order(args.itp_cg)
    
    # Extract force constants from original sections
    force_bond_dict = {}
    force_angle_dict = {}
    force_dihedral_dict = {}
    
    for section_name, section_lines in sections:
        if section_name == 'bonds':
            bonds = parse_bonds_from_section(section_lines)
            for bond in bonds:
                key = (bond['i'], bond['j'])
                force_bond_dict[key] = bond['force_k']
        elif section_name == 'angles':
            angles = parse_angles_from_section(section_lines)
            for angle in angles:
                key = (angle['i'], angle['j'], angle['k'])
                force_angle_dict[key] = angle['force_k']
        elif section_name == 'dihedrals':
            dihedrals = parse_dihedrals_from_section(section_lines)
            for dihedral in dihedrals:
                key = (dihedral['i'], dihedral['j'], dihedral['k'], dihedral['l'])
                force_dihedral_dict[key] = dihedral['force_k']
    
    print(f"Loaded {len(force_bond_dict)} bonds, {len(force_angle_dict)} angles, "
          f"{len(force_dihedral_dict)} dihedrals from ITP")
    
    # Track missing force constants
    missing_bonds = []
    missing_angles = []
    missing_dihedrals = []
    
    # Track forced-to-min/max values
    forced_bonds = []
    forced_angles = []
    forced_dihedrals = []
    
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Generate new ITP preserving exact section order
    with open(args.itp_out, 'w') as f:
        # Write header
        f.write(f";;;;;; {args.molecule_name} - Final topology\n")
        f.write(f"; Generated on {timestamp}\n")
        f.write(f"; Iteration: {current_iter}\n")
        f.write(f"; Distribution points: {args.distribution_points}\n")
        f.write(f"; Minimum forces: Bond={MIN_FORCE_BOND}, Angle={MIN_FORCE_ANGLE}, Dihedral={MIN_FORCE_DIHEDRAL}\n")
        f.write(f"; Maximum forces: Bond={MAX_FORCE_BOND}, Angle={MAX_FORCE_ANGLE}, Dihedral={MAX_FORCE_DIHEDRAL}\n\n")
        
        # Write sections in original order, replacing bonds/angles/dihedrals
        for section_name, section_lines in sections:
            if section_name == '__header__':
                # Write header lines (skip empty ones)
                for line in section_lines:
                    if line.strip():
                        f.write(f"{line}\n")
                continue
            
            if section_name == 'bonds':
                # Write new bonds section
                f.write("[ bonds ]\n")
                f.write("; i    j    funct    length    force.k\n")
                for idx in range(num_bonds):
                    if idx < len(bond_pairs):
                        i, j = bond_pairs[idx]
                    else:
                        i, j = 0, 0
                        warnings.warn(f"Bond index {idx} out of range for bond_pairs")
                    
                    bond_key = (i, j)
                    FORCE_BOND = force_bond_dict.get(bond_key, DEFAULT_FORCE_BOND)
                    
                    if bond_key not in force_bond_dict:
                        missing_bonds.append(f"({i}, {j})")
                    
                    ref_mean = bonds_ref_mean[idx] if idx < len(bonds_ref_mean) and not np.isnan(bonds_ref_mean[idx]) else 0.0
                    ref_std = bonds_ref_std[idx] if idx < len(bonds_ref_std) and not np.isnan(bonds_ref_std[idx]) else 0.0
                    sim_std = bonds_sim_std[idx] if idx < len(bonds_sim_std) and not np.isnan(bonds_sim_std[idx]) else 0.0
                    prev_std = bonds_prev_std[idx] if idx < len(bonds_prev_std) and bonds_prev_std[idx] is not None and not np.isnan(bonds_prev_std[idx]) else None
                    r2 = bonds_r2[idx] if idx < len(bonds_r2) else 0.0
                    
                    F_new, status = choose_best_force(FORCE_BOND, ref_std, sim_std, prev_std, r2, force_type='bond')
                    
                    if "FORCED TO MIN" in status:
                        forced_bonds.append(f"({i}, {j}) -> {F_new:.1f}")
                    elif "FORCED TO MAX" in status:
                        forced_bonds.append(f"({i}, {j}) -> {F_new:.1f}")
                    
                    f.write(format_bond(i, j, ref_mean, F_new) +
                            f" ; σ_ref={ref_std:.4f} σ_sim={sim_std:.4f} "
                            f"σ_prev={prev_std if prev_std is not None else 'NA'} "
                            f"R²={r2:.4f} status={status}\n")
                f.write("\n")
                
            elif section_name == 'angles':
                # Write new angles section
                f.write("[ angles ]\n")
                f.write("; i    j    k    funct    angle    force.k\n")
                for idx in range(num_angles):
                    if idx < len(angle_triplets):
                        i, j, k = angle_triplets[idx]
                    else:
                        i, j, k = 0, 0, 0
                        warnings.warn(f"Angle index {idx} out of range for angle_triplets")
                    
                    angle_key = (i, j, k)
                    FORCE_ANGLE = force_angle_dict.get(angle_key, DEFAULT_FORCE_ANGLE)
                    
                    if angle_key not in force_angle_dict:
                        missing_angles.append(f"({i}, {j}, {k})")
                    
                    ref_mean = angles_ref_mean[idx] if idx < len(angles_ref_mean) and not np.isnan(angles_ref_mean[idx]) else 0.0
                    ref_std = angles_ref_std[idx] if idx < len(angles_ref_std) and not np.isnan(angles_ref_std[idx]) else 0.0
                    sim_std = angles_sim_std[idx] if idx < len(angles_sim_std) and not np.isnan(angles_sim_std[idx]) else 0.0
                    prev_std = angles_prev_std[idx] if idx < len(angles_prev_std) and angles_prev_std[idx] is not None and not np.isnan(angles_prev_std[idx]) else None
                    r2 = angles_r2[idx] if idx < len(angles_r2) else 0.0
                    
                    F_new, status = choose_best_force(FORCE_ANGLE, ref_std, sim_std, prev_std, r2, force_type='angle')
                    
                    if "FORCED TO MIN" in status:
                        forced_angles.append(f"({i}, {j}, {k}) -> {F_new:.1f}")
                    elif "FORCED TO MAX" in status:
                        forced_angles.append(f"({i}, {j}, {k}) -> {F_new:.1f}")
                    
                    f.write(format_angle(i, j, k, ref_mean, F_new) +
                            f" ; σ_ref={ref_std:.4f} σ_sim={sim_std:.4f} "
                            f"σ_prev={prev_std if prev_std is not None else 'NA'} "
                            f"R²={r2:.4f} status={status}\n")
                f.write("\n")
                
            elif section_name == 'dihedrals':
                # Write new dihedrals section
                f.write("[ dihedrals ]\n")
                f.write("; i    j    k    l    funct    angle    force.k\n")
                for idx in range(num_dihedrals):
                    if idx < len(dihedral_quartets):
                        i, j, k, l = dihedral_quartets[idx]
                    else:
                        i, j, k, l = 0, 0, 0, 0
                        warnings.warn(f"Dihedral index {idx} out of range for dihedral_quartets")
                    
                    dihedral_key = (i, j, k, l)
                    FORCE_DIHEDRAL = force_dihedral_dict.get(dihedral_key, DEFAULT_FORCE_DIHEDRAL)
                    
                    if dihedral_key not in force_dihedral_dict:
                        missing_dihedrals.append(f"({i}, {j}, {k}, {l})")
                    
                    ref_mean = dihedrals_ref_mean[idx] if idx < len(dihedrals_ref_mean) and not np.isnan(dihedrals_ref_mean[idx]) else 0.0
                    ref_std = dihedrals_ref_std[idx] if idx < len(dihedrals_ref_std) and not np.isnan(dihedrals_ref_std[idx]) else 0.0
                    sim_std = dihedrals_sim_std[idx] if idx < len(dihedrals_sim_std) and not np.isnan(dihedrals_sim_std[idx]) else 0.0
                    prev_std = dihedrals_prev_std[idx] if idx < len(dihedrals_prev_std) and dihedrals_prev_std[idx] is not None and not np.isnan(dihedrals_prev_std[idx]) else None
                    r2 = dihedrals_r2[idx] if idx < len(dihedrals_r2) else 0.0
                    
                    F_new, status = choose_best_force(FORCE_DIHEDRAL, ref_std, sim_std, prev_std, r2, force_type='dihedral')
                    
                    if "FORCED TO MIN" in status:
                        forced_dihedrals.append(f"({i}, {j}, {k}, {l}) -> {F_new:.1f}")
                    elif "FORCED TO MAX" in status:
                        forced_dihedrals.append(f"({i}, {j}, {k}, {l}) -> {F_new:.1f}")
                    
                    f.write(format_dihedral(i, j, k, l, ref_mean, F_new) +
                            f" ; σ_ref={ref_std:.4f} "
                            f"σ_sim={sim_std:.4f} "
                            f"σ_prev={prev_std if prev_std is not None else 'NA'} "
                            f"R²={r2:.4f} status={status}\n")
                f.write("\n")
                
            else:
                # Preserve all other sections exactly as they were
                # Skip writing the section header if it's already there
                for line_idx, line in enumerate(section_lines):
                    # Skip the header line if it's the section header
                    if line_idx == 0 and line.strip().startswith('[') and line.strip().endswith(']'):
                        f.write(f"{line}\n")
                    else:
                        f.write(f"{line}\n")
                f.write("\n")
    
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
        print(f"\n  Bonds forced to limits: {len(forced_bonds)}")
        if len(forced_bonds) <= 10:
            for fb in forced_bonds:
                print(f"    {fb}")
        else:
            for fb in forced_bonds[:5]:
                print(f"    {fb}")
            print(f"    ... and {len(forced_bonds)-5} more")
    
    if forced_angles:
        print(f"\n  Angles forced to limits: {len(forced_angles)}")
        if len(forced_angles) <= 10:
            for fa in forced_angles:
                print(f"    {fa}")
        else:
            for fa in forced_angles[:5]:
                print(f"    {fa}")
            print(f"    ... and {len(forced_angles)-5} more")
    
    if forced_dihedrals:
        print(f"\n  Dihedrals forced to limits: {len(forced_dihedrals)}")
        if len(forced_dihedrals) <= 10:
            for fd in forced_dihedrals:
                print(f"    {fd}")
        else:
            for fd in forced_dihedrals[:5]:
                print(f"    {fd}")
            print(f"    ... and {len(forced_dihedrals)-5} more")
    
    print(f"\n{'='*60}")
    print(f"ITP file successfully generated: {args.itp_out}")
    print(f"Total bonds processed: {num_bonds}")
    print(f"Total angles processed: {num_angles}")
    print(f"Total dihedrals processed: {num_dihedrals}")
    print(f"Mean R²: {np.mean(all_r2) if all_r2 else 0.0:.6f}")
    print(f"Global R²: {1.0 - global_error:.6f}")
    print(f"{'='*60}\n")

if __name__ == "__main__":
    main()
