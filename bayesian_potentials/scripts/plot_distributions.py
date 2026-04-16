import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import seaborn as sns
import re
import argparse
from scipy.stats import gaussian_kde
import gc


# ============================================
# Argument parser
# ============================================

parser = argparse.ArgumentParser(description="Process XVG files (bonds, angles, dihedrals)")

# Inputs
parser.add_argument("--bonds_ref_dir", default="bonds")
parser.add_argument("--angles_ref_dir", default="angles")
parser.add_argument("--dihedrals_ref_dir", default="dihedrals")

parser.add_argument("--bonds_sim_dir", default="bonds")
parser.add_argument("--angles_sim_dir", default="angles")
parser.add_argument("--dihedrals_sim_dir", default="dihedrals")

# Outputs (figures)
parser.add_argument("--figures_dir", default="figures")

args = parser.parse_args()


# ============================================
# Plot style configuration
# ============================================

# Use white background
plt.style.use('default')
sns.set_style("white")
sns.set_palette("husl")

# Set font sizes
plt.rcParams['font.size'] = 16
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14

# Limit the number of open figures
plt.rcParams['figure.max_open_warning'] = 0


# ============================================
# Create figure directories for each type
# ============================================

# Bonds directories
bonds_dist_joined_dir = os.path.join(args.figures_dir, "bonds", "distr_ref_simulated_joined")
bonds_dist_separated_ref_dir = os.path.join(args.figures_dir, "bonds", "distr_ref_simulated_separated", "ref")
bonds_dist_separated_sim_dir = os.path.join(args.figures_dir, "bonds", "distr_ref_simulated_separated", "simulated")
bonds_ts_joined_dir = os.path.join(args.figures_dir, "bonds", "time_series_ref_simulated_joined")
bonds_ts_separated_ref_dir = os.path.join(args.figures_dir, "bonds", "time_series_ref_simulated_separated", "ref")
bonds_ts_separated_sim_dir = os.path.join(args.figures_dir, "bonds", "time_series_ref_simulated_separated", "simulated")

# Angles directories
angles_dist_joined_dir = os.path.join(args.figures_dir, "angles", "distr_ref_simulated_joined")
angles_dist_separated_ref_dir = os.path.join(args.figures_dir, "angles", "distr_ref_simulated_separated", "ref")
angles_dist_separated_sim_dir = os.path.join(args.figures_dir, "angles", "distr_ref_simulated_separated", "simulated")
angles_ts_joined_dir = os.path.join(args.figures_dir, "angles", "time_series_ref_simulated_joined")
angles_ts_separated_ref_dir = os.path.join(args.figures_dir, "angles", "time_series_ref_simulated_separated", "ref")
angles_ts_separated_sim_dir = os.path.join(args.figures_dir, "angles", "time_series_ref_simulated_separated", "simulated")

# Dihedrals directories
dihedrals_dist_joined_dir = os.path.join(args.figures_dir, "dihedrals", "distr_ref_simulated_joined")
dihedrals_dist_separated_ref_dir = os.path.join(args.figures_dir, "dihedrals", "distr_ref_simulated_separated", "ref")
dihedrals_dist_separated_sim_dir = os.path.join(args.figures_dir, "dihedrals", "distr_ref_simulated_separated", "simulated")
dihedrals_ts_joined_dir = os.path.join(args.figures_dir, "dihedrals", "time_series_ref_simulated_joined")
dihedrals_ts_separated_ref_dir = os.path.join(args.figures_dir, "dihedrals", "time_series_ref_simulated_separated", "ref")
dihedrals_ts_separated_sim_dir = os.path.join(args.figures_dir, "dihedrals", "time_series_ref_simulated_separated", "simulated")

# Create all directories
all_dirs = [
    bonds_dist_joined_dir, bonds_dist_separated_ref_dir, bonds_dist_separated_sim_dir,
    bonds_ts_joined_dir, bonds_ts_separated_ref_dir, bonds_ts_separated_sim_dir,
    angles_dist_joined_dir, angles_dist_separated_ref_dir, angles_dist_separated_sim_dir,
    angles_ts_joined_dir, angles_ts_separated_ref_dir, angles_ts_separated_sim_dir,
    dihedrals_dist_joined_dir, dihedrals_dist_separated_ref_dir, dihedrals_dist_separated_sim_dir,
    dihedrals_ts_joined_dir, dihedrals_ts_separated_ref_dir, dihedrals_ts_separated_sim_dir
]

for directory in all_dirs:
    os.makedirs(directory, exist_ok=True)


# ============================================
# Helper functions
# ============================================

def extract_index(filename):
    """Extract numeric index from filename"""
    match = re.search(r'\d+', os.path.basename(filename))
    return int(match.group()) if match else 0


def read_xvg(filename):
    """Read XVG file and return data as numpy array"""
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(('#', '@')):
                continue
            if line.strip():
                try:
                    values = line.strip().split()
                    if len(values) >= 2:
                        data.append([float(values[0]), float(values[1])])
                    elif len(values) == 1:
                        data.append([float(values[0])])
                except ValueError:
                    continue
    return np.array(data)


def get_distribution_data(data, n_points=1000):
    """Extract distribution data from time series or histogram"""
    if len(data) == 0:
        return np.array([])
    
    if data.shape[1] >= 2:
        # Time series data
        values = data[:, 1]
        return values
    else:
        # Histogram data (bins, density)
        bins = data[:, 0]
        density = data[:, 1] if data.shape[1] > 1 else np.ones_like(bins)
        
        # Reconstruct points from histogram
        if len(bins) > 1:
            bin_width = bins[1] - bins[0]
            # Generate points based on density
            total_points = 10000
            probabilities = density / np.sum(density) if np.sum(density) > 0 else np.ones_like(density) / len(density)
            selected_bins = np.random.choice(len(bins), size=total_points, p=probabilities)
            values = bins[selected_bins] + np.random.uniform(-bin_width/2, bin_width/2, total_points)
            return values
        else:
            return np.array([])


def plot_distribution_joined(values_ref, values_sim, title, output_path, xlabel='Value'):
    """Plot reference and simulated distributions together (curves only, no histograms)"""
    fig, ax = plt.subplots(figsize=(10, 6))
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    
    # Calculate KDE for reference
    if len(values_ref) > 1:
        try:
            kde_ref = gaussian_kde(values_ref)
            x_min = min(values_ref.min(), values_sim.min()) if len(values_sim) > 0 else values_ref.min()
            x_max = max(values_ref.max(), values_sim.max()) if len(values_sim) > 0 else values_ref.max()
            x_range = np.linspace(x_min, x_max, 1000)
            ax.plot(x_range, kde_ref(x_range), 'b-', linewidth=2.5, 
                   label=f'Reference (n={len(values_ref)})', alpha=0.8)
        except Exception as e:
            print(f"    Warning: Could not compute KDE for reference: {e}")
    
    # Calculate KDE for simulated
    if len(values_sim) > 1:
        try:
            kde_sim = gaussian_kde(values_sim)
            x_min = min(values_ref.min(), values_sim.min()) if len(values_ref) > 0 else values_sim.min()
            x_max = max(values_ref.max(), values_sim.max()) if len(values_ref) > 0 else values_sim.max()
            x_range = np.linspace(x_min, x_max, 1000)
            ax.plot(x_range, kde_sim(x_range), 'r-', linewidth=2.5, 
                   label=f'Simulated (n={len(values_sim)})', alpha=0.8)
        except Exception as e:
            print(f"    Warning: Could not compute KDE for simulated: {e}")
    
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel('Probability Density', fontsize=16)
    ax.set_title(title, fontsize=18)
    
    # Only add legend if there are labels
    if len(values_ref) > 1 or len(values_sim) > 1:
        ax.legend(loc='best', fontsize=14, frameon=True, fancybox=True, shadow=True)
    
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(labelsize=14)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}.png", dpi=150, bbox_inches='tight', facecolor='white')
    plt.savefig(f"{output_path}.svg", bbox_inches='tight', facecolor='white')
    
    # Close the figure and free memory
    plt.close(fig)
    fig = None
    gc.collect()


def plot_distribution_separate(values, title, output_path, xlabel='Value', color='blue'):
    """Plot individual distribution (curve only, no histogram)"""
    fig, ax = plt.subplots(figsize=(10, 6))
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    
    if len(values) > 1:
        try:
            kde = gaussian_kde(values)
            x_range = np.linspace(values.min(), values.max(), 1000)
            
            # Map color name to matplotlib color code
            color_code = 'b' if color == 'blue' else 'r' if color == 'red' else 'g'
            ax.plot(x_range, kde(x_range), f'{color_code}-', linewidth=2.5, 
                   label=f'n={len(values)}', alpha=0.8)
        except Exception as e:
            print(f"    Warning: Could not compute KDE: {e}")
    
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel('Probability Density', fontsize=16)
    ax.set_title(title, fontsize=18)
    
    if len(values) > 1:
        ax.legend(loc='best', fontsize=14, frameon=True, fancybox=True, shadow=True)
    
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(labelsize=14)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}.png", dpi=150, bbox_inches='tight', facecolor='white')
    plt.savefig(f"{output_path}.svg", bbox_inches='tight', facecolor='white')
    
    # Close the figure and free memory
    plt.close(fig)
    fig = None
    gc.collect()


def plot_time_series_joined(data_ref, data_sim, title, output_path, ylabel='Value'):
    """Plot reference and simulated time series together"""
    fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
    fig.patch.set_facecolor('white')
    
    has_content = False
    
    # Reference time series
    if len(data_ref) > 0 and data_ref.shape[1] >= 2:
        time_ref = data_ref[:, 0]
        values_ref = data_ref[:, 1]
        
        axes[0].plot(time_ref, values_ref, 'b-', alpha=0.7, linewidth=1.0)
        axes[0].set_ylabel(ylabel, fontsize=16)
        axes[0].set_title(f'{title} - Reference', fontsize=18)
        axes[0].grid(True, alpha=0.3, linestyle='--')
        axes[0].set_facecolor('white')
        
        # Add mean line
        mean_ref = np.mean(values_ref)
        axes[0].axhline(y=mean_ref, color='r', linestyle='--', linewidth=2,
                       label=f'Mean: {mean_ref:.3f}', alpha=0.7)
        axes[0].legend(fontsize=14, frameon=True, fancybox=True, shadow=True)
        axes[0].tick_params(labelsize=14)
        has_content = True
    
    # Simulated time series
    if len(data_sim) > 0 and data_sim.shape[1] >= 2:
        time_sim = data_sim[:, 0]
        values_sim = data_sim[:, 1]
        
        axes[1].plot(time_sim, values_sim, 'r-', alpha=0.7, linewidth=1.0)
        axes[1].set_xlabel('Time (ps)', fontsize=16)
        axes[1].set_ylabel(ylabel, fontsize=16)
        axes[1].set_title(f'{title} - Simulated', fontsize=18)
        axes[1].grid(True, alpha=0.3, linestyle='--')
        axes[1].set_facecolor('white')
        
        # Add mean line
        mean_sim = np.mean(values_sim)
        axes[1].axhline(y=mean_sim, color='b', linestyle='--', linewidth=2,
                       label=f'Mean: {mean_sim:.3f}', alpha=0.7)
        axes[1].legend(fontsize=14, frameon=True, fancybox=True, shadow=True)
        axes[1].tick_params(labelsize=14)
        has_content = True
    
    if has_content:
        plt.tight_layout()
        plt.savefig(f"{output_path}.png", dpi=150, bbox_inches='tight', facecolor='white')
        plt.savefig(f"{output_path}.svg", bbox_inches='tight', facecolor='white')
    
    # Close the figure and free memory
    plt.close(fig)
    fig = None
    gc.collect()


def plot_time_series_separate(data, title, output_path, ylabel='Value', color='blue'):
    """Plot individual time series"""
    fig, ax = plt.subplots(figsize=(14, 6))
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    
    mean_val, std_val = 0, 0
    
    if len(data) > 0 and data.shape[1] >= 2:
        time = data[:, 0]
        values = data[:, 1]
        
        # Map color name to matplotlib color code
        color_code = 'b' if color == 'blue' else 'r' if color == 'red' else 'g'
        ax.plot(time, values, f'{color_code}-', alpha=0.7, linewidth=1.0)
        ax.set_xlabel('Time (ps)', fontsize=16)
        ax.set_ylabel(ylabel, fontsize=16)
        ax.set_title(title, fontsize=18)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.tick_params(labelsize=14)
        
        # Add mean line
        mean_val = np.mean(values)
        std_val = np.std(values)
        ax.axhline(y=mean_val, color='r', linestyle='--', linewidth=2,
                  label=f'Mean: {mean_val:.3f}', alpha=0.7)
        ax.legend(fontsize=14, frameon=True, fancybox=True, shadow=True)
        
        plt.tight_layout()
        plt.savefig(f"{output_path}.png", dpi=150, bbox_inches='tight', facecolor='white')
        plt.savefig(f"{output_path}.svg", bbox_inches='tight', facecolor='white')
    
    # Close the figure and free memory
    plt.close(fig)
    fig = None
    gc.collect()
    
    return mean_val, std_val


def process_type(data_type, ref_dir, sim_dir, 
                 dist_joined_path, dist_sep_ref_path, dist_sep_sim_path,
                 ts_joined_path, ts_sep_ref_path, ts_sep_sim_path,
                 xlabel, ylabel):
    """Process a specific type (bonds, angles, dihedrals)"""
    
    print(f"\nProcessing {data_type.upper()}...")
    print(f"  Reference dir: {ref_dir}")
    print(f"  Simulated dir: {sim_dir}")
    
    # Determine file prefix
    if data_type == 'bond':
        prefix = 'bond'
    elif data_type == 'angle':
        prefix = 'ang'
    elif data_type == 'dihedral':
        prefix = 'dih'
    else:
        prefix = data_type[:3]
    
    # Find all files
    ref_pattern = os.path.join(ref_dir, f"{prefix}_*.xvg")
    sim_pattern = os.path.join(sim_dir, f"{prefix}_*.xvg")
    
    ref_files = sorted(glob.glob(ref_pattern), key=extract_index)
    sim_files = sorted(glob.glob(sim_pattern), key=extract_index)
    
    print(f"  Found {len(ref_files)} reference files, {len(sim_files)} simulated files")
    
    # Create mapping of indices to files
    ref_dict = {extract_index(f): f for f in ref_files}
    sim_dict = {extract_index(f): f for f in sim_files}
    
    # Process common indices
    common_indices = set(ref_dict.keys()) & set(sim_dict.keys())
    
    if not common_indices:
        print(f"  Warning: No common indices found for {data_type}")
        return
    
    total = len(common_indices)
    for i, idx in enumerate(sorted(common_indices)):
        print(f"  Processing {data_type} {idx}... ({i+1}/{total})")
        
        # Read data
        ref_data = read_xvg(ref_dict[idx])
        sim_data = read_xvg(sim_dict[idx])
        
        if len(ref_data) == 0:
            print(f"    Warning: Empty reference data for index {idx}, skipping...")
            continue
        
        if len(sim_data) == 0:
            print(f"    Warning: Empty simulated data for index {idx}, skipping...")
            continue
        
        # Get distribution data
        ref_values = get_distribution_data(ref_data)
        sim_values = get_distribution_data(sim_data)
        
        # Plot distributions
        if len(ref_values) > 1 or len(sim_values) > 1:
            # Joined distribution
            plot_distribution_joined(
                ref_values, sim_values,
                f'{data_type.capitalize()} {idx} - Reference vs Simulated',
                os.path.join(dist_joined_path, f'{data_type}_{idx}'),
                xlabel
            )
            
            # Separate distributions
            if len(ref_values) > 1:
                plot_distribution_separate(
                    ref_values,
                    f'{data_type.capitalize()} {idx} - Reference',
                    os.path.join(dist_sep_ref_path, f'{data_type}_{idx}'),
                    xlabel,
                    'blue'
                )
            
            if len(sim_values) > 1:
                plot_distribution_separate(
                    sim_values,
                    f'{data_type.capitalize()} {idx} - Simulated',
                    os.path.join(dist_sep_sim_path, f'{data_type}_{idx}'),
                    xlabel,
                    'red'
                )
        else:
            print(f"    Warning: Insufficient data for distribution plots for {data_type} {idx}")
        
        # Plot time series (only if time series data is available)
        if ref_data.shape[1] >= 2 and sim_data.shape[1] >= 2:
            # Joined time series
            plot_time_series_joined(
                ref_data, sim_data,
                f'{data_type.capitalize()} {idx}',
                os.path.join(ts_joined_path, f'{data_type}_{idx}'),
                ylabel
            )
            
            # Separate time series
            plot_time_series_separate(
                ref_data,
                f'{data_type.capitalize()} {idx} - Reference Time Series',
                os.path.join(ts_sep_ref_path, f'{data_type}_{idx}'),
                ylabel,
                'blue'
            )
            
            plot_time_series_separate(
                sim_data,
                f'{data_type.capitalize()} {idx} - Simulated Time Series',
                os.path.join(ts_sep_sim_path, f'{data_type}_{idx}'),
                ylabel,
                'red'
            )
        else:
            print(f"    Note: Time series not available for {data_type} {idx} (histogram data only)")


# ============================================
# Main processing
# ============================================

print("="*60)
print("Starting distribution plotting")
print("="*60)
print(f"Font size: 16 (global)")
print(f"Background: White")
print(f"Figures will be organized by type (bonds/angles/dihedrals)")
print("="*60)

# Process Bonds
process_type(
    'bond',
    args.bonds_ref_dir, args.bonds_sim_dir,
    bonds_dist_joined_dir, bonds_dist_separated_ref_dir, bonds_dist_separated_sim_dir,
    bonds_ts_joined_dir, bonds_ts_separated_ref_dir, bonds_ts_separated_sim_dir,
    'Distance (nm)', 'Distance (nm)'
)

# Process Angles
process_type(
    'angle',
    args.angles_ref_dir, args.angles_sim_dir,
    angles_dist_joined_dir, angles_dist_separated_ref_dir, angles_dist_separated_sim_dir,
    angles_ts_joined_dir, angles_ts_separated_ref_dir, angles_ts_separated_sim_dir,
    'Angle (degrees)', 'Angle (degrees)'
)

# Process Dihedrals
process_type(
    'dihedral',
    args.dihedrals_ref_dir, args.dihedrals_sim_dir,
    dihedrals_dist_joined_dir, dihedrals_dist_separated_ref_dir, dihedrals_dist_separated_sim_dir,
    dihedrals_ts_joined_dir, dihedrals_ts_separated_ref_dir, dihedrals_ts_separated_sim_dir,
    'Dihedral Angle (degrees)', 'Dihedral Angle (degrees)'
)


# ============================================
# Final summary
# ============================================

print("\n" + "="*60)
print("Processing completed!")
print("="*60)

print("\nGenerated directories:")
print("\nBONDS:")
print(f"  - {bonds_dist_joined_dir}")
print(f"  - {bonds_dist_separated_ref_dir}")
print(f"  - {bonds_dist_separated_sim_dir}")
print(f"  - {bonds_ts_joined_dir}")
print(f"  - {bonds_ts_separated_ref_dir}")
print(f"  - {bonds_ts_separated_sim_dir}")

print("\nANGLES:")
print(f"  - {angles_dist_joined_dir}")
print(f"  - {angles_dist_separated_ref_dir}")
print(f"  - {angles_dist_separated_sim_dir}")
print(f"  - {angles_ts_joined_dir}")
print(f"  - {angles_ts_separated_ref_dir}")
print(f"  - {angles_ts_separated_sim_dir}")

print("\nDIHEDRALS:")
print(f"  - {dihedrals_dist_joined_dir}")
print(f"  - {dihedrals_dist_separated_ref_dir}")
print(f"  - {dihedrals_dist_separated_sim_dir}")
print(f"  - {dihedrals_ts_joined_dir}")
print(f"  - {dihedrals_ts_separated_ref_dir}")
print(f"  - {dihedrals_ts_separated_sim_dir}")

# Count figures
n_figures = 0
for directory in all_dirs:
    if os.path.exists(directory):
        n_figures += len(glob.glob(os.path.join(directory, "*.png")))

print(f"\nTotal figures generated: {n_figures}")
print("="*60)
