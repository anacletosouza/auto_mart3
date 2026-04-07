import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import seaborn as sns
import re
import argparse


# ============================================
# Argument parser
# ============================================

parser = argparse.ArgumentParser(description="Process XVG files (bonds, angles, dihedrals)")

# Inputs
parser.add_argument("--bonds_dir", default="bonds_mapped")
parser.add_argument("--angles_dir", default="angles_mapped")
parser.add_argument("--dihedrals_dir", default="dihedrals_mapped")

# Outputs (figures)
parser.add_argument("--figures_dir", default="figures")

# Outputs (tables)
parser.add_argument("--bond_out", default="bond_statistics.tsv")
parser.add_argument("--angle_out", default="angle_statistics.tsv")
parser.add_argument("--dihedral_out", default="dihedral_statistics.tsv")

args = parser.parse_args()


# ============================================
# Plot style configuration
# ============================================

plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")


# ============================================
# Create figure directories
# ============================================

dist_dir = os.path.join(args.figures_dir, "distances")
ang_dir = os.path.join(args.figures_dir, "angles")
dih_dir = os.path.join(args.figures_dir, "dihedrals")

os.makedirs(dist_dir, exist_ok=True)
os.makedirs(ang_dir, exist_ok=True)
os.makedirs(dih_dir, exist_ok=True)


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


def plot_time_series(data, title, output_prefix, data_type):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    if data.shape[1] >= 2:
        time = data[:, 0]
        values = data[:, 1]

        axes[0].plot(time, values, 'b-', alpha=0.7, linewidth=0.5)
        axes[0].set_xlabel('Time (ps)')
        axes[0].set_ylabel('Value')
        axes[0].set_title(f'{title} - Time Series')
        axes[0].grid(True, alpha=0.3)

        mean_val = np.mean(values)
        std_val = np.std(values)

        axes[0].axhline(y=mean_val, color='r', linestyle='--',
                        label=f'Mean: {mean_val:.3f}')
        axes[0].axhline(y=mean_val + std_val, color='orange',
                        linestyle=':', alpha=0.7)
        axes[0].axhline(y=mean_val - std_val, color='orange',
                        linestyle=':', alpha=0.7)
        axes[0].legend()

        axes[1].hist(values, bins=50, density=True, alpha=0.7,
                     color='skyblue', edgecolor='black')
        axes[1].set_xlabel('Value')
        axes[1].set_ylabel('Probability Density')
        axes[1].set_title(f'{title} - Distribution')
        axes[1].grid(True, alpha=0.3)

        x = np.linspace(min(values), max(values), 100)
        pdf = stats.norm.pdf(x, mean_val, std_val)
        axes[1].plot(x, pdf, 'r-', linewidth=2,
                     label=f'Normal fit\nμ={mean_val:.3f}, σ={std_val:.3f}')
        axes[1].legend()

    else:
        bins = data[:, 0]
        density = data[:, 1]

        axes[1].bar(
            bins,
            density,
            width=bins[1] - bins[0] if len(bins) > 1 else 0.1,
            alpha=0.7,
            color='skyblue',
            edgecolor='black'
        )

        axes[1].set_xlabel('Value')
        axes[1].set_ylabel('Probability Density')
        axes[1].set_title(f'{title} - Distribution')
        axes[1].grid(True, alpha=0.3)

        mean_val = np.sum(bins * density) / np.sum(density) if np.sum(density) > 0 else 0
        var_val = np.sum(((bins - mean_val)**2) * density) / np.sum(density) if np.sum(density) > 0 else 0
        std_val = np.sqrt(var_val)

        axes[1].text(
            0.05,
            0.95,
            f'Mean: {mean_val:.3f}\nStd: {std_val:.3f}',
            transform=axes[1].transAxes,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        )

    plt.tight_layout()

    plt.savefig(f"{output_prefix}.png", dpi=150, bbox_inches='tight')
    plt.savefig(f"{output_prefix}.svg", bbox_inches='tight')
    plt.close()

    return mean_val, std_val


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

            mean_val, std_val = plot_time_series(
                data,
                f'Bond {idx}',
                f'{dist_dir}/bond_{idx}',
                'bond'
            )

            bond_stats.append({
                'index': idx,
                'mean': mean_val,
                'sd': std_val,
                'type': 'bond'
            })

            print(f"  Bond {idx}: mean={mean_val:.4f}, sd={std_val:.4f}")

        else:

            dist_file = os.path.join(args.bonds_dir, f"distr_bond_{idx}.xvg")

            if os.path.exists(dist_file):

                dist_data = read_xvg(dist_file)

                if len(dist_data) > 0:

                    mean_val, std_val = plot_time_series(
                        dist_data,
                        f'Bond {idx} Distribution',
                        f'{dist_dir}/bond_{idx}_dist',
                        'bond'
                    )

                    bond_stats.append({
                        'index': idx,
                        'mean': mean_val,
                        'sd': std_val,
                        'type': 'bond'
                    })

                    print(f"  Bond {idx} (dist): mean={mean_val:.4f}, sd={std_val:.4f}")


if bond_stats:
    df_bonds = pd.DataFrame(bond_stats)
    df_bonds.to_csv(args.bond_out, sep='\t', index=False)
    print(f"Saved {args.bond_out} with {len(bond_stats)} bonds")


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

            mean_val, std_val = plot_time_series(
                data,
                f'Angle {idx}',
                f'{ang_dir}/angle_{idx}',
                'angle'
            )

            angle_stats.append({
                'index': idx,
                'mean': mean_val,
                'sd': std_val,
                'type': 'angle'
            })

            print(f"  Angle {idx}: mean={mean_val:.4f}, sd={std_val:.4f}")

        else:

            dist_file = os.path.join(args.angles_dir, f"distr_ang_{idx}.xvg")

            if os.path.exists(dist_file):

                dist_data = read_xvg(dist_file)

                if len(dist_data) > 0:

                    mean_val, std_val = plot_time_series(
                        dist_data,
                        f'Angle {idx} Distribution',
                        f'{ang_dir}/angle_{idx}_dist',
                        'angle'
                    )

                    angle_stats.append({
                        'index': idx,
                        'mean': mean_val,
                        'sd': std_val,
                        'type': 'angle'
                    })

                    print(f"  Angle {idx} (dist): mean={mean_val:.4f}, sd={std_val:.4f}")


if angle_stats:
    df_angles = pd.DataFrame(angle_stats)
    df_angles.to_csv(args.angle_out, sep='\t', index=False)
    print(f"Saved {args.angle_out} with {len(angle_stats)} angles")


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

            mean_val, std_val = plot_time_series(
                data,
                f'Dihedral {idx}',
                f'{dih_dir}/dihedral_{idx}',
                'dihedral'
            )

            dihedral_stats.append({
                'index': idx,
                'mean': mean_val,
                'sd': std_val
            })

            print(f"  Dihedral {idx}: mean={mean_val:.4f}, sd={std_val:.4f}")

        else:

            dist_file = os.path.join(args.dihedrals_dir, f"distr_dih_{idx}.xvg")

            if os.path.exists(dist_file):

                dist_data = read_xvg(dist_file)

                if len(dist_data) > 0:

                    mean_val, std_val = plot_time_series(
                        dist_data,
                        f'Dihedral {idx} Distribution',
                        f'{dih_dir}/dihedral_{idx}_dist',
                        'dihedral'
                    )

                    dihedral_stats.append({
                        'index': idx,
                        'mean': mean_val,
                        'sd': std_val
                    })

                    print(f"  Dihedral {idx} (dist): mean={mean_val:.4f}, sd={std_val:.4f}")


if dihedral_stats:
    df_dihedrals = pd.DataFrame(dihedral_stats)
    df_dihedrals.to_csv(
        args.dihedral_out,
        sep='\t',
        index=False,
        columns=['index', 'mean', 'sd']
    )

    print(f"Saved {args.dihedral_out} with {len(dihedral_stats)} dihedrals")


# ============================================
# Final summary
# ============================================

print("\nProcessing completed!")

print("Generated files:")
print(f"  - {args.bond_out}")
print(f"  - {args.angle_out}")
print(f"  - {args.dihedral_out}")

print("\nFigures saved in:")
print(f"  - {dist_dir}")
print(f"  - {ang_dir}")
print(f"  - {dih_dir}")

n_figures = len(glob.glob(os.path.join(args.figures_dir, "**/*.png"), recursive=True))
print(f"\nTotal figures generated: {n_figures}")
