#!/usr/bin/env python3
"""
Generate Martini-style .itp with automatic force constant adjustment
based on reference and simulated statistics.

FINAL VERSION WITH SIMULATED ANNEALING:
- Bayesian force update (Gibbs-like sampling)
- Simulated Annealing acceptance for global error
- Reuse previous .itp if iteration is rejected
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

# ============================================
# Helper functions
# ============================================

def read_tsv(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} not found.")
    return pd.read_csv(file_path, sep="\t")

def read_ndx(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} not found.")
    
    indices = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('[') or line.startswith(';') or not line:
                continue
            parts = line.split()
            indices.append(tuple(int(x) for x in parts))
    return indices

def format_bond(i,j,length,force):
    return f"  {i:2d}  {j:2d}   1    {length:.3f}    {force:.1f}   ;"

def format_angle(i,j,k,angle,force):
    return f"  {i:2d}  {j:2d}  {k:2d}   2    {angle:.1f}    {force:.1f}   ;"

def format_dihedral(i,j,k,l,angle,force):
    return f"  {i:2d}  {j:2d}  {k:2d}  {l:2d}   2    {angle:.1f}    {force:.1f}   ;"

# ============================================
# Bayesian update (Gibbs-like)
# ============================================

def choose_best_force(F_orig, sigma_ref, sigma_curr, sigma_prev):
    if sigma_curr == 0:
        return F_orig, "NOT CHANGED (zero sigma)"

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
        F_candidates.append(F_orig * (s / sigma_ref))

    if not F_candidates:
        return F_orig, "NOT CHANGED (no samples)"

    best_F = F_orig
    best_error = abs(sigma_ref - sigma_curr)

    for F_cand in F_candidates:
        sigma_pred = sigma_curr * (F_orig / F_cand)
        err = abs(sigma_ref - sigma_pred)
        if err < best_error:
            best_error = err
            best_F = F_cand

    if best_F == F_orig:
        return F_orig, "NOT CHANGED (Bayes reject)"
    else:
        return best_F, "CHANGED (Bayes)"

# ============================================
# Argument parsing
# ============================================

parser = argparse.ArgumentParser()
parser.add_argument("--bonds_ref", required=True)
parser.add_argument("--angles_ref", required=True)
parser.add_argument("--dihedrals_ref", required=True)
parser.add_argument("--bonds_sim", required=True)
parser.add_argument("--angles_sim", required=True)
parser.add_argument("--dihedrals_sim", required=True)
parser.add_argument("--atoms_json", required=True)
parser.add_argument("--ndx_bounds", required=True)
parser.add_argument("--ndx_angles", required=True)
parser.add_argument("--ndx_dihedrals", required=True)
parser.add_argument("--molecule_name", default="FA2")
parser.add_argument("--dihedrals_target", action='store_true')
parser.add_argument("--itp_out", default="FA2_final.itp")

args = parser.parse_args()

# ============================================
# Detect iteration
# ============================================

current_iter = int(os.path.basename(args.bonds_sim).split("_")[-1].split(".")[0])
prev_iter = current_iter - 1

base_dir = os.path.abspath(os.path.join(os.path.dirname(args.bonds_sim), "../../"))
prev_stat_dir = os.path.join(base_dir, f"iter_{prev_iter}", "stat")

# ============================================
# Load data
# ============================================

bonds_ref_df = read_tsv(args.bonds_ref)
angles_ref_df = read_tsv(args.angles_ref)
dihedrals_ref_df = read_tsv(args.dihedrals_ref)

bonds_sim_df = read_tsv(args.bonds_sim)
angles_sim_df = read_tsv(args.angles_sim)
dihedrals_sim_df = read_tsv(args.dihedrals_sim)

def safe_load(file):
    return read_tsv(file) if os.path.exists(file) else None

bonds_prev_df = safe_load(f"{prev_stat_dir}/bond_{prev_iter}.tsv")
angles_prev_df = safe_load(f"{prev_stat_dir}/angle_{prev_iter}.tsv")
dihedrals_prev_df = safe_load(f"{prev_stat_dir}/dihedral_{prev_iter}.tsv")

# ============================================
# Load topology
# ============================================

with open(args.atoms_json) as f:
    atoms_data = json.load(f)

bond_pairs = read_ndx(args.ndx_bounds)
angle_triplets = read_ndx(args.ndx_angles)
dihedral_quartets = read_ndx(args.ndx_dihedrals)

# ============================================
# Global error
# ============================================

def compute_global_error(df_ref, df_sim):
    return ((df_ref['sd'] - df_sim['sd']).abs()).mean()

global_error = (
    compute_global_error(bonds_ref_df, bonds_sim_df) +
    compute_global_error(angles_ref_df, angles_sim_df) +
    compute_global_error(dihedrals_ref_df, dihedrals_sim_df)
)

# ============================================
# Save iteration evolution
# ============================================

inter_dir = os.path.abspath(os.path.join(os.path.dirname(args.itp_out), "../inter"))
os.makedirs(inter_dir, exist_ok=True)

iteration_file = os.path.join(inter_dir, "dihedral_iteration.tsv")

new_row = pd.DataFrame({
    "iteration": [current_iter],
    "global_error": [global_error]
})

# ============================================
# Simulated Annealing acceptance
# ============================================

T0 = 10.0
alpha = 0.85
T = T0 * (alpha ** current_iter)

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

        shutil.copy(prev_itp, args.itp_out)
        print(f"Reused previous ITP: {prev_itp}")
        print(f"Output written to: {args.itp_out}")
        exit(0)

    merged = pd.concat([old, new_row]).drop_duplicates("iteration").sort_values("iteration")
else:
    merged = new_row

merged.to_csv(iteration_file, sep="\t", index=False)

# ============================================
# Plot evolution
# ============================================

plt.figure()
plt.plot(
    merged["iteration"],
    merged["global_error"],
    marker='o',
    markeredgecolor='black'
)
plt.xlabel("iteration")
plt.ylabel("global error")
plt.grid()
plt.savefig(os.path.join(inter_dir, "dihedral_evolution.png"))
plt.savefig(os.path.join(inter_dir, "dihedral_evolution.svg"))
plt.close()

# ============================================
# Generate ITP
# ============================================

FORCE_BOND = 1250
FORCE_ANGLE = 25
FORCE_DIHEDRAL = 25

timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

with open(args.itp_out, 'w') as f:
    f.write(f";;;;;; {args.molecule_name} - Final topology\n")
    f.write(f"; Generated on {timestamp}\n\n")

    f.write("[ moleculetype ]\n")
    f.write(f"  {args.molecule_name} 1\n\n")

    f.write("[ atoms ]\n")
    for bead in atoms_data:
        f.write(f"  {bead['nr']:2d}  {bead['type']}  {bead['resnr']}  {bead['residue']}  {bead['atom']}  {bead['cgnr']}  {bead['charge']}  {bead['mass']}\n")

    # Bonds (FIX)
    f.write("\n[ bonds ]\n")
    for idx, row in bonds_ref_df.iterrows():
        if idx < len(bond_pairs):
            i, j = bond_pairs[idx]
        else:
            i, j = 0, 0

        sigma_sim = bonds_sim_df.loc[idx,'sd'] if idx < len(bonds_sim_df) else row['sd']
        prev_sigma = bonds_prev_df.loc[idx,'sd'] if bonds_prev_df is not None and idx < len(bonds_prev_df) else None

        F_new, status = choose_best_force(FORCE_BOND, row['sd'], sigma_sim, prev_sigma)

        f.write(format_bond(i,j,row['mean'],F_new) +
                f" ; σ_ref={row['sd']:.4f} σ_sim={sigma_sim:.4f} "
                f"σ_prev={prev_sigma if prev_sigma is not None else 'NA'} status={status}\n")

    # Angles (FIX)
    f.write("\n[ angles ]\n")
    for idx, row in angles_ref_df.iterrows():
        if idx < len(angle_triplets):
            i, j, k = angle_triplets[idx]
        else:
            i, j, k = 0, 0, 0

        sigma_sim = angles_sim_df.loc[idx,'sd'] if idx < len(angles_sim_df) else row['sd']
        prev_sigma = angles_prev_df.loc[idx,'sd'] if angles_prev_df is not None and idx < len(angles_prev_df) else None

        F_new, status = choose_best_force(FORCE_ANGLE, row['sd'], sigma_sim, prev_sigma)

        f.write(format_angle(i,j,k,row['mean'],F_new) +
                f" ; σ_ref={row['sd']:.4f} σ_sim={sigma_sim:.4f} "
                f"σ_prev={prev_sigma if prev_sigma is not None else 'NA'} status={status}\n")

    # Dihedrals (FIX)
    f.write("\n[ dihedrals ]\n")
    f.write("; i  j  k  l  funct angle force.k\n")
    for idx, row in dihedrals_ref_df.iterrows():
        if idx < len(dihedral_quartets):
            i, j, k, l = dihedral_quartets[idx]
        else:
            i, j, k, l = 0, 0, 0, 0

        sigma_sim = dihedrals_sim_df.loc[idx,'sd'] if idx < len(dihedrals_sim_df) else row['sd']
        prev_sigma = dihedrals_prev_df.loc[idx,'sd'] if dihedrals_prev_df is not None and idx < len(dihedrals_prev_df) else None

        angle_val = row['mean'] if args.dihedrals_target else 0

        F_new, status = choose_best_force(FORCE_DIHEDRAL, row['sd'], sigma_sim, prev_sigma)

        f.write(format_dihedral(i,j,k,l,angle_val,F_new) +
                f" ; σ_ref={row['sd']:.4f} "
                f"σ_sim={sigma_sim:.4f} "
                f"σ_prev={prev_sigma if prev_sigma is not None else 'NA'} "
                f"status={status}\n")

print(f"\nITP file generated: {args.itp_out}")
