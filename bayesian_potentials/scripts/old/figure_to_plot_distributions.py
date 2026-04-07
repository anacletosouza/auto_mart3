#!/usr/bin/env python3
"""
Script para processar distribuições de bonds, angles e dihedrals do GROMACS
e gerar figuras e estatísticas.

Uso: python3 plot_distributions.py
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import seaborn as sns

# Configuração de estilo dos gráficos
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

# Criar diretório figures se não existir
os.makedirs("figures/distances", exist_ok=True)
os.makedirs("figures/angles", exist_ok=True)
os.makedirs("figures/dihedrals", exist_ok=True)

def read_xvg(filename):
    """Lê arquivo .xvg do GROMACS ignorando comentários e headers."""
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(('#', '@')):
                continue
            if line.strip():
                try:
                    values = line.strip().split()
                    # Pode ser 1 ou 2 colunas (tempo, valor ou bin, densidade)
                    if len(values) == 2:
                        data.append([float(values[0]), float(values[1])])
                    elif len(values) == 1:
                        data.append([float(values[0])])
                except ValueError:
                    continue
    return np.array(data)

def plot_time_series(data, title, output_prefix, data_type):
    """Plota série temporal e distribuição."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Série temporal
    if data.shape[1] >= 2:  # Tem dados de tempo
        time = data[:, 0]
        values = data[:, 1]
        axes[0].plot(time, values, 'b-', alpha=0.7, linewidth=0.5)
        axes[0].set_xlabel('Time (ps)')
        axes[0].set_ylabel('Value')
        axes[0].set_title(f'{title} - Time Series')
        axes[0].grid(True, alpha=0.3)
        
        # Estatísticas
        mean_val = np.mean(values)
        std_val = np.std(values)
        axes[0].axhline(y=mean_val, color='r', linestyle='--', 
                       label=f'Mean: {mean_val:.3f}')
        axes[0].axhline(y=mean_val + std_val, color='orange', 
                       linestyle=':', alpha=0.7)
        axes[0].axhline(y=mean_val - std_val, color='orange', 
                       linestyle=':', alpha=0.7)
        axes[0].legend()
        
        # Histograma/Distribuição
        axes[1].hist(values, bins=50, density=True, alpha=0.7, 
                    color='skyblue', edgecolor='black')
        axes[1].set_xlabel('Value')
        axes[1].set_ylabel('Probability Density')
        axes[1].set_title(f'{title} - Distribution')
        axes[1].grid(True, alpha=0.3)
        
        # Fit de distribuição normal
        x = np.linspace(min(values), max(values), 100)
        pdf = stats.norm.pdf(x, mean_val, std_val)
        axes[1].plot(x, pdf, 'r-', linewidth=2, 
                    label=f'Normal fit\nμ={mean_val:.3f}, σ={std_val:.3f}')
        axes[1].legend()
        
    else:
        # Arquivo de distribuição (2 colunas: bin, densidade)
        bins = data[:, 0]
        density = data[:, 1]
        axes[1].bar(bins, density, width=bins[1]-bins[0] if len(bins) > 1 else 0.1, 
                   alpha=0.7, color='skyblue', edgecolor='black')
        axes[1].set_xlabel('Value')
        axes[1].set_ylabel('Probability Density')
        axes[1].set_title(f'{title} - Distribution')
        axes[1].grid(True, alpha=0.3)
        
        # Estatísticas aproximadas
        mean_val = np.sum(bins * density) / np.sum(density) if np.sum(density) > 0 else 0
        var_val = np.sum(((bins - mean_val)**2) * density) / np.sum(density) if np.sum(density) > 0 else 0
        std_val = np.sqrt(var_val)
        
        # Colocar texto com estatísticas
        axes[1].text(0.05, 0.95, f'Mean: {mean_val:.3f}\nStd: {std_val:.3f}', 
                    transform=axes[1].transAxes, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    
    # Salvar em PNG e SVG
    plt.savefig(f"{output_prefix}.png", dpi=150, bbox_inches='tight')
    plt.savefig(f"{output_prefix}.svg", bbox_inches='tight')
    plt.close()
    
    return mean_val, std_val

# ============================================
# Processar Bonds (Distances)
# ============================================
print("Processando Bonds...")
bond_stats = []

# Procurar arquivos bond_X.xvg no diretório bonds_mapped
bond_files = sorted(glob.glob("bonds_mapped/bond_*.xvg"))

for bond_file in bond_files:
    # Extrair índice do nome do arquivo
    idx = int(os.path.basename(bond_file).replace('bond_', '').replace('.xvg', ''))
    
    # Ler dados
    data = read_xvg(bond_file)
    
    if len(data) > 0:
        # Verificar se é série temporal (primeira coluna varia)
        if data.shape[1] == 2 and np.std(data[:,0]) > 0.1:
            mean_val, std_val = plot_time_series(
                data, 
                f'Bond {idx}', 
                f'figures/distances/bond_{idx}', 
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
            # Tentar arquivo de distribuição
            dist_file = f"bonds_mapped/distr_bond_{idx}.xvg"
            if os.path.exists(dist_file):
                dist_data = read_xvg(dist_file)
                if len(dist_data) > 0:
                    mean_val, std_val = plot_time_series(
                        dist_data, 
                        f'Bond {idx} Distribution', 
                        f'figures/distances/bond_{idx}_dist', 
                        'bond'
                    )
                    bond_stats.append({
                        'index': idx,
                        'mean': mean_val,
                        'sd': std_val,
                        'type': 'bond'
                    })
                    print(f"  Bond {idx} (dist): mean={mean_val:.4f}, sd={std_val:.4f}")

# Salvar estatísticas de bonds
if bond_stats:
    df_bonds = pd.DataFrame(bond_stats)
    df_bonds.to_csv('bond_statistics.tsv', sep='\t', index=False)
    print(f"Salvo bond_statistics.tsv com {len(bond_stats)} bonds")

# ============================================
# Processar Angles
# ============================================
print("\nProcessando Angles...")
angle_stats = []

# Procurar arquivos ang_X.xvg no diretório angles_mapped
angle_files = sorted(glob.glob("angles_mapped/ang_*.xvg"))

for angle_file in angle_files:
    idx = int(os.path.basename(angle_file).replace('ang_', '').replace('.xvg', ''))
    
    data = read_xvg(angle_file)
    
    if len(data) > 0:
        if data.shape[1] == 2 and np.std(data[:,0]) > 0.1:
            mean_val, std_val = plot_time_series(
                data, 
                f'Angle {idx}', 
                f'figures/angles/angle_{idx}', 
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
            dist_file = f"angles_mapped/distr_ang_{idx}.xvg"
            if os.path.exists(dist_file):
                dist_data = read_xvg(dist_file)
                if len(dist_data) > 0:
                    mean_val, std_val = plot_time_series(
                        dist_data, 
                        f'Angle {idx} Distribution', 
                        f'figures/angles/angle_{idx}_dist', 
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
    df_angles.to_csv('angle_statistics.tsv', sep='\t', index=False)
    print(f"Salvo angle_statistics.tsv com {len(angle_stats)} angles")

# ============================================
# Processar Dihedrals
# ============================================
print("\nProcessando Dihedrals...")
dihedral_stats = []

# Procurar arquivos dih_X.xvg no diretório dihedrals_mapped
dihedral_files = sorted(glob.glob("dihedrals_mapped/dih_*.xvg"))

for dihedral_file in dihedral_files:
    idx = int(os.path.basename(dihedral_file).replace('dih_', '').replace('.xvg', ''))
    
    data = read_xvg(dihedral_file)
    
    if len(data) > 0:
        if data.shape[1] == 2 and np.std(data[:,0]) > 0.1:
            mean_val, std_val = plot_time_series(
                data, 
                f'Dihedral {idx}', 
                f'figures/dihedrals/dihedral_{idx}', 
                'dihedral'
            )
            dihedral_stats.append({
                'index': idx,
                'mean': mean_val,
                'sd': std_val
            })
            print(f"  Dihedral {idx}: mean={mean_val:.4f}, sd={std_val:.4f}")
        else:
            dist_file = f"dihedrals_mapped/distr_dih_{idx}.xvg"
            if os.path.exists(dist_file):
                dist_data = read_xvg(dist_file)
                if len(dist_data) > 0:
                    mean_val, std_val = plot_time_series(
                        dist_data, 
                        f'Dihedral {idx} Distribution', 
                        f'figures/dihedrals/dihedral_{idx}_dist', 
                        'dihedral'
                    )
                    dihedral_stats.append({
                        'index': idx,
                        'mean': mean_val,
                        'sd': std_val
                    })
                    print(f"  Dihedral {idx} (dist): mean={mean_val:.4f}, sd={std_val:.4f}")

# Salvar dihedral_statistics.tsv (formato específico solicitado)
if dihedral_stats:
    df_dihedrals = pd.DataFrame(dihedral_stats)
    df_dihedrals.to_csv('dihedral_statistics.tsv', sep='\t', index=False, 
                        columns=['index', 'mean', 'sd'])
    print(f"Salvo dihedral_statistics.tsv com {len(dihedral_stats)} dihedrals")

print("\nProcessamento concluído!")
print("Arquivos gerados:")
print("  - bond_statistics.tsv")
print("  - angle_statistics.tsv")
print("  - dihedral_statistics.tsv")
print("\nFiguras salvas em:")
print("  - figures/distances/")
print("  - figures/angles/")
print("  - figures/dihedrals/")

# Verificar se figuras foram criadas
n_figures = len(glob.glob("figures/**/*.png", recursive=True))
print(f"\nTotal de figuras geradas: {n_figures}")
