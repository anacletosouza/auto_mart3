#!/usr/bin/env python3
"""
Command-line interface for Auto_Mart3 package.

This module provides a unified CLI for the Auto_Mart3 package, offering tools for:
- Mapping atomistic trajectories to coarse-grained (CG) resolution
- Generating CG topologies and force field files
- Analyzing bond, angle, and dihedral distributions
- Optimizing CG force constants using Bayesian methods
- Setting up CG simulation systems

"""

import argparse
import sys
import os
import subprocess
from pathlib import Path

try:
    import setproctitle
    setproctitle.setproctitle("auto_mart3")
except ImportError:
    pass 


def get_package_dir():
    """Get the package installation directory.
    
    Returns:
        Path: Path object pointing to the package directory
    """
    return Path(__file__).parent.resolve()

def get_scripts_dir():
    """Get the scripts directory inside the package.
    
    Returns:
        Path: Path object pointing to the scripts directory
    """
    return get_package_dir() / "scripts"
    
def get_data_dir():
    """Get the data directory inside the package.
    
    Returns:
        Path: Path object pointing to the scripts directory
    """
    return get_package_dir() / "data"

def get_bin_dir():
    """Get the bin directory inside the package.
    
    Returns:
        Path: Path object pointing to the bin directory
    """
    return get_package_dir() / "bin"


def get_default_mdp():
    """Get the default MDP file from the package data directory.
    
    Returns:
        str or None: Path to default MDP file if found, None otherwise
    """
    data_dir = get_data_dir()
    if data_dir:
        mdp_file = data_dir / "mdp" / "minimization.mdp"
        if mdp_file.exists():
            return str(mdp_file)
    
    return None

def get_default_ff_dir():
    """Get default force field directory from data directory.
    
    Returns:
        str or None: Path to force field directory if found, None otherwise
    """
    data_dir = get_data_dir()
    if data_dir:
        ff_dir = data_dir / "ff_files"
        if ff_dir.exists():
            return str(ff_dir)
    return None

def find_script(script_name):
    """Find a script in multiple possible locations.
    
    Args:
        script_name (str): Name of the script to find
        
    Returns:
        Path or None: Path to script if found, None otherwise
    """
    # List of possible locations
    possible_locations = [
        Path.cwd() / "bin" / script_name,
        Path.cwd() / script_name,
        get_bin_dir() / script_name,
        Path(__file__).parent.parent / "bin" / script_name,
        Path(__file__).parent.parent / script_name,
        Path.home() / ".local" / "bin" / script_name,
    ]
    
    # Also check if script is in PATH
    import shutil
    which_script = shutil.which(script_name)
    if which_script:
        possible_locations.append(Path(which_script))
    
    for location in possible_locations:
        if location.exists() and os.access(location, os.X_OK):
            return location
    
    return None

def run_mapping_script(args):
    """Run the mapping script with proper module imports.
    
    This function maps atomistic trajectories to coarse-grained resolution
    using predefined bead definitions.
    
    Args:
        args: Parsed command-line arguments containing mapping parameters
    """
    from auto_mart3.scripts.map_aa_to_cg import main as map_main
    
    sys_argv = ["map_aa_to_cg.py"]
    
    # Map arguments for the mapping script
    for key, value in vars(args).items():
        if key == "command" or value is None:
            continue
        
        # Handle the special case: cg_ndx should become index_cg for mapping script
        if key == "cg_ndx":
            sys_argv.append("--index_cg")
            sys_argv.append(str(value))
        elif key == "gmx_cmd":
            sys_argv.append("--gmx_cmd")   
            sys_argv.append(str(value))
        else:
            cmd_key = "--" + key
            if isinstance(value, bool):
                if value:
                    sys_argv.append(cmd_key)
            else:
                sys_argv.append(cmd_key)
                sys_argv.append(str(value))
    
    original_argv = sys.argv
    sys.argv = sys_argv
    
    try:
        map_main()
    finally:
        sys.argv = original_argv

def run_topology_script(args):
    """Run the topology generation script.
    
    Generates GROMACS topology files (.top) for CG simulations including
    force field definitions, molecule types, and system setup.
    
    Args:
        args: Parsed command-line arguments containing topology parameters
    """
    from auto_mart3.scripts.generate_cg_top import main as top_main
    
    sys_argv = ["generate_cg_top.py"]
    
    for key, value in vars(args).items():
        if key == "command" or value is None:
            continue
        
        cmd_key = "--" + key
        
        if isinstance(value, bool):
            if value:
                sys_argv.append(cmd_key)
        else:
            sys_argv.append(cmd_key)
            sys_argv.append(str(value))
    
    original_argv = sys.argv
    sys.argv = sys_argv
    
    try:
        top_main()
    finally:
        sys.argv = original_argv

def run_analysis_script(args):
    """Run the bonds/angles/dihedrals analysis script.
    
    Analyzes MD trajectories to compute distributions of bonds, angles,
    and dihedrals for parameterization.
    
    Args:
        args: Parsed command-line arguments containing analysis parameters
    """
    from auto_mart3.scripts.generate_bonds_angles_dihedrals import main as analysis_main
    
    sys_argv = ["generate_bonds_angles_dihedrals.py"]
    
    for key, value in vars(args).items():
        if key == "command" or value is None:
            continue
        
        # Skip internal flags
        if key in ["skip_grompp", "verbose"]:
            continue
        
        # Generate_bonds_angles_dihedrals.py expects underscores, not hyphens
        cmd_key = "--" + key
        
        if isinstance(value, bool):
            if value:
                sys_argv.append(cmd_key)
        else:
            sys_argv.append(cmd_key)
            sys_argv.append(str(value))
    
    original_argv = sys.argv
    sys.argv = sys_argv
    
    try:
        analysis_main()
    finally:
        sys.argv = original_argv

def run_distribution_script(args):
    """Run the bond/angle/dihedral distribution statistics script.
    
    Computes statistical properties (mean, std dev, etc.) from distribution
    XVG files generated by trajectory analysis.
    
    Args:
        args: Parsed command-line arguments containing distribution parameters
    """
    from auto_mart3.scripts.bp_distributions import main as bp_main
    
    sys_argv = ["bp_distributions.py"]
    
    for key, value in vars(args).items():
        if key == "command" or value is None:
            continue
        
        # Skip internal flags
        if key in ["skip_grompp", "verbose"]:
            continue
        
        # bp_distributions.py expects underscores, not hyphens
        cmd_key = "--" + key
        
        if isinstance(value, bool):
            if value:
                sys_argv.append(cmd_key)
        else:
            sys_argv.append(cmd_key)
            sys_argv.append(str(value))
    
    original_argv = sys.argv
    sys.argv = sys_argv
    
    try:
        bp_main()
    finally:
        sys.argv = original_argv

def run_adaptation_script(args):
    """Run the ITP adaptation script to match GRO atom names.
    
    Updates atom names in ITP topology files to match the naming convention
    used in GRO structure files, ensuring compatibility.
    
    Args:
        args: Parsed command-line arguments containing adaptation parameters
    """
    from auto_mart3.scripts.adaptation_gro_itp import main as adapt_main
    
    sys_argv = ["adaptation_gro_itp.py"]
    
    for key, value in vars(args).items():
        if key == "command" or value is None:
            continue
        
        # Skip internal flags
        if key in ["verbose"]:
            continue
        
        cmd_key = "--" + key
        
        if isinstance(value, bool):
            if value:
                sys_argv.append(cmd_key)
        else:
            sys_argv.append(cmd_key)
            sys_argv.append(str(value))
    
    # Add verbose flag if needed
    if getattr(args, 'verbose', False):
        sys_argv.append("--verbose")
    
    original_argv = sys.argv
    sys.argv = sys_argv
    
    try:
        adapt_main()
    finally:
        sys.argv = original_argv

def run_cg_setup_script(args):
    """Run the CG system setup pipeline script.
    
    Prepares a complete CG simulation system including solvation,
    ion addition, and energy minimization setup.
    
    Args:
        args: Parsed command-line arguments containing system setup parameters
    """
    from auto_mart3.scripts.bp_prep import main as cg_setup_main
    
    sys_argv = ["bp_prep.py"]
    
    for key, value in vars(args).items():
        if key == "command" or value is None:
            continue
        
        # Skip internal flags
        if key in ["verbose"]:
            continue
        
        cmd_key = "--" + key
        
        if isinstance(value, bool):
            if value:
                sys_argv.append(cmd_key)
        else:
            sys_argv.append(cmd_key)
            sys_argv.append(str(value))
    
    # Add verbose flag if needed
    if getattr(args, 'verbose', False):
        sys_argv.append("--verbose")
    
    original_argv = sys.argv
    sys.argv = sys_argv
    
    try:
        cg_setup_main()
    finally:
        sys.argv = original_argv

def run_potential_adjustment_script(args):
    """Run the force constant adjustment script using Bayesian update and simulated annealing with R² correction.
    
    Optimizes CG force constants by comparing reference distributions from
    AA simulations with CG simulations using Bayesian inference and
    simulated annealing optimization.
    
    Args:
        args: Parsed command-line arguments containing optimization parameters
    """
    from auto_mart3.scripts.potential_adjustment import main as force_main
    
    sys_argv = ["potential_adjustment.py"]
    
    # Map CLI arguments (with underscores) to script arguments (with underscores)
    # The CLI uses underscores, the script uses underscores internally
    arg_mapping = {
        # Required XVG directories
        'bonds_ref_xvg_dir': '--bonds_ref_xvg_dir',
        'angles_ref_xvg_dir': '--angles_ref_xvg_dir',
        'dihedrals_ref_xvg_dir': '--dihedrals_ref_xvg_dir',
        'bonds_sim_xvg_dir': '--bonds_sim_xvg_dir',
        'angles_sim_xvg_dir': '--angles_sim_xvg_dir',
        'dihedrals_sim_xvg_dir': '--dihedrals_sim_xvg_dir',
        
        # Required input files
        'itp_cg': '--itp_cg',
        'ndx_bounds': '--ndx_bounds',
        'ndx_angles': '--ndx_angles',
        'ndx_dihedrals': '--ndx_dihedrals',
        
        # Optional parameters
        'molecule_name': '--molecule_name',
        'distribution_points': '--distribution_points',
        'T0': '--T0',
        'alpha': '--alpha',
        'itp_out': '--itp_out',
        
        # XVG prefix options
        'prefix_xvg_bond_ref': '--prefix_xvg_bond_ref',
        'prefix_xvg_angle_ref': '--prefix_xvg_angle_ref',
        'prefix_xvg_dihedral_ref': '--prefix_xvg_dihedral_ref',
        'prefix_xvg_bond_sim': '--prefix_xvg_bond_sim',
        'prefix_xvg_angle_sim': '--prefix_xvg_angle_sim',
        'prefix_xvg_dihedral_sim': '--prefix_xvg_dihedral_sim',
        
        # Force limits
        'min_force_bond': '--min_force_bond',
        'min_force_angle': '--min_force_angle',
        'min_force_dihedral': '--min_force_dihedral',
        'max_force_bond': '--max_force_bond',
        'max_force_angle': '--max_force_angle',
        'max_force_dihedral': '--max_force_dihedral',
        'default_force_bond': '--default_force_bond',
        'default_force_angle': '--default_force_angle',
        'default_force_dihedral': '--default_force_dihedral',
    }
    
    for key, value in vars(args).items():
        if key == "command" or value is None:
            continue
        
        # Skip internal flags
        if key in ["verbose"]:
            continue
        
        if key not in arg_mapping:
            if getattr(args, 'verbose', False):
                print(f"Warning: Unknown argument '{key}' will be ignored")
            continue
        
        cmd_key = arg_mapping[key]
        
        if isinstance(value, bool):
            if value:
                sys_argv.append(cmd_key)
        else:
            sys_argv.append(cmd_key)
            sys_argv.append(str(value))
    
    original_argv = sys.argv
    sys.argv = sys_argv
    
    try:
        force_main()
    finally:
        sys.argv = original_argv
        
def run_plot_distributions_script(args):
    """Run the distribution plotting script.
    
    Generates publication-quality plots comparing reference and simulated
    distributions for bonds, angles, and dihedrals.
    
    Args:
        args: Parsed command-line arguments containing plotting parameters
    """
    from auto_mart3.scripts.plot_distributions import main as plot_main
    
    sys_argv = ["plot_distributions.py"]
    
    for key, value in vars(args).items():
        if key == "command" or value is None:
            continue
        
        # Skip internal flags
        if key in ["verbose"]:
            continue
        
        cmd_key = "--" + key
        
        if isinstance(value, bool):
            if value:
                sys_argv.append(cmd_key)
        else:
            sys_argv.append(cmd_key)
            sys_argv.append(str(value))
    
    original_argv = sys.argv
    sys.argv = sys_argv
    
    try:
        plot_main()
    finally:
        sys.argv = original_argv
        

def run_shell_script_AA(args):
    """Run the main shell pipeline script for AA to CG.
    
    Executes the complete atomistic to coarse-grained parameterization
    pipeline including mapping, topology generation, and system setup.
    
    Args:
        args: Parsed command-line arguments containing pipeline parameters
    """
    script_path = find_script("auto_mart_AA.sh")
    
    if not script_path:
        print(f"Error: auto_mart_AA.sh not found")
        print("Please ensure the script is in one of these locations:")
        print("  - ./bin/auto_mart_AA.sh")
        print("  - ./auto_mart_AA.sh")
        print(f"  - {get_bin_dir()}/auto_mart_AA.sh")
        sys.exit(1)
    
    print(f"Using script: {script_path}")
    
    # Make sure script is executable
    os.chmod(script_path, 0o755)
    
    # Build command arguments for shell script
    cmd_args = [str(script_path)]
    
    # Map Python argument names (with underscores) to shell script flags (with underscores)
    arg_mapping = {
        'aa_tpr': '--aa_tpr',
        'aa_xtc': '--aa_xtc',
        'aa_gro': '--aa_gro',
        'aa_itp': '--aa_itp',
        'beads_json': '--beads_json',
        'input_mdp': '--input_mdp',
        'path_ff': '--path_ff',
        'output_dir': '--output_dir',
        'name_molecule': '--name_molecule',
        'force_application': '--force_application',
        'beads_position': '--beads_position',
        'cycle_restr': '--cycle_restr',
        'maxwarn': '--maxwarn',
        'distance_from_atom': '--distance_from_atom',
        'salt': '--salt',
    }
    
    # Arguments to ignore
    ignore_args = {
        'verbose', 'skip_adapt_itp', 'skip_grompp', 'skip_analysis',
        'skip_distributions', 'keep_intermediate', 'keep_temp',
        'analyze_remove_pbc', 'analyze_group_1', 'analyze_group_2',
        'run_distributions', 'dist_bonds_dir', 'dist_angles_dir',
        'dist_dihedrals_dir', 'dist_output_dir', 'dist_bond_out',
        'dist_angle_out', 'dist_dihedral_out', 'number_molecule',
        'ff', 'ions', 'solvent', 'title_comments', 'title_system',
        'output_topol', 'default_martini', 'remove_pbc', 'no_pbc',
        'command'
    }
    
    for key, value in vars(args).items():
        if key in ignore_args or value is None:
            continue
        
        if key not in arg_mapping:
            if getattr(args, 'verbose', False):
                print(f"Warning: Unknown argument '{key}' will be ignored")
            continue
        
        cmd_key = arg_mapping[key]
        
        if isinstance(value, bool):
            if value:
                cmd_args.append(cmd_key)
        else:
            cmd_args.append(cmd_key)
            cmd_args.append(str(value))
    
    # Debug: print the command if verbose
    if getattr(args, 'verbose', False):
        print(f"\nRunning command:")
        print(' '.join(cmd_args))
        print()
    
    # Run the script
    try:
        result = subprocess.run(cmd_args, env=os.environ.copy())
        sys.exit(result.returncode)
    except FileNotFoundError as e:
        print(f"Error: Failed to execute {script_path}")
        print(f"Details: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nInterrupted by user")
        sys.exit(1)

def run_shell_script_CG(args):
    """Run the CG parameter optimization pipeline script.
    
    Executes iterative Bayesian optimization of CG force constants by
    running CG simulations, comparing with reference distributions,
    and updating parameters automatically.
    
    Args:
        args: Parsed command-line arguments containing optimization pipeline parameters
    """
    # First, try to find auto_mart_CG.sh in various locations
    script_path = find_script("auto_mart_CG.sh")
    
    # If not found, check if user provided scripts_path
    if not script_path and hasattr(args, 'scripts_path') and args.scripts_path:
        # Check if scripts_path is a directory containing the script
        scripts_path = Path(args.scripts_path)
        if scripts_path.is_dir():
            potential_script = scripts_path / "auto_mart_CG.sh"
            if potential_script.exists():
                script_path = potential_script
        elif scripts_path.is_file() and scripts_path.name == "auto_mart_CG.sh":
            script_path = scripts_path
    
    if not script_path:
        print(f"Error: auto_mart_CG.sh not found")
        print("Please ensure the script is in one of these locations:")
        print("  - ./bin/auto_mart_CG.sh")
        print("  - ./auto_mart_CG.sh")
        print(f"  - {get_bin_dir()}/auto_mart_CG.sh")
        print("  - Or provide --scripts_path pointing to the directory containing auto_mart_CG.sh")
        sys.exit(1)
    
    print(f"Using script: {script_path}")
    
    # Make sure script is executable
    os.chmod(script_path, 0o755)
    
    # Build command arguments for shell script
    cmd_args = [str(script_path)]
    
    # Map CLI argument names to shell script flags
    arg_mapping = {
        # Required arguments
        'input_auto_mart_aa_dir': '--INPUT_AUTO_MART_AA_DIR',
        'output_auto_mart_cg_dir': '--OUTPUT_AUTO_MART_CG_DIR',
        
        # Simulation Parameters
        'ntomp': '--ntomp',
        'ntmpi': '--ntmpi',
        'ref_t': '--ref_t',
        'ref_p': '--ref_p',
        'dt_nvt_ps': '--dt_nvt_ps',
        'time_nvt_ps': '--time_nvt_ps',
        'dt_npt_ps': '--dt_npt_ps',
        'time_npt_ps': '--time_npt_ps',
        'dt_md_ps': '--dt_md_ps',
        'time_md_ps': '--time_md_ps',
        
        # File Names
        'topol_cg_file': '--topol_cg_file',
        'solv_ions_gro': '--solv_ions_gro',
        'em_tpr': '--em_tpr',
        'itp_to_optimize': '--itp_to_optimize',
        
        # Force Constants Limits
        'min_force_bond': '--min_force_bond',
        'max_force_bond': '--max_force_bond',
        'default_force_bond': '--default_force_bond',
        'min_force_angle': '--min_force_angle',
        'max_force_angle': '--max_force_angle',
        'default_force_angle': '--default_force_angle',
        'min_force_dihedral': '--min_force_dihedral',
        'max_force_dihedral': '--max_force_dihedral',
        'default_force_dihedral': '--default_force_dihedral',
        
        # Optimization Parameters
        'T0': '--T0',
        'alpha': '--alpha',
        'distribution_points': '--distribution_points',
        'n_iter': '--n_iter',
        
        # Input/Output Paths
        'bonds_ndx': '--bonds_ndx',
        'angles_ndx': '--angles_ndx',
        'dihedrals_ndx': '--dihedrals_ndx',
        'bonds_ref_xvg_dir': '--bonds_ref_xvg_dir',
        'angles_ref_xvg_dir': '--angles_ref_xvg_dir',
        'dihedrals_ref_xvg_dir': '--dihedrals_ref_xvg_dir',
        
        # Prefixes for XVG Files
        'prefix_xvg_bond_ref': '--prefix_xvg_bond_ref',
        'prefix_xvg_angle_ref': '--prefix_xvg_angle_ref',
        'prefix_xvg_dihedral_ref': '--prefix_xvg_dihedral_ref',
        'prefix_xvg_bond_sim': '--prefix_xvg_bond_sim',
        'prefix_xvg_angle_sim': '--prefix_xvg_angle_sim',
        'prefix_xvg_dihedral_sim': '--prefix_xvg_dihedral_sim',
        
        # GROMACS Processing
        'group_2': '--group_2',
        'group_1': '--group_1',
        'index': '--index',
        
        # GROMACS MDRUN Options (GPU/CPU configuration)
        'pin_on': '--pin_on',
        'pin_off': '--pin_off',
        'nb_gpu': '--nb_gpu',
        'nb_cpu': '--nb_cpu',
        'nb_none': '--nb_none',
        'pme_gpu': '--pme_gpu',
        'pme_cpu': '--pme_cpu',
        'pme_none': '--pme_none',
        'bonded_gpu': '--bonded_gpu',
        'bonded_cpu': '--bonded_cpu',
        'bonded_none': '--bonded_none',
        'npme': '--npme',
        'cuda_visible_devices': '--cuda_visible_devices',
        'gpu_id': '--gpu_id',
        
        # Other
        'molecule_name': '--molecule_name',
    }
    
    # Boolean flags (no value)
    bool_flags = [
        'no_remove_gmx_files_in_iter',
        'pin_on',
        'pin_off',
        'nb_gpu',
        'nb_cpu',
        'nb_none',
        'pme_gpu',
        'pme_cpu',
        'pme_none',
        'bonded_gpu',
        'bonded_cpu',
        'bonded_none',
    ]
    
    # Arguments to ignore (internal to CLI)
    ignore_args = {
        'command', 'scripts_path', 'skip_md', 'keep_temp', 'verbose'
    }
    
    # Process regular arguments
    for key, value in vars(args).items():
        if key in ignore_args or value is None:
            continue
        
        # Skip if not in mapping
        if key not in arg_mapping and key not in bool_flags:
            if getattr(args, 'verbose', False):
                print(f"Warning: Argument '{key}' not recognized, ignoring")
            continue
        
        # Handle boolean flags
        if key in bool_flags:
            if value:
                cmd_args.append(arg_mapping[key])
        else:
            cmd_key = arg_mapping[key]
            cmd_args.append(cmd_key)
            cmd_args.append(str(value))
    
    # Debug: print the command if verbose
    if getattr(args, 'verbose', False):
        print(f"\nRunning command:")
        print(' '.join(cmd_args))
        print()
    
    # Run the script
    try:
        result = subprocess.run(cmd_args, env=os.environ.copy())
        sys.exit(result.returncode)
    except FileNotFoundError as e:
        print(f"Error: Failed to execute {script_path}")
        print(f"Details: {e}")
        print("\nTroubleshooting:")
        print("1. Check if the script has execute permissions: chmod +x auto_mart_CG.sh")
        print("2. Check if the script has a valid shebang (#!/bin/bash)")
        print("3. Check if all dependencies (gmx, python3) are in PATH")
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nInterrupted by user")
        sys.exit(1)

def main():
    """Main CLI entry point.
    
    Parses command-line arguments and dispatches to the appropriate
    subcommand handler functions.
    
    Usage Examples:
        # Map AA trajectory to CG resolution
        auto_mart3 auto_map --index_cg beads.ndx --aa_tpr topol.tpr --aa_xtc traj.xtc
        
        # Generate CG topology
        auto_mart3 auto_gen_top --path_ff ./ff --itp_ligand cg.itp --name_molecule ligand --number_molecule 1 --output_topol topol.top
        
        # Analyze bonds, angles, dihedrals from trajectory
        auto_mart3 auto_analyze --bonds_ndx bonds.ndx --angles_ndx angles.ndx --dihedrals_ndx dihedrals.ndx --xtc_file traj.xtc --tpr_file topol.tpr --output_all_files ./analysis
        
        # Compute distribution statistics
        auto_mart3 auto_distributions --bonds_dir ./bonds --angles_dir ./angles --dihedrals_dir ./dihedrals --dir_to_output ./stats
        
        # Adapt ITP to match GRO atom names
        auto_mart3 auto_adapt_itp --input_itp input.itp --input_gro_ref reference.gro --output_itp_adapted output.itp
        
        # Setup CG system with solvent and ions
        auto_mart3 auto_prep --input_ref_dir ./input --output_dir ./output --input_gro molecule.gro --input_itp molecule.itp
        
        # Adjust potentials using Bayesian optimization
        auto_mart3 bayes_potential_adjust --bonds_ref_xvg_dir ./ref_bonds --angles_ref_xvg_dir ./ref_angles --dihedrals_ref_xvg_dir ./ref_dihedrals --bonds_sim_xvg_dir ./sim_bonds --angles_sim_xvg_dir ./sim_angles --dihedrals_sim_xvg_dir ./sim_dihedrals --itp_cg molecule.itp --ndx_bounds bonds.ndx --ndx_angles angles.ndx --ndx_dihedrals dihedrals.ndx
        
        # Plot distributions comparison
        auto_mart3 auto_plot_distributions --bonds_ref_dir ./ref_bonds --bonds_sim_dir ./sim_bonds --figures_dir ./figures
        
        # Run full AA to CG parametrization pipeline
        auto_mart3 auto_mart_aa --aa_tpr topol.tpr --aa_xtc traj.xtc --aa_gro structure.gro --aa_itp molecule.itp --beads_json beads.json --input_mdp minimization.mdp --path_ff ./ff --output_dir ./output --name_molecule ligand
        
        # Run CG parameter optimization pipeline with GPU acceleration
        auto_mart3 auto_mart_cg --INPUT_AUTO_MART_AA_DIR ./aa_results --OUTPUT_AUTO_MART_CG_DIR ./cg_optimization --cuda_visible_devices "0,1,2,3" --nb_gpu --pme_gpu --bonded_gpu --pin_on --npme 4
        
        # Run CG parameter optimization pipeline with CPU only
        auto_mart3 auto_mart_cg --INPUT_AUTO_MART_AA_DIR ./aa_results --OUTPUT_AUTO_MART_CG_DIR ./cg_optimization --nb_cpu --pme_cpu --bonded_cpu --pin_on --ntomp 20
    """
    parser = argparse.ArgumentParser(
        prog="auto_mart3",
        description="Auto_Mart3 - A pipeline to map atomistic (AA) trajectories to coarse-grained (CG) representations for Martini 3 force field parameters."
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Commands")
    
    # auto_map command
    map_parser = subparsers.add_parser("auto_map", 
                                   help="Map atomistic trajectory to coarse-grained (CG) resolution using GROMACS",
                                   description="""Map atomistic trajectory to coarse-grained (CG) resolution using GROMACS.""",
                                   formatter_class=argparse.RawDescriptionHelpFormatter)

    map_parser.add_argument("--index_cg", "-i", required=True,
                       help="CG index file (.ndx) containing bead definitions")
    map_parser.add_argument("--aa_tpr", "-t", required=True,
                       help="Atomistic topology file (.tpr) from AA simulation")
    map_parser.add_argument("--aa_xtc", "-x", required=True,
                       help="Atomistic trajectory file (.xtc) from AA simulation")
    map_parser.add_argument("--output_mapped", "-o", default="mapped.xtc",
                       help="Output CG trajectory filename (default: mapped.xtc)")
    map_parser.add_argument("--output_cg_gro", default="molecule.gro",
                       help="Output CG .gro file from first frame (default: molecule.gro)")
    map_parser.add_argument("--remove_pbc", action="store_true", default=True,
                       help="Remove PBC using 'gmx trjconv -pbc whole' (default: True)")
    map_parser.add_argument("--corrected_pbc", action="store_false", dest="remove_pbc",
                       help="Use when trajectory already has corrected PBC (skip trjconv step)")
    map_parser.add_argument("--verbose", "-v", action="store_true",
                       help="Enable verbose output with detailed progress messages")
    map_parser.add_argument("--dry_run", action="store_true",
                       help="Perform dry run (print commands without executing)")
    map_parser.add_argument("--keep_temp", action="store_true",
                       help="Keep temporary files for debugging")
    map_parser.add_argument("--gmx_cmd", default="gmx",
                       help="GROMACS command to use (default: gmx)")

    # auto_gen_top command
    top_parser = subparsers.add_parser("auto_gen_top", 
                                   help="Generate GROMACS topology (.top) file for CG simulations",
                                   description="""Generate a GROMACS topology (.top) file including selected force field.""",
                                   formatter_class=argparse.RawDescriptionHelpFormatter)

    top_parser.add_argument("--path_ff", type=str, required=True,
                       help="Directory containing force field ITP files")
    top_parser.add_argument("--itp_ligand", type=str, required=True,
                       help="Ligand ITP file to include (use 'None' or 'none' to skip)")
    top_parser.add_argument("--name_molecule", type=str, required=True,
                       help="Molecule name for the ligand")
    top_parser.add_argument("--number_molecule", type=int, required=True,
                       help="Number of molecules of this type in the system")
    top_parser.add_argument("--output_topol", type=str, required=True,
                       help="Output topology filename")
    top_parser.add_argument("--ff", type=str, default="martini_v3.0.0.itp",
                       help="Force field ITP filename (default: martini_v3.0.0.itp)")
    top_parser.add_argument("--ions", type=str, default="martini_v3.0.0_ions_v1.itp",
                       help="Ions ITP filename (default: martini_v3.0.0_ions_v1.itp)")
    top_parser.add_argument("--solvent", type=str, default="martini_v3.0.0_solvents_v1.itp",
                       help="Solvent ITP filename (default: martini_v3.0.0_solvents_v1.itp)")
    top_parser.add_argument("--title_comments", type=str, default="",
                       help="Optional comment header at the top of the topology file")
    top_parser.add_argument("--title_system", type=str, default="",
                       help="Optional system name (used in [ system ] section)")    
    
    # auto_analyze command
    analyze_parser = subparsers.add_parser("auto_analyze", 
                                       help="Calculate bonds, angles, and dihedrals from MD trajectories",
                                       description="""Compute bond distances, angles, and dihedral distributions from MD trajectories.""",
                                       formatter_class=argparse.RawDescriptionHelpFormatter)

    analyze_parser.add_argument("--bonds_ndx", required=True,
                           help="Index file defining bonds")
    analyze_parser.add_argument("--angles_ndx", required=True,
                           help="Index file defining angles")
    analyze_parser.add_argument("--dihedrals_ndx", required=True,
                           help="Index file defining dihedrals")
    analyze_parser.add_argument("--xtc_file", required=True,
                           help="Input trajectory file (.xtc)")
    analyze_parser.add_argument("--tpr_file", required=True,
                           help="Input topology file (.tpr)")
    analyze_parser.add_argument("--output_all_files", required=True,
                           help="Base output directory path")
    analyze_parser.add_argument("--index", default=None,
                           help="Index file for selecting groups")
    analyze_parser.add_argument("--remove_pbc", action="store_true", default=False,
                           help="Remove PBC and align trajectory")
    analyze_parser.add_argument("--group_1", default=None,
                           help="Group for fitting/alignment")
    analyze_parser.add_argument("--group_2", default="System",
                           help="Group for output/trajectory writing")
    analyze_parser.add_argument("--keep_intermediate", action="store_true", default=False,
                           help="Keep all intermediate files for debugging")

    # auto_distributions command
    dist_parser = subparsers.add_parser("auto_distributions", 
                                    help="Calculate statistical properties from XVG distribution files",
                                    description="""Compute statistical properties from XVG distribution files.""",
                                    formatter_class=argparse.RawDescriptionHelpFormatter)

    dist_parser.add_argument("--bonds_dir", default="bonds",
                        help="Directory containing bond XVG files (default: bonds/)")
    dist_parser.add_argument("--angles_dir", default="angles",
                        help="Directory containing angle XVG files (default: angles/)")
    dist_parser.add_argument("--dihedrals_dir", default="dihedrals",
                        help="Directory containing dihedral XVG files (default: dihedrals/)")
    dist_parser.add_argument("--dir_to_output", default="TSV_statistics",
                        help="Output directory for TSV statistics files (default: TSV_statistics/)")
    dist_parser.add_argument("--bond_out", default="bond_statistics.tsv",
                        help="Output filename for bond statistics (default: bond_statistics.tsv)")
    dist_parser.add_argument("--angle_out", default="angle_statistics.tsv",
                        help="Output filename for angle statistics (default: angle_statistics.tsv)")
    dist_parser.add_argument("--dihedral_out", default="dihedral_statistics.tsv",
                        help="Output filename for dihedral statistics (default: dihedral_statistics.tsv)")
    
    # auto_adapt_itp command
    adapt_parser = subparsers.add_parser("auto_adapt_itp", 
                                     help="Update atom names in ITP file to match GRO reference",
                                     description="""Update atom names in a GROMACS ITP file to match those from a reference GRO file.""",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    adapt_parser.add_argument("--input_itp", required=True,
                         help="Input ITP file to be modified")
    adapt_parser.add_argument("--input_gro_ref", required=True,
                         help="Reference GRO file containing correct atom names")
    adapt_parser.add_argument("--output_itp_adapted", required=True,
                         help="Output ITP file with updated atom names")
    adapt_parser.add_argument("--verbose", action="store_true",
                         help="Print detailed debug information")    

    # auto_prep command
    cg_setup_parser = subparsers.add_parser("auto_prep", 
                                        help="Setup CG system with solvent and ions using Martini 3 force field",
                                        description="""Prepare a coarse-grained (CG) molecular dynamics system.""",
                                        formatter_class=argparse.RawDescriptionHelpFormatter)

    cg_setup_parser.add_argument("--input_ref_dir", required=True,
                            help="Directory containing input files")
    cg_setup_parser.add_argument("--output_dir", required=True,
                            help="Output directory for final system files")
    cg_setup_parser.add_argument("--input_gro", default="cg.gro",
                            help="Input CG structure filename (default: cg.gro)")
    cg_setup_parser.add_argument("--input_itp", default="cg.itp",
                            help="Input CG topology ITP filename (default: cg.itp)")
    cg_setup_parser.add_argument("--input_topol", default="topol_cg.top",
                            help="Input system topology filename (default: topol_cg.top)")
    cg_setup_parser.add_argument("--input_ff_dir", default=None,
                            help="Directory containing force field ITP files")
    cg_setup_parser.add_argument("--ff", default="martini_v3.0.0.itp",
                            help="Force field ITP filename (default: martini_v3.0.0.itp)")
    cg_setup_parser.add_argument("--ions", default="martini_v3.0.0_ions_v1.itp",
                            help="Ions ITP filename (default: martini_v3.0.0_ions_v1.itp)")
    cg_setup_parser.add_argument("--solvent", default="martini_v3.0.0_solvents_v1.itp",
                            help="Solvent ITP filename (default: martini_v3.0.0_solvents_v1.itp)")
    cg_setup_parser.add_argument("--water_dir", default=None,
                            help="Directory containing water bead GRO file")
    cg_setup_parser.add_argument("--water_file_gro", default="water.gro",
                            help="Water structure filename (default: water.gro)")
    cg_setup_parser.add_argument("--input_mdp_dir", default=None,
                            help="MDP files directory")
    cg_setup_parser.add_argument("--input_name_file_mdp", default="minimization.mdp",
                            help="MDP filename for minimization (default: minimization.mdp)")
    cg_setup_parser.add_argument("--ions_mdp_dir", default=None,
                            help="Directory for ion addition MDP files")
    cg_setup_parser.add_argument("--ions_file_mdp", default="ions.mdp",
                            help="MDP filename for ion addition (default: ions.mdp)")
    cg_setup_parser.add_argument("--use_distance_from_atom", action="store_true",
                            help="Use distance from atom to calculate box size")
    cg_setup_parser.add_argument("--distance_from_atom", type=float, default=2.0,
                            help="Distance (nm) from molecule to box edge (default: 2.0)")
    cg_setup_parser.add_argument("--box_size", type=float, default=None,
                            help="Fixed cubic box size in nm")
    cg_setup_parser.add_argument("--max_solvent", type=int, default=None,
                            help="Maximum number of solvent molecules to add")
    cg_setup_parser.add_argument("--solvent_radius", type=float, default=0.21,
                            help="Radius of solvent beads in nm (default: 0.21)")
    cg_setup_parser.add_argument("--salt", type=float, default=0.15,
                            help="Salt concentration in Molar (default: 0.15 M)")
    cg_setup_parser.add_argument("--skip_grompp", action="store_true",
                            help="Skip grompp step (do not generate em.tpr)")
    cg_setup_parser.add_argument("--keep_temp", action="store_true",
                            help="Keep temporary files for debugging")
    cg_setup_parser.add_argument("--verbose", action="store_true",
                            help="Verbose output with detailed progress information")

    # bayes_potential_adjust command
    potential_parser = subparsers.add_parser("bayes_potential_adjust", 
                                 help="Adjust CG force constants using Bayesian update with R² correction",
                                 description="""Automated optimization script for Martini coarse-grained force field parameters.""",
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

    required_group = potential_parser.add_argument_group("Required arguments")
    required_group.add_argument("--bonds_ref_xvg_dir", required=True,
                         help="Directory containing reference bond XVG files")
    required_group.add_argument("--angles_ref_xvg_dir", required=True,
                         help="Directory containing reference angle XVG files")
    required_group.add_argument("--dihedrals_ref_xvg_dir", required=True,
                         help="Directory containing reference dihedral XVG files")
    required_group.add_argument("--bonds_sim_xvg_dir", required=True,
                         help="Directory containing simulated bond XVG files")
    required_group.add_argument("--angles_sim_xvg_dir", required=True,
                         help="Directory containing simulated angle XVG files")
    required_group.add_argument("--dihedrals_sim_xvg_dir", required=True,
                         help="Directory containing simulated dihedral XVG files")
    required_group.add_argument("--itp_cg", required=True,
                         help="CG ITP file with topology information")
    required_group.add_argument("--ndx_bounds", required=True,
                         help="NDX file with bond indices")
    required_group.add_argument("--ndx_angles", required=True,
                         help="NDX file with angle indices")
    required_group.add_argument("--ndx_dihedrals", required=True,
                         help="NDX file with dihedral indices")

    optional_group = potential_parser.add_argument_group("Optional arguments")
    optional_group.add_argument("--molecule_name", default="molecule",
                         help="Molecule name for ITP header (default: molecule)")
    optional_group.add_argument("--distribution_points", type=int, default=200,
                         help="Number of points for distribution interpolation (default: 200)")
    optional_group.add_argument("--T0", type=float, default=10.0,
                         help="Initial temperature for simulated annealing (default: 10.0)")
    optional_group.add_argument("--alpha", type=float, default=0.85,
                         help="Cooling factor for simulated annealing (default: 0.85)")
    optional_group.add_argument("--itp_out", default="cg.itp",
                         help="Output ITP filename (default: cg.itp)")

    prefix_group = potential_parser.add_argument_group("XVG prefix options")
    prefix_group.add_argument("--prefix_xvg_bond_ref", default="bond_",
                         help="Prefix for reference bond XVG files (default: bond_)")
    prefix_group.add_argument("--prefix_xvg_angle_ref", default="ang_",
                         help="Prefix for reference angle XVG files (default: ang_)")
    prefix_group.add_argument("--prefix_xvg_dihedral_ref", default="dih_",
                         help="Prefix for reference dihedral XVG files (default: dih_)")
    prefix_group.add_argument("--prefix_xvg_bond_sim", default="bond_",
                         help="Prefix for simulated bond XVG files (default: bond_)")
    prefix_group.add_argument("--prefix_xvg_angle_sim", default="ang_",
                         help="Prefix for simulated angle XVG files (default: ang_)")
    prefix_group.add_argument("--prefix_xvg_dihedral_sim", default="dih_",
                         help="Prefix for simulated dihedral XVG files (default: dih_)")

    force_limits_group = potential_parser.add_argument_group("Force limits options")
    force_limits_group.add_argument("--min_force_bond", type=float, default=None,
                             help="Minimum force for bonds (default: 750.0)")
    force_limits_group.add_argument("--min_force_angle", type=float, default=None,
                             help="Minimum force for angles (default: 15.0)")
    force_limits_group.add_argument("--min_force_dihedral", type=float, default=None,
                             help="Minimum force for dihedrals (default: 15.0)")
    force_limits_group.add_argument("--max_force_bond", type=float, default=None,
                             help="Maximum force for bonds (default: 10000.0)")
    force_limits_group.add_argument("--max_force_angle", type=float, default=None,
                             help="Maximum force for angles (default: 150.0)")
    force_limits_group.add_argument("--max_force_dihedral", type=float, default=None,
                             help="Maximum force for dihedrals (default: 150.0)")
    force_limits_group.add_argument("--default_force_bond", type=float, default=None,
                             help="Default force for bonds (default: 1250.0)")
    force_limits_group.add_argument("--default_force_angle", type=float, default=None,
                             help="Default force for angles (default: 25.0)")
    force_limits_group.add_argument("--default_force_dihedral", type=float, default=None,
                             help="Default force for dihedrals (default: 25.0)")

    control_group = potential_parser.add_argument_group("Control options")
    control_group.add_argument("--verbose", action="store_true",
                         help="Verbose output")    
    
    # auto_plot_distributions command
    plot_parser = subparsers.add_parser("auto_plot_distributions", 
                                    help="Plot reference vs simulated distributions from XVG files",
                                    description="""Plot reference and simulated distributions for bonds, angles, and dihedrals.""",
                                    formatter_class=argparse.RawDescriptionHelpFormatter)

    plot_parser.add_argument("--bonds_ref_dir", default="bonds",
                        help="Directory containing reference bond XVG files (default: bonds)")
    plot_parser.add_argument("--angles_ref_dir", default="angles",
                        help="Directory containing reference angle XVG files (default: angles)")
    plot_parser.add_argument("--dihedrals_ref_dir", default="dihedrals",
                        help="Directory containing reference dihedral XVG files (default: dihedrals)")
    plot_parser.add_argument("--bonds_sim_dir", default="bonds",
                        help="Directory containing simulated bond XVG files (default: bonds)")
    plot_parser.add_argument("--angles_sim_dir", default="angles",
                        help="Directory containing simulated angle XVG files (default: angles)")
    plot_parser.add_argument("--dihedrals_sim_dir", default="dihedrals",
                        help="Directory containing simulated dihedral XVG files (default: dihedrals)")
    plot_parser.add_argument("--figures_dir", default="figures",
                        help="Output directory for figures (default: figures)")
    
    # auto_mart_aa command
    pipeline_aa_parser = subparsers.add_parser("auto_mart_aa", 
        help="Run full AA to CG parametrization pipeline",
        description="""Run the complete AA to CG parametrization pipeline.""")

    essential_group_aa = pipeline_aa_parser.add_argument_group("Required arguments")
    essential_group_aa.add_argument("--aa_tpr", required=True, help="AA TPR file")
    essential_group_aa.add_argument("--aa_xtc", required=True, help="AA trajectory XTC file")
    essential_group_aa.add_argument("--aa_gro", required=True, help="AA GRO file")
    essential_group_aa.add_argument("--aa_itp", required=True, help="AA ITP file")
    essential_group_aa.add_argument("--beads_json", required=True, help="Beads definition JSON file")
    essential_group_aa.add_argument(
        "--input_mdp",
        default=get_default_mdp(),
        help="MDP file for grompp (default: package data)"
    )

    essential_group_aa.add_argument(
        "--path_ff",
        default=get_default_ff_dir(),
        help="Force field directory (default: package data)"
    )
    
    essential_group_aa.add_argument("--output_dir", required=True, help="Output directory")
    essential_group_aa.add_argument("--name_molecule", required=True, help="Molecule name")

    optional_group_aa = pipeline_aa_parser.add_argument_group("Optional arguments")
    optional_group_aa.add_argument("--force_application", default="random=[1250,30;30,1]",
                              help="Force application method")
    optional_group_aa.add_argument("--beads_position", default="geom",
                              choices=["com", "geom", "min"],
                              help="Beads positioning method (default: 'geom')")
    optional_group_aa.add_argument("--cycle_restr", default="fix=3,mode=cycle",
                              help='Cycle restraints configuration')
    optional_group_aa.add_argument("--maxwarn", type=int, default=2,
                              help="Max warnings for gmx grompp (default: 2)")
    optional_group_aa.add_argument("--distance_from_atom", type=float, default=2.0,
                              help="Distance for solvent placement in nm (default: 2.0)")
    optional_group_aa.add_argument("--salt", type=float, default=0.15,
                              help="Salt concentration in M (default: 0.15)")

    control_group_aa = pipeline_aa_parser.add_argument_group("Control options")
    control_group_aa.add_argument("--verbose", action="store_true",
                             help="Verbose output")
    control_group_aa.add_argument("--keep_temp", action="store_true",
                             help="Keep temporary files (for debugging)")
    
    dist_group_aa = pipeline_aa_parser.add_argument_group("Distribution statistics options")
    dist_group_aa.add_argument("--run_distributions", action="store_true", 
                           help="Run distribution statistics after analysis")
    dist_group_aa.add_argument("--dist_bonds_dir", default="bonds",
                           help="Directory for bond XVG files (default: bonds)")
    dist_group_aa.add_argument("--dist_angles_dir", default="angles",
                           help="Directory for angle XVG files (default: angles)")
    dist_group_aa.add_argument("--dist_dihedrals_dir", default="dihedrals",
                           help="Directory for dihedral XVG files (default: dihedrals)")
    dist_group_aa.add_argument("--dist_output_dir", default="STATISTICS",
                           help="Output directory for statistics TSV files (default: STATISTICS)")
    dist_group_aa.add_argument("--dist_bond_out", default="bond_statistics.tsv",
                           help="Bond statistics output filename (default: bond_statistics.tsv)")
    dist_group_aa.add_argument("--dist_angle_out", default="angle_statistics.tsv",
                           help="Angle statistics output filename (default: angle_statistics.tsv)")
    dist_group_aa.add_argument("--dist_dihedral_out", default="dihedral_statistics.tsv",
                           help="Dihedral statistics output filename (default: dihedral_statistics.tsv)")

    # auto_mart_cg command 
    pipeline_cg_parser = subparsers.add_parser("auto_mart_cg", 
                                               help="Run CG parameter optimization pipeline (iterative Bayesian optimization)",
                                               description="""Automated parameter optimization for Coarse-Grained (CG) molecular dynamics 
                                               simulations using Bayesian potentials. Performs iterative optimization of bonds, angles, 
                                               and dihedrals force constants by comparing reference distributions from AA simulations 
                                               with CG simulations.""",
                                               formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # Required arguments
    essential_group_cg = pipeline_cg_parser.add_argument_group("Required arguments")
    essential_group_cg.add_argument("--INPUT_AUTO_MART_AA_DIR", required=True, dest="input_auto_mart_aa_dir",
                                   help="Directory containing folders from auto_mart_aa results (BONDS_ANGLES_DIHEDRALS_XVG_REF, CG_MARTINI3, GMX, MDRUN_CG, STATISTICS)")
    essential_group_cg.add_argument("--OUTPUT_AUTO_MART_CG_DIR", required=True, dest="output_auto_mart_cg_dir",
                                   help="Directory for CG optimization output")
    
    # Simulation Parameters
    sim_group = pipeline_cg_parser.add_argument_group("Simulation Parameters")
    sim_group.add_argument("--ntomp", type=int, default=10, dest="ntomp",
                          help="OpenMP threads (default: 10). Set to 'None' for auto-detection.")
    sim_group.add_argument("--ntmpi", type=int, default=1, dest="ntmpi",
                          help="MPI threads (default: 1). Set to 'None' for auto-detection.")
    sim_group.add_argument("--ref_t", type=float, default=310, dest="ref_t",
                          help="Reference temperature in K (default: 310)")
    sim_group.add_argument("--ref_p", type=float, default=1.0, dest="ref_p",
                          help="Reference pressure in bar (default: 1.0)")
    sim_group.add_argument("--dt_nvt_ps", type=float, default=0.001, dest="dt_nvt_ps",
                          help="NVT timestep in ps (default: 0.001)")
    sim_group.add_argument("--time_nvt_ps", type=float, default=1000, dest="time_nvt_ps",
                          help="NVT duration in ps (default: 1000)")
    sim_group.add_argument("--dt_npt_ps", type=float, default=0.001, dest="dt_npt_ps",
                          help="NPT timestep in ps (default: 0.001)")
    sim_group.add_argument("--time_npt_ps", type=float, default=5000, dest="time_npt_ps",
                          help="NPT duration in ps (default: 5000)")
    sim_group.add_argument("--dt_md_ps", type=float, default=0.002, dest="dt_md_ps",
                          help="MD production timestep in ps (default: 0.002)")
    sim_group.add_argument("--time_md_ps", type=float, default=10000, dest="time_md_ps",
                          help="MD production duration in ps (default: 10000)")
    
    # File Names
    files_group = pipeline_cg_parser.add_argument_group("File Names")
    files_group.add_argument("--topol_cg_file", default="topol_cg.top", dest="topol_cg_file",
                            help="CG topology file name (default: topol_cg.top)")
    files_group.add_argument("--solv_ions_gro", default="solv_ions_CG.gro", dest="solv_ions_gro",
                            help="Solvated ions GRO file (default: solv_ions_CG.gro)")
    files_group.add_argument("--em_tpr", default="em.tpr", dest="em_tpr",
                            help="Energy minimization TPR file (default: em.tpr)")
    files_group.add_argument("--itp_to_optimize", default="cg.itp", dest="itp_to_optimize",
                            help="ITP file to optimize (default: cg.itp)")
    
    # Force Constants Limits
    force_group = pipeline_cg_parser.add_argument_group("Force Constants Limits")
    force_group.add_argument("--min_force_bond", type=float, default=500.0, dest="min_force_bond",
                            help="Minimum bond force constant (default: 500.0)")
    force_group.add_argument("--max_force_bond", type=float, default=50000.0, dest="max_force_bond",
                            help="Maximum bond force constant (default: 50000.0)")
    force_group.add_argument("--default_force_bond", type=float, default=1250.0, dest="default_force_bond",
                            help="Default bond force constant (default: 1250.0)")
    force_group.add_argument("--min_force_angle", type=float, default=10.0, dest="min_force_angle",
                            help="Minimum angle force constant (default: 10.0)")
    force_group.add_argument("--max_force_angle", type=float, default=150.0, dest="max_force_angle",
                            help="Maximum angle force constant (default: 150.0)")
    force_group.add_argument("--default_force_angle", type=float, default=25.0, dest="default_force_angle",
                            help="Default angle force constant (default: 25.0)")
    force_group.add_argument("--min_force_dihedral", type=float, default=10.0, dest="min_force_dihedral",
                            help="Minimum dihedral force constant (default: 10.0)")
    force_group.add_argument("--max_force_dihedral", type=float, default=150.0, dest="max_force_dihedral",
                            help="Maximum dihedral force constant (default: 150.0)")
    force_group.add_argument("--default_force_dihedral", type=float, default=25.0, dest="default_force_dihedral",
                            help="Default dihedral force constant (default: 25.0)")
    
    # Optimization Parameters
    opt_group = pipeline_cg_parser.add_argument_group("Optimization Parameters")
    opt_group.add_argument("--T0", type=float, default=10.0, dest="T0",
                          help="Initial simulated annealing temperature (default: 10.0)")
    opt_group.add_argument("--alpha", type=float, default=0.85, dest="alpha",
                          help="Cooling factor (default: 0.85)")
    opt_group.add_argument("--distribution_points", type=int, default=100, dest="distribution_points",
                          help="Points for R² calculation (default: 100)")
    opt_group.add_argument("--n_iter", type=int, default=30, dest="n_iter",
                          help="Number of optimization iterations (default: 30)")
    
    # Input/Output Paths
    paths_group = pipeline_cg_parser.add_argument_group("Input/Output Paths")
    paths_group.add_argument("--bonds_ndx", default=None, dest="bonds_ndx",
                            help="Bonds index file (default: INPUT_AUTO_MART_AA_DIR/CG_MARTINI3/NDX/bonds.ndx)")
    paths_group.add_argument("--angles_ndx", default=None, dest="angles_ndx",
                            help="Angles index file (default: INPUT_AUTO_MART_AA_DIR/CG_MARTINI3/NDX/angles.ndx)")
    paths_group.add_argument("--dihedrals_ndx", default=None, dest="dihedrals_ndx",
                            help="Dihedrals index file (default: INPUT_AUTO_MART_AA_DIR/CG_MARTINI3/NDX/dihedrals.ndx)")
    paths_group.add_argument("--bonds_ref_xvg_dir", default=None, dest="bonds_ref_xvg_dir",
                            help="Bonds reference XVG directory (default: INPUT_AUTO_MART_AA_DIR/BONDS_ANGLES_DIHEDRALS_XVG_REF/bonds)")
    paths_group.add_argument("--angles_ref_xvg_dir", default=None, dest="angles_ref_xvg_dir",
                            help="Angles reference XVG directory (default: INPUT_AUTO_MART_AA_DIR/BONDS_ANGLES_DIHEDRALS_XVG_REF/angles)")
    paths_group.add_argument("--dihedrals_ref_xvg_dir", default=None, dest="dihedrals_ref_xvg_dir",
                            help="Dihedrals reference XVG directory (default: INPUT_AUTO_MART_AA_DIR/BONDS_ANGLES_DIHEDRALS_XVG_REF/dihedrals)")
    
    # Prefixes for XVG Files
    prefix_group_cg = pipeline_cg_parser.add_argument_group("Prefixes for XVG Files")
    prefix_group_cg.add_argument("--prefix_xvg_bond_ref", default="bond_", dest="prefix_xvg_bond_ref",
                                help="Prefix for bond reference XVG files (default: bond_)")
    prefix_group_cg.add_argument("--prefix_xvg_angle_ref", default="ang_", dest="prefix_xvg_angle_ref",
                                help="Prefix for angle reference XVG files (default: ang_)")
    prefix_group_cg.add_argument("--prefix_xvg_dihedral_ref", default="dih_", dest="prefix_xvg_dihedral_ref",
                                help="Prefix for dihedral reference XVG files (default: dih_)")
    prefix_group_cg.add_argument("--prefix_xvg_bond_sim", default="bond_", dest="prefix_xvg_bond_sim",
                                help="Prefix for bond simulation XVG files (default: bond_)")
    prefix_group_cg.add_argument("--prefix_xvg_angle_sim", default="ang_", dest="prefix_xvg_angle_sim",
                                help="Prefix for angle simulation XVG files (default: ang_)")
    prefix_group_cg.add_argument("--prefix_xvg_dihedral_sim", default="dih_", dest="prefix_xvg_dihedral_sim",
                                help="Prefix for dihedral simulation XVG files (default: dih_)")
    
    # GROMACS Processing
    gmx_group = pipeline_cg_parser.add_argument_group("GROMACS Processing")
    gmx_group.add_argument("--group_2", default="System", dest="group_2",
                          help="Output group for trajectory processing (default: System)")
    gmx_group.add_argument("--group_1", default="System", dest="group_1",
                          help="Center group for trajectory processing (default: System)")
    gmx_group.add_argument("--index", default="None", dest="index",
                          help="Index file for GROMACS selections (default: None)")
    
    # GROMACS MDRUN GPU/CPU Configuration
    mdrun_group = pipeline_cg_parser.add_argument_group("GROMACS MDRUN Configuration (GPU/CPU)")
    mdrun_group.add_argument("--pin_on", action="store_true", dest="pin_on",
                            help="Enable thread pinning for better performance (default: on)")
    mdrun_group.add_argument("--pin_off", action="store_true", dest="pin_off",
                            help="Disable thread pinning")
    mdrun_group.add_argument("--nb_gpu", action="store_true", dest="nb_gpu",
                            help="Use GPU for non-bonded interactions")
    mdrun_group.add_argument("--nb_cpu", action="store_true", dest="nb_cpu",
                            help="Use CPU for non-bonded interactions")
    mdrun_group.add_argument("--nb_none", action="store_true", dest="nb_none",
                            help="Auto-detect for non-bonded interactions (default)")
    mdrun_group.add_argument("--pme_gpu", action="store_true", dest="pme_gpu",
                            help="Use GPU for PME calculations")
    mdrun_group.add_argument("--pme_cpu", action="store_true", dest="pme_cpu",
                            help="Use CPU for PME calculations")
    mdrun_group.add_argument("--pme_none", action="store_true", dest="pme_none",
                            help="Auto-detect for PME calculations (default)")
    mdrun_group.add_argument("--bonded_gpu", action="store_true", dest="bonded_gpu",
                            help="Use GPU for bonded interactions")
    mdrun_group.add_argument("--bonded_cpu", action="store_true", dest="bonded_cpu",
                            help="Use CPU for bonded interactions")
    mdrun_group.add_argument("--bonded_none", action="store_true", dest="bonded_none",
                            help="Auto-detect for bonded interactions (default)")
    mdrun_group.add_argument("--npme", type=str, default=None, dest="npme",
                            help="Number of PME threads. Set to 'None' for auto-detection (default: auto)")
    mdrun_group.add_argument("--cuda_visible_devices", type=str, default=None, dest="cuda_visible_devices",
                            help="CUDA visible devices (e.g., '0,1,2,3')")
    mdrun_group.add_argument("--gpu_id", type=str, default=None, dest="gpu_id",
                            help="Alias for CUDA_VISIBLE_DEVICES")
    
    # Other
    other_group = pipeline_cg_parser.add_argument_group("Other")
    other_group.add_argument("--molecule_name", default="molecule", dest="molecule_name",
                            help="Molecule name for potential adjustment (default: molecule)")
    
    # Boolean flags
    bool_group = pipeline_cg_parser.add_argument_group("Boolean flags")
    bool_group.add_argument("--no_remove_gmx_files_in_iter", action="store_true", dest="no_remove_gmx_files_in_iter",
                           help="Preserve GMX temporary files (*.xtc, *.edr, etc.)")
    
    # Control options (internal to CLI)
    control_group_cg = pipeline_cg_parser.add_argument_group("Control options (internal)")
    control_group_cg.add_argument("--scripts_path", default=None, dest="scripts_path",
                                 help="Path to scripts directory (auto-detected if not provided)")
    control_group_cg.add_argument("--verbose", action="store_true", dest="verbose",
                                 help="Verbose output")
    
    # Show version
    parser.add_argument("--version", action="version", version="Auto_Mart3 1.0.0")
    
    args = parser.parse_args()
    if args.command == "auto_mart_aa":
        if args.input_mdp is None:
            print("Error: --input_mdp not provided and no default found in package.")
            sys.exit(1)

        if args.path_ff is None:
            print("Error: --path_ff not provided and no default found in package.")
            sys.exit(1)
    
    # Handle mutually exclusive GPU/CPU options for auto_mart_cg
    if args.command == "auto_mart_cg":
        # Handle pinning options (default to on if neither specified)
        if not args.pin_on and not args.pin_off:
            args.pin_on = True
        
        # Handle nb options (default to none if none specified)
        if not args.nb_gpu and not args.nb_cpu and not args.nb_none:
            args.nb_none = True
        
        # Handle pme options (default to none if none specified)
        if not args.pme_gpu and not args.pme_cpu and not args.pme_none:
            args.pme_none = True
        
        # Handle bonded options (default to none if none specified)
        if not args.bonded_gpu and not args.bonded_cpu and not args.bonded_none:
            args.bonded_none = True
        
        # Handle CUDA visible devices from gpu_id alias
        if args.gpu_id and not args.cuda_visible_devices:
            args.cuda_visible_devices = args.gpu_id
        
        # Set default values for paths
        if args.bonds_ndx is None:
            args.bonds_ndx = f"{args.input_auto_mart_aa_dir}/CG_MARTINI3/NDX/bonds.ndx"
        if args.angles_ndx is None:
            args.angles_ndx = f"{args.input_auto_mart_aa_dir}/CG_MARTINI3/NDX/angles.ndx"
        if args.dihedrals_ndx is None:
            args.dihedrals_ndx = f"{args.input_auto_mart_aa_dir}/CG_MARTINI3/NDX/dihedrals.ndx"
        if args.bonds_ref_xvg_dir is None:
            args.bonds_ref_xvg_dir = f"{args.input_auto_mart_aa_dir}/BONDS_ANGLES_DIHEDRALS_XVG_REF/bonds"
        if args.angles_ref_xvg_dir is None:
            args.angles_ref_xvg_dir = f"{args.input_auto_mart_aa_dir}/BONDS_ANGLES_DIHEDRALS_XVG_REF/angles"
        if args.dihedrals_ref_xvg_dir is None:
            args.dihedrals_ref_xvg_dir = f"{args.input_auto_mart_aa_dir}/BONDS_ANGLES_DIHEDRALS_XVG_REF/dihedrals"
    
    if args.command == "auto_map":
        run_mapping_script(args)
    elif args.command == "auto_gen_top":
        run_topology_script(args)
    elif args.command == "auto_analyze":
        run_analysis_script(args)
    elif args.command == "auto_distributions":
        run_distribution_script(args)
    elif args.command == "auto_adapt_itp":
        run_adaptation_script(args)
    elif args.command == "bayes_potential_adjust":
        run_potential_adjustment_script(args)
    elif args.command == "auto_prep":
        run_cg_setup_script(args)
    elif args.command == "auto_plot_distributions":
        run_plot_distributions_script(args) 
    elif args.command == "auto_mart_aa":
        run_shell_script_AA(args)
    elif args.command == "auto_mart_cg":
        run_shell_script_CG(args)
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()
