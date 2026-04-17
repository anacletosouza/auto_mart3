#!/usr/bin/env python3
"""
Command-line interface for Bayesian Potentials package.
"""

import argparse
import sys
import os
import subprocess
from pathlib import Path


def get_package_dir():
    """Get the package installation directory."""
    return Path(__file__).parent.resolve()

def get_scripts_dir():
    """Get the scripts directory inside the package."""
    return get_package_dir() / "scripts"

def get_bin_dir():
    """Get the bin directory inside the package."""
    return get_package_dir() / "bin"

def get_data_dir():
    """Get the data directory (relative to package or installed)."""
    pkg_dir = get_package_dir()
    
    # Try relative to package
    data_dir = pkg_dir.parent / "data"
    if data_dir.exists():
        return data_dir
    
    # Try installed location
    import site
    for site_dir in site.getsitepackages():
        data_dir = Path(site_dir) / "bayesian_potentials_data"
        if data_dir.exists():
            return data_dir
    
    return None

def get_default_mdp():
    """Get the default MDP file from the package data directory."""
    data_dir = get_data_dir()
    if data_dir:
        mdp_file = data_dir / "mdp" / "minimization.mdp"
        if mdp_file.exists():
            return str(mdp_file)
    
    # Try to find in common locations
    possible_paths = [
        Path.cwd() / "data" / "mdp" / "minimization.mdp",
        Path(__file__).parent.parent / "data" / "mdp" / "minimization.mdp",
    ]
    
    for path in possible_paths:
        if path.exists():
            return str(path)
    
    return None

def get_default_ff_dir():
    """Get default force field directory from data directory."""
    data_dir = get_data_dir()
    if data_dir:
        ff_dir = data_dir / "ff_files"
        if ff_dir.exists():
            return str(ff_dir)
    return None

def find_script(script_name):
    """Find a script in multiple possible locations."""
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
    """Run the mapping script with proper module imports."""
    from bayesian_potentials.scripts.map_aa_to_cg import main as map_main
    
    sys_argv = ["map_aa_to_cg.py"]
    
    # Map arguments for the mapping script
    for key, value in vars(args).items():
        if key == "command" or value is None:
            continue
        
        # Handle the special case: cg_ndx should become index_cg for mapping script
        if key == "cg_ndx":
            sys_argv.append("--index_cg")
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
    """Run the topology generation script."""
    from bayesian_potentials.scripts.generate_cg_top import main as top_main
    
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
    """Run the bonds/angles/dihedrals analysis script."""
    from bayesian_potentials.scripts.generate_bonds_angles_dihedrals import main as analysis_main
    
    sys_argv = ["generate_bonds_angles_dihedrals.py"]
    
    for key, value in vars(args).items():
        if key == "command" or value is None:
            continue
        
        # Skip internal flags
        if key in ["skip_grompp", "verbose"]:
            continue
        
        # Generate_bonds_angles_dihedrals.py expects underscores, not hyphens
        cmd_key = "--" + key  # Keep underscores
        
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
    """Run the bond/angle/dihedral distribution statistics script."""
    from bayesian_potentials.scripts.bp_distributions import main as bp_main
    
    sys_argv = ["bp_distributions.py"]
    
    for key, value in vars(args).items():
        if key == "command" or value is None:
            continue
        
        # Skip internal flags
        if key in ["skip_grompp", "verbose"]:
            continue
        
        # bp_distributions.py expects underscores, not hyphens
        cmd_key = "--" + key  # Keep underscores
        
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
    """Run the ITP adaptation script to match GRO atom names."""
    from bayesian_potentials.scripts.adaptation_gro_itp import main as adapt_main
    
    sys_argv = ["adaptation_gro_itp.py"]
    
    for key, value in vars(args).items():
        if key == "command" or value is None:
            continue
        
        # Skip internal flags
        if key in ["verbose"]:
            continue
        
        cmd_key = "--" + key.replace("_", "_")  # Keep underscores
        
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
    """Run the CG system setup pipeline script."""
    from bayesian_potentials.scripts.bp_prep import main as cg_setup_main
    
    sys_argv = ["bp_prep.py"]
    
    for key, value in vars(args).items():
        if key == "command" or value is None:
            continue
        
        # Skip internal flags
        if key in ["verbose"]:
            continue
        
        # Convert underscores to hyphens for command line arguments
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

def run_force_adjustment_script(args):
    """Run the force constant adjustment script using Bayesian update and simulated annealing."""
    from bayesian_potentials.scripts.force_adjustment import main as force_main
    
    sys_argv = ["force_adjustment.py"]
    
    for key, value in vars(args).items():
        if key == "command" or value is None:
            continue
        
        # Skip internal flags
        if key in ["verbose"]:
            continue
        
        # Build command line arguments
        cmd_key = "--" + key
        
        if isinstance(value, bool):
            if value:
                sys_argv.append(cmd_key)
            # If False, don't add the flag (it's handled by the complementary flag)
        else:
            sys_argv.append(cmd_key)
            sys_argv.append(str(value))
    
    original_argv = sys.argv
    sys.argv = sys_argv
    
    try:
        force_main()
    finally:
        sys.argv = original_argv

def run_shell_script_AA(args):
    """Run the main shell pipeline script for AA to CG."""
    script_path = find_script("autoparam_AA.sh")
    
    if not script_path:
        print(f"Error: autoparam_AA.sh not found")
        print("Please ensure the script is in one of these locations:")
        print("  - ./bin/autoparam_AA.sh")
        print("  - ./autoparam_AA.sh")
        print(f"  - {get_bin_dir()}/autoparam_AA.sh")
        sys.exit(1)
    
    print(f"Using script: {script_path}")
    
    # Make sure script is executable
    os.chmod(script_path, 0o755)
    
    # Build command arguments for shell script
    cmd_args = [str(script_path)]
    
    for key, value in vars(args).items():
        if key == "command" or value is None:
            continue
        
        # Skip None values and empty lists
        if isinstance(value, list) and len(value) == 0:
            continue
        
        # IMPORTANT: Keep underscores for shell script (don't convert to hyphens)
        # The shell script expects underscores, not hyphens
        cmd_key = "--" + key  # Keep original underscores
        
        # Handle boolean flags
        if isinstance(value, bool):
            if value:
                # For boolean True, add the flag
                cmd_args.append(cmd_key)
        else:
            # For non-boolean, add flag and value
            cmd_args.append(cmd_key)
            cmd_args.append(str(value))
    
    # Debug: print the command if verbose
    if getattr(args, 'verbose', False):
        print(f"Running: {' '.join(cmd_args)}")
    
    # Run the script
    try:
        result = subprocess.run(cmd_args)
        sys.exit(result.returncode)
    except KeyboardInterrupt:
        print("\nInterrupted by user")
        sys.exit(1)

def run_shell_script_CG(args):
    """Run the CG parameter optimization pipeline script."""
    # First, try to find autoparam_CG.sh in various locations
    script_path = find_script("autoparam_CG.sh")
    
    # If not found, check if user provided scripts_path
    if not script_path and hasattr(args, 'scripts_path') and args.scripts_path:
        # Check if scripts_path is a directory containing the script
        scripts_path = Path(args.scripts_path)
        if scripts_path.is_dir():
            potential_script = scripts_path / "autoparam_CG.sh"
            if potential_script.exists():
                script_path = potential_script
        elif scripts_path.is_file() and scripts_path.name == "autoparam_CG.sh":
            script_path = scripts_path
    
    if not script_path:
        print(f"Error: autoparam_CG.sh not found")
        print("Please ensure the script is in one of these locations:")
        print("  - ./bin/autoparam_CG.sh")
        print("  - ./autoparam_CG.sh")
        print(f"  - {get_bin_dir()}/autoparam_CG.sh")
        print("  - Or provide --scripts_path pointing to the directory containing autoparam_CG.sh")
        sys.exit(1)
    
    print(f"Using script: {script_path}")
    
    # Make sure script is executable
    os.chmod(script_path, 0o755)
    
    # Build command arguments for shell script
    cmd_args = [str(script_path)]
    
    for key, value in vars(args).items():
        if key == "command" or value is None:
            continue
        
        # Skip internal flags
        if key in ["verbose"]:
            continue
        
        # Skip scripts_path as it's handled separately
        if key == "scripts_path":
            continue
        
        # Skip None values and empty lists
        if isinstance(value, list) and len(value) == 0:
            continue
        
        # IMPORTANT: Keep underscores for shell script (don't convert to hyphens)
        # The shell script expects underscores, not hyphens
        cmd_key = "--" + key.replace("_", "_")  # Keep underscores
        
        # Handle boolean flags
        if isinstance(value, bool):
            if value:
                # For boolean True, add the flag
                cmd_args.append(cmd_key)
            # If False, don't add anything (the script handles defaults)
        else:
            # For non-boolean, add flag and value
            cmd_args.append(cmd_key)
            cmd_args.append(str(value))
    
    # Add scripts_path if provided
    if hasattr(args, 'scripts_path') and args.scripts_path:
        cmd_args.append("--scripts_path")
        cmd_args.append(str(args.scripts_path))
    
    # Debug: print the command if verbose
    if getattr(args, 'verbose', False):
        print(f"\nRunning command:")
        print(' '.join(cmd_args))
        print()
    
    # Run the script
    try:
        # Use shell=False for security, but pass arguments properly
        result = subprocess.run(cmd_args, env=os.environ.copy())
        sys.exit(result.returncode)
    except FileNotFoundError as e:
        print(f"Error: Failed to execute {script_path}")
        print(f"Details: {e}")
        print("\nTroubleshooting:")
        print("1. Check if the script has execute permissions: chmod +x autoparam_CG.sh")
        print("2. Check if the script has a valid shebang (#!/bin/bash)")
        print("3. Check if all dependencies (gmx, python3) are in PATH")
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nInterrupted by user")
        sys.exit(1)

def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="bayesian-potentials",
        description="Bayesian Potentials - Tools for coarse-grained MD"
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Commands")
    
    # Map trajectory command
    map_parser = subparsers.add_parser("map", help="Map AA trajectory to CG")
    map_parser.add_argument("--index_cg", "-i", required=True,
                           help="CG index file for mapping")
    map_parser.add_argument("--aa_tpr", "-t", required=True)
    map_parser.add_argument("--aa_xtc", "-x", required=True)
    map_parser.add_argument("--output_mapped", "-o", default="mapped.xtc")
    map_parser.add_argument("--output_cg_gro", default="molecule.gro")
    map_parser.add_argument("--remove_pbc", action="store_true", default=True)
    map_parser.add_argument("--corrected_pbc", action="store_false", dest="remove_pbc")
    map_parser.add_argument("--verbose", "-v", action="store_true")
    map_parser.add_argument("--dry-run", action="store_true")
    map_parser.add_argument("--keep-temp", action="store_true")
    map_parser.add_argument("--gmx-cmd", default="gmx")
    
    # Generate topology command
    top_parser = subparsers.add_parser("gen-top", help="Generate CG topology")
    top_parser.add_argument("--path_ff", type=str, required=True,
                           help="Directory containing force field ITP files")
    top_parser.add_argument("--ff", type=str, default="martini_v3.0.0.itp",
                           help="Force field ITP filename (default: martini_v3.0.0.itp)")
    top_parser.add_argument("--ions", type=str, default="martini_v3.0.0_ions_v1.itp",
                           help="Ions ITP filename (default: martini_v3.0.0_ions_v1.itp)")
    top_parser.add_argument("--solvent", type=str, default="martini_v3.0.0_solvents_v1.itp",
                           help="Solvent ITP filename (default: martini_v3.0.0_solvents_v1.itp)")
    top_parser.add_argument("--itp_ligand", type=str, required=True,
                           help="Ligand ITP file to include (use 'none' to skip)")
    top_parser.add_argument("--name_molecule", type=str, required=True,
                           help="Molecule name for the ligand")
    top_parser.add_argument("--number_molecule", type=int, required=True,
                           help="Number of molecules for the ligand")
    top_parser.add_argument("--title_comments", type=str, default="",
                           help="Optional comment title at the top of the file")
    top_parser.add_argument("--title_system", type=str, default="",
                           help="Optional system title")
    top_parser.add_argument("--output_topol", type=str, required=True,
                           help="Output topology filename")
    
    # Analyze bonds/angles/dihedrals command
    analyze_parser = subparsers.add_parser("analyze", help="Calculate bonds, angles, and dihedrals from trajectory")
    analyze_parser.add_argument("--bonds_ndx", required=True,
                               help="Bonds index file")
    analyze_parser.add_argument("--angles_ndx", required=True,
                               help="Angles index file")
    analyze_parser.add_argument("--dihedrals_ndx", required=True,
                               help="Dihedrals index file")
    analyze_parser.add_argument("--xtc_file", required=True,
                               help="XTC trajectory file (can be mapped.xtc or processed trajectory)")
    analyze_parser.add_argument("--tpr_file", required=True,
                               help="TPR topology file (CG.tpr)")
    analyze_parser.add_argument("--index", default=None,
                               help="Index file for selecting groups (optional)")
    analyze_parser.add_argument("--remove_pbc", action="store_true", default=False,
                               help="Remove PBC and align trajectory (recommended for non-aligned systems)")
    analyze_parser.add_argument("--group_1", default="System",
                               help="Group for fitting (default: 'System')")
    analyze_parser.add_argument("--group_2", default="System",
                               help="Group for output (default: 'System')")
    analyze_parser.add_argument("--keep_intermediate", action="store_true", default=False,
                               help="Keep all intermediate files (whole.xtc, center_fit.xtc, reference.pdb)")
    
    # Distribution statistics command
    dist_parser = subparsers.add_parser("distributions", help="Calculate distribution statistics from XVG files")
    dist_parser.add_argument("--bonds_dir", default="bonds",
                            help="Directory containing bond XVG files (default: bonds)")
    dist_parser.add_argument("--angles_dir", default="angles",
                            help="Directory containing angle XVG files (default: angles)")
    dist_parser.add_argument("--dihedrals_dir", default="dihedrals",
                            help="Directory containing dihedral XVG files (default: dihedrals)")
    dist_parser.add_argument("--dir_to_output", default="TSV_statistics",
                            help="Directory to save statistics TSV files (default: TSV_statistics)")
    dist_parser.add_argument("--bond_out", default="bond_statistics.tsv",
                            help="Bond statistics output filename (default: bond_statistics.tsv)")
    dist_parser.add_argument("--angle_out", default="angle_statistics.tsv",
                            help="Angle statistics output filename (default: angle_statistics.tsv)")
    dist_parser.add_argument("--dihedral_out", default="dihedral_statistics.tsv",
                            help="Dihedral statistics output filename (default: dihedral_statistics.tsv)")
    
    # Adapt ITP command
    adapt_parser = subparsers.add_parser("adapt-itp", help="Adapt ITP atom names to match GRO reference")
    adapt_parser.add_argument("--input_itp", required=True,
                             help="Input ITP file to adapt")
    adapt_parser.add_argument("--input_gro_ref", required=True,
                             help="Reference GRO file with correct atom names")
    adapt_parser.add_argument("--output_itp_adapted", required=True,
                             help="Output adapted ITP file")
    adapt_parser.add_argument("--verbose", action="store_true",
                             help="Verbose output showing name changes")
    
    # Force adjustment command
    force_parser = subparsers.add_parser("force-adjust", help="Adjust force constants using Bayesian update and simulated annealing")
    force_parser.add_argument("--bonds_ref", required=True,
                             help="Reference bonds TSV file")
    force_parser.add_argument("--angles_ref", required=True,
                             help="Reference angles TSV file")
    force_parser.add_argument("--dihedrals_ref", required=True,
                             help="Reference dihedrals TSV file")
    force_parser.add_argument("--bonds_sim", required=True,
                             help="Simulated bonds TSV file")
    force_parser.add_argument("--angles_sim", required=True,
                             help="Simulated angles TSV file")
    force_parser.add_argument("--dihedrals_sim", required=True,
                             help="Simulated dihedrals TSV file")
    force_parser.add_argument("--itp_cg", required=True,
                             help="CG ITP file with topology information")
    force_parser.add_argument("--ndx_bounds", required=True,
                             help="NDX file with bond indices")
    force_parser.add_argument("--ndx_angles", required=True,
                             help="NDX file with angle indices")
    force_parser.add_argument("--ndx_dihedrals", required=True,
                             help="NDX file with dihedral indices")
    force_parser.add_argument("--molecule_name", default="molecule",
                             help="Molecule name (default: molecule)")
    
    # Boolean flags for multimodal mode (default: True)
    force_parser.add_argument("--multimodal_mode", action="store_true", default=True,
                             help="Use peak density mode instead of mean (default: True)")
    force_parser.add_argument("--no_multimodal_mode", action="store_false", dest="multimodal_mode",
                             help="Disable multimodal mode (use mean instead of peak)")
    
    # Boolean flags for variance multimodal (default: True)
    force_parser.add_argument("--variance_multimodal", action="store_true", default=True,
                             help="Use variance from all data instead of peak (default: True)")
    force_parser.add_argument("--no_variance_multimodal", action="store_false", dest="variance_multimodal",
                             help="Use variance from peak only")
    
    force_parser.add_argument("--T0", type=float, default=10.0,
                             help="Initial temperature for simulated annealing (default: 10.0)")
    force_parser.add_argument("--alpha", type=float, default=0.85,
                             help="Cooling factor for simulated annealing (default: 0.85)")
    force_parser.add_argument("--itp_out", default="cg.itp",
                             help="Output ITP file (default: cg.itp)")
    
    # CG System Setup command (bp_prep)
    cg_setup_parser = subparsers.add_parser("bp-prep", help="Setup CG system with solvent and ions")
    cg_setup_parser.add_argument("--input_ref_dir", required=True,
                                help="Directory containing input files (cg.gro, cg.itp, topol_cg.top)")
    cg_setup_parser.add_argument("--input_ff_dir", required=True,
                                help="Force field directory with Martini ITP files")
    cg_setup_parser.add_argument("--output_dir", default="cg_system",
                                help="Output directory (default: cg_system)")
    cg_setup_parser.add_argument("--input_mdp_dir", default=None,
                                help="MDP files directory (optional, will use default if not provided)")
    cg_setup_parser.add_argument("--input_name_file_mdp", default="minimization.mdp",
                                help="MDP filename (default: minimization.mdp)")
    cg_setup_parser.add_argument("--input_gro", default="cg.gro",
                                help="Input GRO filename (default: cg.gro)")
    cg_setup_parser.add_argument("--input_itp", default="cg.itp",
                                help="Input ITP filename (default: cg.itp)")
    cg_setup_parser.add_argument("--input_topol", default="topol_cg.top",
                                help="Input topology filename (default: topol_cg.top)")
    cg_setup_parser.add_argument("--pbc", default="cubic",
                                help="PBC type (default: cubic)")
    cg_setup_parser.add_argument("--sol", default="W",
                                help="Solvent type (default: W)")
    cg_setup_parser.add_argument("--salt", type=float, default=0.15,
                                help="Salt concentration in M (default: 0.15)")
    cg_setup_parser.add_argument("--skip_grompp", action="store_true",
                                help="Skip grompp step")
    cg_setup_parser.add_argument("--verbose", action="store_true",
                                help="Verbose output")
    
    # Full pipeline command for AA to CG (autoparam_aa)
    pipeline_aa_parser = subparsers.add_parser("autoparam-aa", help="Run full AA to CG parametrization pipeline")
    
    # Essential arguments (required)
    essential_group_aa = pipeline_aa_parser.add_argument_group("Essential arguments")
    essential_group_aa.add_argument("--aa_tpr", required=True, help="AA TPR file")
    essential_group_aa.add_argument("--aa_xtc", required=True, help="AA trajectory XTC file")
    essential_group_aa.add_argument("--aa_gro", required=True, help="AA GRO file")
    essential_group_aa.add_argument("--aa_itp", required=True, help="AA ITP file")
    essential_group_aa.add_argument("--beads_json", required=True, help="Beads definition JSON file")
    essential_group_aa.add_argument("--force_application", required=True, help="Force application (e.g., 'random=[1250,30;30,4]')")
    essential_group_aa.add_argument("--beads_position", required=True, choices=["com", "geom"], help="Beads position ('com' or 'geom')")
    essential_group_aa.add_argument("--input_mdp", required=True, help="MDP file for grompp")
    essential_group_aa.add_argument("--path_ff", required=True, help="Force field directory")
    
    # Optional arguments with defaults
    optional_group_aa = pipeline_aa_parser.add_argument_group("Optional arguments")
    optional_group_aa.add_argument("--output_dir", default="results", help="Output directory (default: results)")
    optional_group_aa.add_argument("--name_molecule", default="molecule", help="Molecule name (default: molecule)")
    optional_group_aa.add_argument("--number_molecule", type=int, default=1, help="Number of molecules (default: 1)")
    optional_group_aa.add_argument("--ff", default="martini_v3.0.0.itp", help="Force field ITP (default: martini_v3.0.0.itp)")
    optional_group_aa.add_argument("--ions", default="martini_v3.0.0_ions_v1.itp", help="Ions ITP (default: martini_v3.0.0_ions_v1.itp)")
    optional_group_aa.add_argument("--solvent", default="martini_v3.0.0_solvents_v1.itp", help="Solvent ITP (default: martini_v3.0.0_solvents_v1.itp)")
    optional_group_aa.add_argument("--title_comments", default="Topology system in Martini 3", help="Topology comments")
    optional_group_aa.add_argument("--title_system", default="molecule in aqueous solution", help="System title")
    optional_group_aa.add_argument("--output_topol", default="topol_cg.top", help="Output topology (default: topol_cg.top)")
    optional_group_aa.add_argument("--cycle_restr", default="fix=3,mode=cycle", 
                               help='Cycle constraints: "none", "fix=3,mode=cycle" (default), or "mode=linear"')
    optional_group_aa.add_argument("--default_martini", action="store_true", help="Use default Martini 3 masses (72) and zero charges")
    optional_group_aa.add_argument("--maxwarn", type=int, default=1, help="Max warnings for grompp (default: 1)")
    optional_group_aa.add_argument("--remove_pbc", action="store_true", default=True, help="Remove PBC (default: true)")
    optional_group_aa.add_argument("--no_pbc", action="store_false", dest="remove_pbc", help="Skip PBC removal")
    optional_group_aa.add_argument("--skip_adapt_itp", action="store_true", help="Skip ITP adaptation step")
    optional_group_aa.add_argument("--skip_grompp", action="store_true", help="Skip grompp step")
    optional_group_aa.add_argument("--skip_analysis", action="store_true", help="Skip bonds/angles/dihedrals analysis")
    optional_group_aa.add_argument("--skip_distributions", action="store_true", help="Skip distribution statistics calculation")
    optional_group_aa.add_argument("--analyze_remove_pbc", action="store_true", help="Remove PBC before analysis")
    optional_group_aa.add_argument("--analyze_group_1", default="System", help="Group for fitting in analysis (default: System)")
    optional_group_aa.add_argument("--analyze_group_2", default="System", help="Group for output in analysis (default: System)")
    optional_group_aa.add_argument("--keep_intermediate", action="store_true", help="Keep intermediate analysis files")
    optional_group_aa.add_argument("--keep_temp", action="store_true", help="Keep temporary files")
    optional_group_aa.add_argument("--verbose", action="store_true", help="Verbose output")
    
    # Distribution statistics options for pipeline
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

    # Full pipeline command for CG optimization (autoparam_cg)
    pipeline_cg_parser = subparsers.add_parser("autoparam-cg", help="Run CG parameter optimization pipeline (iterative Bayesian optimization)")
    
    # Essential arguments (required)
    essential_group_cg = pipeline_cg_parser.add_argument_group("Essential arguments")
    essential_group_cg.add_argument("--cg_md_dir", required=True,
                                   help="Directory containing CG simulation files")
    essential_group_cg.add_argument("--bonds_ref", required=True,
                                   help="Reference bonds TSV file (AA reference)")
    essential_group_cg.add_argument("--angles_ref", required=True,
                                   help="Reference angles TSV file (AA reference)")
    essential_group_cg.add_argument("--dihedrals_ref", required=True,
                                   help="Reference dihedrals TSV file (AA reference)")
    essential_group_cg.add_argument("--aa_xvg_bonds", required=True,
                                   help="Directory with AA bond XVG files (for plotting)")
    essential_group_cg.add_argument("--aa_xvg_angles", required=True,
                                   help="Directory with AA angle XVG files (for plotting)")
    essential_group_cg.add_argument("--aa_xvg_dihedrals", required=True,
                                   help="Directory with AA dihedral XVG files (for plotting)")
    essential_group_cg.add_argument("--ndx_bonds", required=True,
                                   help="NDX file with bond indices")
    essential_group_cg.add_argument("--ndx_angles", required=True,
                                   help="NDX file with angle indices")
    essential_group_cg.add_argument("--ndx_dihedrals", required=True,
                                   help="NDX file with dihedral indices")
    
    # Optional arguments with defaults
    optional_group_cg = pipeline_cg_parser.add_argument_group("Optional arguments")
    optional_group_cg.add_argument("--iterations", type=int, default=20,
                                   help="Number of iterations (default: 20)")
    optional_group_cg.add_argument("--workdir", default="calibration",
                                   help="Working directory (default: calibration)")
    optional_group_cg.add_argument("--ntomp", type=int, default=12,
                                   help="Number of OpenMP threads (default: 12)")
    optional_group_cg.add_argument("--ntmpi", type=int, default=1,
                                   help="Number of MPI threads (default: 1)")
    optional_group_cg.add_argument("--mdp_dir", default=None,
                                   help="MDP files directory (required for MD runs)")
    optional_group_cg.add_argument("--solv_ions_gro", default="solv_ions_CG.gro",
                                   help="Solvated ions GRO filename (default: solv_ions_CG.gro)")
    optional_group_cg.add_argument("--itp_to_optimize", default="cg.itp",
                                   help="ITP file to optimize (default: cg.itp)")
    optional_group_cg.add_argument("--topol_cg_file", default="topol_cg.top",
                                   help="Topology filename (default: topol_cg.top)")
    optional_group_cg.add_argument("--ff_files_dir", default="ff_files",
                                   help="Force field files directory (default: ff_files)")
    optional_group_cg.add_argument("--scripts_path", default=None,
                                   help="Path to scripts directory (auto-detected if not provided)")
    optional_group_cg.add_argument("--index_ndx", default="index.ndx",
                                   help="Index file for analysis (default: index.ndx)")
    optional_group_cg.add_argument("--group_out", default="System",
                                   help="Group for output in analysis (default: System)")
    
    # Multimodal options
    multimodal_group_cg = pipeline_cg_parser.add_argument_group("Multimodal options")
    multimodal_group_cg.add_argument("--multimodal_mode", action="store_true", default=True,
                                     help="Use peak density mode (default: True)")
    multimodal_group_cg.add_argument("--no_multimodal_mode", action="store_false", dest="multimodal_mode",
                                     help="Disable multimodal mode (use mean instead of peak)")
    multimodal_group_cg.add_argument("--variance_multimodal", action="store_true", default=True,
                                     help="Use variance from all data (default: True)")
    multimodal_group_cg.add_argument("--no_variance_multimodal", action="store_false", dest="variance_multimodal",
                                     help="Use variance from peak only")
    
    # Simulated annealing options
    sa_group_cg = pipeline_cg_parser.add_argument_group("Simulated annealing options")
    sa_group_cg.add_argument("--T0", type=float, default=10.0,
                            help="Initial temperature (default: 10.0)")
    sa_group_cg.add_argument("--alpha", type=float, default=0.85,
                            help="Cooling factor (default: 0.85)")
    
    # Output options
    output_group_cg = pipeline_cg_parser.add_argument_group("Output options")
    output_group_cg.add_argument("--output_final", default="FINAL_OPTIMIZED",
                                help="Final optimized output directory (default: FINAL_OPTIMIZED)")
    output_group_cg.add_argument("--keep_temp", action="store_true",
                                help="Keep temporary files")
    output_group_cg.add_argument("--skip_md", action="store_true",
                                help="Skip MD simulation (use existing XTC files)")
    output_group_cg.add_argument("--verbose", action="store_true",
                                help="Verbose output")
    
    # Show version
    parser.add_argument("--version", action="version", version="bayesian_potentials 0.1.0")
    
    args = parser.parse_args()
    
    if args.command == "map":
        run_mapping_script(args)
    elif args.command == "gen-top":
        run_topology_script(args)
    elif args.command == "analyze":
        run_analysis_script(args)
    elif args.command == "distributions":
        run_distribution_script(args)
    elif args.command == "adapt-itp":
        run_adaptation_script(args)
    elif args.command == "force-adjust":
        run_force_adjustment_script(args)
    elif args.command == "bp-prep":
        run_cg_setup_script(args)
    elif args.command == "autoparam-aa":
        run_shell_script_AA(args)
    elif args.command == "autoparam-cg":
        run_shell_script_CG(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
