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

def run_shell_script(args):
    """Run the main shell pipeline script."""
    script_path = get_package_dir() / "bin" / "mapping_aa_to_cg.sh"
    
    if not script_path.exists():
        print(f"Error: Shell script not found at {script_path}")
        sys.exit(1)
    
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
        
        # Keep original underscore format for the shell script
        cmd_key = "--" + key
        
        if isinstance(value, list):
            cmd_args.append(cmd_key)
            for v in value:
                cmd_args.append(str(v))
        elif isinstance(value, bool):
            if value:
                cmd_args.append(cmd_key)
        else:
            cmd_args.append(cmd_key)
            cmd_args.append(str(value))
    
    # Run the script
    try:
        result = subprocess.run(cmd_args)
        sys.exit(result.returncode)
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
    
    # Full pipeline command
    pipeline_parser = subparsers.add_parser("pipeline", help="Run full pipeline")
    
    # Required arguments
    pipeline_parser.add_argument("--cg_ndx", required=True,
                                help="CG index file (used for both mapping and ITP generation)")
    pipeline_parser.add_argument("--aa_tpr", required=True,
                                help="AA TPR file")
    pipeline_parser.add_argument("--aa_xtc", required=True,
                                help="AA trajectory XTC file")
    pipeline_parser.add_argument("--output_mapped", required=True,
                                help="Output mapped trajectory")
    pipeline_parser.add_argument("--output_cg_gro", required=True,
                                help="Output CG GRO file")
    pipeline_parser.add_argument("--aa_itp", required=True,
                                help="AA ITP file")
    pipeline_parser.add_argument("--output_cg_itp", required=True,
                                help="Output CG ITP file")
    
    # Optional arguments with defaults
    pipeline_parser.add_argument("--python_scripts_dir", 
                                default=None,
                                help="Path to scripts directory (auto-detected if not provided)")
    pipeline_parser.add_argument("--input_mdp", 
                                default=None,
                                help="MDP file for grompp (uses default from package if not provided)")
    
    # Topology generation options
    pipeline_parser.add_argument("--path_ff", 
                                type=str,
                                default=None,
                                help="Directory containing force field ITP files (auto-detects if not provided)")
    pipeline_parser.add_argument("--ff", 
                                type=str, 
                                default="martini_v3.0.0.itp",
                                help="Force field ITP filename (default: martini_v3.0.0.itp)")
    pipeline_parser.add_argument("--ions", 
                                type=str, 
                                default="martini_v3.0.0_ions_v1.itp",
                                help="Ions ITP filename (default: martini_v3.0.0_ions_v1.itp)")
    pipeline_parser.add_argument("--solvent", 
                                type=str, 
                                default="martini_v3.0.0_solvents_v1.itp",
                                help="Solvent ITP filename (default: martini_v3.0.0_solvents_v1.itp)")
    pipeline_parser.add_argument("--name_molecule", 
                                type=str, 
                                default="molecule",
                                help="Molecule name for the ligand (default: molecule)")
    pipeline_parser.add_argument("--number_molecule", 
                                type=int, 
                                default=1,
                                help="Number of molecules for the ligand (default: 1)")
    pipeline_parser.add_argument("--title_comments", 
                                default="Topology system in Martini 3",
                                help="Topology comments")
    pipeline_parser.add_argument("--title_system", 
                                default="molecule in aqueous solution",
                                help="System title")
    pipeline_parser.add_argument("--output_topol", 
                                default="topol_cg.top",
                                help="Output topology filename (default: topol_cg.top)")
    pipeline_parser.add_argument("--default_martini", 
                                action="store_true",
                                help="Use default Martini 3 masses (72) and zero charges")
    
    # Other options
    pipeline_parser.add_argument("--maxwarn", 
                                type=int, 
                                default=1,
                                help="Max warnings for grompp (default: 1)")
    pipeline_parser.add_argument("--output_dir", 
                                default="results",
                                help="Output directory (default: results)")
    pipeline_parser.add_argument("--remove_pbc", 
                                action="store_true", 
                                default=True,
                                help="Remove PBC (default)")
    pipeline_parser.add_argument("--no_pbc", 
                                action="store_false", 
                                dest="remove_pbc",
                                help="Skip PBC removal")
    pipeline_parser.add_argument("--skip_grompp", 
                                action="store_true",
                                help="Skip grompp step (only do mapping and topology)")
    pipeline_parser.add_argument("--keep_temp", 
                                action="store_true",
                                help="Keep temporary files")
    pipeline_parser.add_argument("--verbose", 
                                action="store_true",
                                help="Verbose output")
    
    # Show version
    parser.add_argument("--version", action="version", version="bayesian_potentials 0.1.0")
    
    args = parser.parse_args()
    
    # Post-processing for pipeline command
    if args.command == "pipeline":
        # Auto-detect python scripts directory
        if args.python_scripts_dir is None:
            args.python_scripts_dir = str(get_scripts_dir())
        
        # Auto-detect MDP file from package data
        if args.input_mdp is None and not args.skip_grompp:
            default_mdp = get_default_mdp()
            if default_mdp:
                if args.verbose:
                    print(f"Using default MDP file: {default_mdp}")
                args.input_mdp = default_mdp
            else:
                if args.verbose:
                    print("Warning: No default MDP file found. Skipping grompp step.")
                args.skip_grompp = True
        
        # Auto-detect force field directory if not provided
        if args.path_ff is None and not args.skip_grompp:
            default_ff_dir = get_default_ff_dir()
            if default_ff_dir:
                if args.verbose:
                    print(f"Auto-detecting force field directory: {default_ff_dir}")
                args.path_ff = default_ff_dir
            else:
                if args.verbose:
                    print("Warning: Could not auto-detect force field directory.")
                    print("Topology generation will be skipped.")
                args.skip_grompp = True
    
    if args.command == "map":
        run_mapping_script(args)
    elif args.command == "gen-top":
        run_topology_script(args)
    elif args.command == "pipeline":
        run_shell_script(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
