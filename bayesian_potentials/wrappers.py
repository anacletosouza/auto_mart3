#!/usr/bin/env python3
"""
Wrappers for individual scripts to be used as console commands
"""

import sys
import subprocess
import argparse
from pathlib import Path

def get_scripts_dir():
    """Get the directory containing the scripts"""
    # Get the package root directory (where this wrapper is located)
    package_root = Path(__file__).parent
    return package_root / "scripts"

def get_data_dir():
    """Find the data directory containing the Martini3 mapping file"""
    package_root = Path(__file__).parent
    
    # Look in package_root/data/
    mapping_file = package_root / "data" / "definitions_atoms_ff_martini3.json"
    
    if mapping_file.exists():
        return str(mapping_file)
    
    # Try other possible locations
    possible_paths = [
        package_root.parent / "data" / "definitions_atoms_ff_martini3.json",
        Path.cwd() / "data" / "definitions_atoms_ff_martini3.json",
    ]
    
    for path in possible_paths:
        if path.exists():
            return str(path)
    
    return None

def run_script(script_name, args):
    """Run a script with given arguments"""
    script_path = get_scripts_dir() / script_name

    if not script_path.exists():
        print(f"Error: Script {script_name} not found at {script_path}")
        print(f"Looking in: {get_scripts_dir()}")
        sys.exit(1)

    # Build command
    cmd = [sys.executable, str(script_path)] + args

    # Run the script
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(e.returncode)
    except KeyboardInterrupt:
        print("\nInterrupted by user")
        sys.exit(1)

def gro_to_json():
    """Wrapper for 1-gro_to_beads_json.py"""
    run_script("1-gro_to_beads_json.py", sys.argv[1:])

def json_to_gro():
    """Wrapper for 2-json_to_beads_gro.py"""
    run_script("2-json_to_beads_gro.py", sys.argv[1:])

def generate_index():
    """Wrapper for 3-index_map.py"""
    run_script("3-index_map.py", sys.argv[1:])

def generate_itp():
    """Wrapper for 4-defining_atoms_type.py"""
    # Check if --def_json was provided
    args = sys.argv[1:]
    def_json_provided = False
    def_json_index = -1
    
    for i, arg in enumerate(args):
        if arg == "--def_json" and i+1 < len(args):
            def_json_provided = True
            break
    
    # If --def_json not provided, try to auto-detect
    if not def_json_provided:
        mapping_file = get_data_dir()
        if mapping_file:
            print(f"Auto-detecting Martini3 mapping file: {mapping_file}")
            # Insert --def_json and the mapping file into arguments
            args.extend(["--def_json", mapping_file])
        else:
            print("Warning: Could not auto-detect definitions_atoms_ff_martini3.json")
            print("Please provide it manually with --def_json")
    
    run_script("4-defining_atoms_type.py", args)

def generate_bonds():
    """Wrapper for 5-bonds_general.py"""
    run_script("5-bonds_general.py", sys.argv[1:])

def generate_angles():
    """Wrapper for 6-angles.py"""
    run_script("6-angles.py", sys.argv[1:])

def generate_dihedrals():
    """Wrapper for 7-dihedrals.py"""
    run_script("7-dihedrals.py", sys.argv[1:])

def generate_final():
    """Wrapper for 8-itp_final.py"""
    run_script("8-itp_final.py", sys.argv[1:])

def cg_martini3_pipeline():
    """
    Complete pipeline wrapper that runs all steps in sequence
    This is the main entry point for the cg-martini3 command
    """
    parser = argparse.ArgumentParser(
        description="Complete CG-Martini3 pipeline for carbohydrates and glycoproteins",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  cg-martini3 --input_gro carb.gro --input_beads_definitions beads_config.json \\
              --output_dir OUTPUT --beads_position geom --aa_itp carb.itp
  
  cg-martini3 --input_gro carb.gro --input_beads_definitions beads_config.json \\
              --output_dir OUTPUT --beads_position com --aa_itp carb.itp \\
              --max_distance 6.0 --tolerance 1.0
        """
    )
    
    # Required arguments
    parser.add_argument("--input_gro", required=True,
                        help="Input all-atom GRO file")
    parser.add_argument("--input_beads_definitions", required=True,
                        help="JSON configuration file for bead definitions")
    parser.add_argument("--output_dir", required=True,
                        help="Output directory for all generated files")
    parser.add_argument("--beads_position", required=True, choices=["com", "geom"],
                        help="Position calculation method: com (center of mass) or geom (geometric center)")
    parser.add_argument("--aa_itp", required=True,
                        help="All-atom ITP file with [ atoms ] section")
    
    # Optional arguments
    parser.add_argument("--definition_martini3_bead_type_to_itp", 
                        help="Optional: Path to Martini3 mapping file (auto-detected if not provided)")
    parser.add_argument("--max_distance", type=float, default=5.0,
                        help="Maximum distance for external bonds in Å (default: 5.0)")
    parser.add_argument("--tolerance", type=float, default=0.5,
                        help="Tolerance for external bonds in Å (default: 0.5)")
    parser.add_argument("--python_dir", 
                        help="Python directory (not used, kept for compatibility)")
    parser.add_argument("--keep_intermediate", action="store_true",
                        help="Keep all intermediate files (default: False)")
    parser.add_argument("--exclusion", default="all",
                        help="Exclusions: 'all', 'none', or list like '1-5,7-9,20-26' (default: 'all')")
    parser.add_argument("--force_application", default="fix=[1250;25]",
                        help="Force constants: 'fix=[bond_k;angle_k]' or 'random=[mean_bond,sd_bond;mean_angle,sd_angle]' (default: 'fix=[1250;25]')")
    
    args = parser.parse_args()
    
    # Auto-detect Martini3 mapping file if not provided
    if not args.definition_martini3_bead_type_to_itp:
        mapping_file = get_data_dir()
        if mapping_file:
            print(f"✓ Auto-detected Martini3 mapping file: {mapping_file}")
            args.definition_martini3_bead_type_to_itp = mapping_file
        else:
            print("✗ Error: Could not find definitions_atoms_ff_martini3.json")
            print("  Please provide it manually with --definition_martini3_bead_type_to_itp")
            sys.exit(1)
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Define file paths
    beads_json = output_dir / "beads.json"
    cg_gro = output_dir / "cg_beads.gro"
    ndx_file = output_dir / "cgbuilder.ndx"
    map_file = output_dir / "cgbuilder.map"
    itp_file = output_dir / "carb_cg.itp"
    bonds_file = output_dir / "bonds.ndx"
    angles_file = output_dir / "angles.ndx"
    dihedrals_file = output_dir / "dihedrals.ndx"
    
    print("\n" + "="*60)
    print("CG-Martini3 Pipeline Started")
    print("="*60)
    print(f"Input GRO:           {args.input_gro}")
    print(f"Beads definitions:   {args.input_beads_definitions}")
    print(f"Martini3 mapping:    {args.definition_martini3_bead_type_to_itp}")
    print(f"Output directory:    {args.output_dir}")
    print(f"Beads position:      {args.beads_position}")
    print(f"AA ITP:              {args.aa_itp}")
    print(f"Max distance:        {args.max_distance} Å")
    print(f"Tolerance:           {args.tolerance} Å")
    print("="*60)
    
    # Step 1: Convert GRO to JSON
    print("\nStep 1/8: Converting GRO to JSON bead definitions...")
    run_script("1-gro_to_beads_json.py", [
        "--input_gro", args.input_gro,
        "--output_json", str(beads_json),
        "--config_json", args.input_beads_definitions,
        "--position_beads", args.beads_position
    ])
    print(f"   ✓ Generated: {beads_json}")
    
    # Step 2: Convert JSON to CG GRO
    print("\nStep 2/8: Converting JSON to CG bead GRO file...")
    run_script("2-json_to_beads_gro.py", [
        "--json_beads_definition", str(beads_json),
        "--input_gro", args.input_gro,
        "--output_beads_gro", str(cg_gro),
        "--coordinate_option", args.beads_position
    ])
    print(f"   ✓ Generated: {cg_gro}")
    
    # Step 3: Generate index and map files
    print("\nStep 3/8: Generating index and map files...")
    run_script("3-index_map.py", [
        "--json_beads_definition", str(beads_json),
        "--input_gro", args.input_gro,
        "--output_ndx", str(ndx_file),
        "--output_map", str(map_file)
    ])
    print(f"   ✓ Generated: {ndx_file}")
    print(f"   ✓ Generated: {map_file}")
    
    # Step 4: Generate CG topology
    print("\nStep 4/8: Building CG topology (ITP)...")
    run_script("4-defining_atoms_type.py", [
        "--cg_gro", str(cg_gro),
        "--cg_ndx", str(ndx_file),
        "--aa_itp", args.aa_itp,
        "--def_json", args.definition_martini3_bead_type_to_itp,
        "--output", str(itp_file)
    ])
    print(f"   ✓ Generated: {itp_file}")
    
    # Step 5: Generate bonds
    print("\nStep 5/8: Detecting bonds...")
    run_script("5-bonds_general.py", [
        "--cg_gro", str(cg_gro),
        "--output_ndx", str(bonds_file),
        "--max_distance", str(args.max_distance),
        "--tolerance", str(args.tolerance)
    ])
    print(f"   ✓ Generated: {bonds_file}")
    
    # Step 6: Generate angles
    print("\nStep 6/8: Calculating angles...")
    run_script("6-angles.py", [
        "--gro", str(cg_gro),
        "--bonds", str(bonds_file),
        "--output", str(angles_file)
    ])
    print(f"   ✓ Generated: {angles_file}")
    
    # Step 7: Generate dihedrals
    print("\nStep 7/8: Calculating dihedrals...")
    run_script("7-dihedrals.py", [
        "-g", str(cg_gro),
        "-b", str(bonds_file),
        "-o", str(dihedrals_file)
    ])
    print(f"   ✓ Generated: {dihedrals_file}")
    
    # Step 8: Generate final ITP
    print("\nStep 8/8: Generating final ITP file...")
    run_script("8-itp_final.py", [
        "--input_bonds", str(bonds_file),
        "--input_angles", str(angles_file),
        "--input_dihedrals", str(dihedrals_file),
        "--itp_incomplete", str(itp_file),
        "--cg_gro", str(cg_gro),
        "--itp_complete", str(output_dir / "carb_final.itp"),
        "--exclusion", args.exclusion,
        "--force_application", args.force_application
    ])
    print(f"   ✓ Generated: {output_dir / 'carb_final.itp'}")
    
    # Summary
    print("\n" + "="*60)
    print("✓ Pipeline Completed Successfully!")
    print("="*60)
    print(f"\nOutput directory: {args.output_dir}")
    print("\nGenerated files:")
    print(f"  • {beads_json.name}        - Bead definitions")
    print(f"  • {cg_gro.name}            - CG bead coordinates")
    print(f"  • {ndx_file.name}          - Index file")
    print(f"  • {map_file.name}          - Mapping file")
    print(f"  • {itp_file.name}          - CG topology")
    print(f"  • {bonds_file.name}        - Bond definitions")
    print(f"  • {angles_file.name}       - Angle definitions")
    print(f"  • {dihedrals_file.name}    - Dihedral definitions")
    print(f"  • carb_final.itp          - Final complete topology")
    
    if not args.keep_intermediate:
        print("\nTip: Use --keep_intermediate to preserve all files")
    
    print("\nNext steps:")
    print("  1. Include carb_cg.itp in your GROMACS topology")
    print("  2. Use cg_beads.gro as the starting structure")
    print("  3. Include bonds.ndx, angles.ndx, dihedrals.ndx in your .mdp file")
    print("="*60)


if __name__ == "__main__":
    # If called directly, show available commands
    if len(sys.argv) > 1 and sys.argv[1] in ['cg-martini3', 'pipeline']:
        # Remove the first argument and call pipeline
        sys.argv.pop(0) if sys.argv[0].endswith('cg-martini3') else None
        cg_martini3_pipeline()
    else:
        print("CG-Martini3 Wrapper Module")
        print("\nAvailable commands:")
        print("  cg-martini3              - Run complete pipeline")
        print("  gro-to-json              - Step 1 only")
        print("  json-to-gro              - Step 2 only")
        print("  generate-index           - Step 3 only")
        print("  generate-itp             - Step 4 only")
        print("  generate-bonds           - Step 5 only")
        print("  generate-angles          - Step 6 only")
        print("  generate-dihedrals       - Step 7 only")
        print("  generate_final           - Step 8 only")
        print("\nExample:")
        print("  cg-martini3 --input_gro carb.gro --input_beads_definitions beads_config.json \\")
        print("              --output_dir OUTPUT --beads_position geom --aa_itp carb.itp")
