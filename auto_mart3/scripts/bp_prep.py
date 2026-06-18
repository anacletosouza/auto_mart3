#!/usr/bin/env python3
"""
Prepare a coarse-grained (CG) molecular dynamics system using GROMACS
and the Martini 3 force field.

This pipeline automates the full setup of a CG simulation system starting
from input structure and topology files. It performs box definition,
solvation with Martini water beads, ion addition for charge neutralization
and target salt concentration, and generates input files for energy minimization.

Workflow:
    1. Copy input structure, topology, force field, and MDP files
    2. Center the molecule in the simulation box
    3. Define the simulation box:
        - Fixed size (--box_size), or
        - Distance from solute (--distance_from_atom), or
        - Default padding
    4. Solvate the system with Martini water beads
    5. Add ions (Na+/Cl−) to neutralize and reach desired concentration
    6. Update topology file with correct molecule counts
    7. Generate input file (.tpr) for energy minimization using grompp
    8. Export final system to the output directory

Usage:
    python bp-prep.py [options]

Examples:
    # Standard Martini setup (recommended)
    python bp-prep.py --use_distance_from_atom --distance_from_atom 1.2

    # Fixed cubic box
    python bp-prep.py --box_size 12.0

    # Limit number of solvent beads
    python bp-prep.py --use_distance_from_atom --max_solvent 768

Arguments (key options):
    --input_ref_dir     Directory with input GRO, ITP, and TOP files
    --input_gro         Input CG structure (.gro)
    --input_itp         Molecule topology (.itp)
    --input_topol       System topology (.top)
    --input_ff_dir      Directory containing Martini force field files

    --output_dir        Directory for final system files
    --water_dir         Directory containing water bead GRO file
    --water_file_gro    Water structure file (Martini bead)

    --ff                Force field file (default: martini_v3.0.0.itp)
    --ions              Ions file (default: martini_v3.0.0_ions_v1.itp)
    --solvent           Solvent file (default: martini_v3.0.0_solvents_v1.itp)

    --distance_from_atom  Distance (nm) between solute and box edge
    --box_size            Fixed box size (nm)
    --max_solvent         Maximum number of solvent beads (default : 2000 beads of waters)
    --solvent_radius      Martini bead radius (default: 0.21 nm)

    --salt              Salt concentration in mol/L (default: 0.15 M)

Output:
    The output directory contains:
        - solv_ions_CG.gro   Solvated and neutralized system
        - topol_cg.top       Updated topology file
        - em.tpr             Energy minimization input
        - ff_files/          Force field files
        - mdp/               Simulation parameter files
        - ions_mdp/          Ion addition parameters

Notes:
    - Requires GROMACS (gmx) available in PATH.
    - Assumes Martini 3-compatible input files.
    - Distance-based box definition (~1.2 nm) is recommended for Martini systems.
    - Temporary files are automatically removed after execution.
"""
from importlib.resources import files
import argparse
import os
import shutil
import subprocess
import re
from pathlib import Path
from pathlib import Path

data_path = Path(__file__).resolve().parents[1] / "data"

def get_package_data_path():
    return files("auto_mart3.data")

def parse_arguments():
    
    """Parse command line arguments for the CG system preparation pipeline"""
    parser = argparse.ArgumentParser(
        description='CG system preparation pipeline using GROMACS with Martini 3 parameters',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLES:
  # Basic usage with distance from atom (1.2 nm is standard for Martini)
  python3 bp-prep.py --use_distance_from_atom --distance_from_atom 1.2
  
  # Fixed box size of 12 nm cubic
  python3 bp-prep.py --box_size 12.0
  
  # With maximum solvent limit (like Martini tutorial)
  python3 bp-prep.py --use_distance_from_atom --max_solvent 768 --solvent_radius 0.21
  
  # Using custom force field file names
  python3 bp-prep.py --ff martini_v3.0.0_custom.itp --ions martini_ions.itp
  
  # Complete example with all parameters
  python3 bp-prep.py \\
    --input_ref_dir /path/to/input \\
    --input_gro cg.gro \\
    --input_itp cg.itp \\
    --input_topol topol_cg.top \\
    --input_ff_dir /path/to/ff_files \\
    --ff martini_v3.0.0.itp \\
    --ions martini_v3.0.0_ions_v1.itp \\
    --solvent martini_v3.0.0_solvents_v1.itp \\
    --output_dir /path/to/output \\
    --water_dir /path/to/water \\
    --water_file_gro water.gro \\
    --use_distance_from_atom \\
    --distance_from_atom 1.2 \\
    --solvent_radius 0.21 \\
    --max_solvent 768 \\
    --salt 0.15
        """
    )
    
    # Input/Output paths
    parser.add_argument('--input_ref_dir', 
                        default='/home/anacleto/projects/github/AA_TO_PARAMETRIZATION_CARB/MDRUN/A2/results_A2/GMX',
                        help='Directory containing input files (GRO, ITP, TOP)')
    
    parser.add_argument('--input_gro', default='cg.gro',
                        help='Input GRO filename containing CG coordinates')
    
    parser.add_argument('--input_itp', default='cg.itp',
                        help='Input ITP filename with CG topology parameters')
    
    parser.add_argument('--input_topol', default='topol_cg.top',
                        help='Input topology filename')
    
    parser.add_argument('--input_ff_dir', default=None,
                        help='Directory containing force field files (.itp). If not specified, '
                             'will look in input_ref_dir/ff_files or current directory/ff_files')
    
    # Force field file names (flexible)
    parser.add_argument('--ff', default='martini_v3.0.0.itp',
                        help='Force field ITP filename (default: martini_v3.0.0.itp)')
    
    parser.add_argument('--ions', default='martini_v3.0.0_ions_v1.itp',
                        help='Ions ITP filename (default: martini_v3.0.0_ions_v1.itp)')
    
    parser.add_argument('--solvent', default='martini_v3.0.0_solvents_v1.itp',
                        help='Solvent ITP filename (default: martini_v3.0.0_solvents_v1.itp)')
    
    parser.add_argument('--output_dir', 
                        default='/home/anacleto/projects/github/AA_TO_PARAMETRIZATION_CARB/MDRUN/A2/results_A2/OPTIMIZATION_OF_POTENTIALS/ITER_0',
                        help='Output directory for final system files')
    
    # MDP files
    parser.add_argument('--input_mdp_dir', default=str(data_path / "mdp"))
    
    parser.add_argument('--input_name_file_mdp', default='minimization.mdp',
                        help='MDP filename for energy minimization')
    
    # Water parameters
    parser.add_argument('--water_dir', default=str(data_path / "gro_files"),
                        help='Directory containing water bead GRO file')
    
    parser.add_argument('--water_file_gro', default='water.gro',
                        help='Water GRO filename (single bead representation for Martini)')
    
    # Ion parameters
    parser.add_argument('--ions_mdp_dir', default=str(data_path / "mdp"),
                        help='Directory containing MDP files for ion addition')
    
    parser.add_argument('--ions_file_mdp', default='ions.mdp',
                        help='MDP filename for ion addition step')
    
    parser.add_argument('--salt', type=float, default=0.15,
                        help='Salt concentration in Molar (default: 0.15 M, typical physiological)')
    
    # Box and solvation parameters
    parser.add_argument('--distance_from_atom', type=float, default=2.0,
                        help='Distance in nm from the molecule to add solvent (default: 2.0 nm, standard for Martini)')
    
    parser.add_argument('--use_distance_from_atom', action='store_true',
                        help='Use distance_from_atom to calculate box size (recommended for Martini)')
    
    parser.add_argument('--box_size', type=float, default=None,
                        help='Fixed box size in nm (overrides distance_from_atom if set)')
    
    parser.add_argument('--max_solvent', type=int, default=2000,
                        help='Maximum number of solvent molecules to add (maxsol flag in gmx solvate). Default : 2000 '
                             'Useful for controlling system size')
    
    parser.add_argument('--solvent_radius', type=float, default=0.21,
                        help='Radius of solvent beads in nm for Martini CG (default: 0.21 nm, standard Martini bead size)')
    
    return parser.parse_args()

def setup_directories(output_dir):
    """
    Create output directory and temporary working directory
    
    Args:
        output_dir (str): Path to output directory
        
    Returns:
        tuple: (output_path, temp_dir) as Path objects
    """
    output_path = Path(output_dir).resolve()
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Create temp directory for intermediate manipulations
    temp_dir = output_path / "temp_work"
    temp_dir.mkdir(exist_ok=True)
    
    return output_path, temp_dir

def copy_files_to_temp(args, temp_dir):
    """
    Copy necessary input files to temporary working directory
    
    Args:
        args: Command line arguments
        temp_dir (Path): Temporary directory path
        
    Returns:
        bool: True if successful, False otherwise
    """
    ref_source = Path(args.input_ref_dir).resolve()
    
    if not ref_source.exists():
        print(f"Error: {ref_source} not found")
        return False
    
    # Copy CG system files
    files_to_copy = [args.input_gro, args.input_itp, args.input_topol]
    for filename in files_to_copy:
        src = ref_source / filename
        if src.exists():
            shutil.copy2(src, temp_dir / filename)
            print(f"Copied {src} to {temp_dir / filename}")
        else:
            print(f"Warning: {src} not found")
    
    # Copy water bead file
    water_source = Path(args.water_dir).resolve() / args.water_file_gro
    if water_source.exists():
        shutil.copy2(water_source, temp_dir / 'water.gro')
        print(f"Copied water file {water_source} to {temp_dir / 'water.gro'}")
    else:
        print(f"Error: Water file {water_source} not found")
        return False
    
    return True

def copy_ff_and_mdp_to_output(args, output_dir, temp_dir):
    """
    Copy force field and MDP files to output and temporary directories
    
    Args:
        args: Command line arguments
        output_dir (Path): Output directory path
        temp_dir (Path): Temporary directory path
        
    Returns:
        bool: True if successful, False otherwise
    """
    # Locate ff_files directory in multiple possible locations
    ff_source = None
    
    # Check user-specified location
    if args.input_ff_dir:
        ff_source = Path(args.input_ff_dir).resolve()
        if not ff_source.exists():
            print(f"Warning: {ff_source} not found")
            ff_source = None
    
    # Check input reference directory
    if not ff_source:
        ff_source = Path(args.input_ref_dir).resolve() / 'ff_files'
        if ff_source.exists():
            print(f"Found ff_files in {ff_source}")
    
    # Check current working directory
    if not ff_source or not ff_source.exists():
        ff_source = Path.cwd() / 'ff_files'
        if ff_source.exists():
            print(f"Found ff_files in {ff_source}")
    
    if ff_source and ff_source.exists():
        # Copy to output directory
        ff_dest = output_dir / "ff_files"
        if ff_dest.exists():
            shutil.rmtree(ff_dest)
        shutil.copytree(ff_source, ff_dest)
        print(f"\nCopied force field directory to {ff_dest}")
        
        # Copy to temp directory for grompp
        ff_temp_dest = temp_dir / "ff_files"
        if ff_temp_dest.exists():
            shutil.rmtree(ff_temp_dest)
        shutil.copytree(ff_source, ff_temp_dest)
        print(f"Copied force field directory to {ff_temp_dest}")
        
        # Verify required force field files exist
        ff_files_required = [args.ff, args.ions, args.solvent]
        for ff_file in ff_files_required:
            if not (ff_temp_dest / ff_file).exists():
                print(f"Warning: {ff_file} not found in {ff_temp_dest}")
                print(f"  Please ensure the file exists or specify correct names with --ff, --ions, --solvent")
    else:
        print(f"\nERROR: Could not find ff_files directory!")
        print(f"Please ensure ff_files directory exists with Martini .itp files")
        return False
    
    # Copy MDP directory for minimization
    mdp_source = Path(args.input_mdp_dir).resolve()
    if mdp_source.exists():
        mdp_dest = output_dir / "mdp"
        if mdp_dest.exists():
            shutil.rmtree(mdp_dest)
        shutil.copytree(mdp_source, mdp_dest)
        print(f"Copied mdp directory to {mdp_dest}")
        
        # Copy minimization MDP to temp
        min_mdp = mdp_source / args.input_name_file_mdp
        if min_mdp.exists():
            shutil.copy2(min_mdp, temp_dir / 'minimization.mdp')
            print(f"Copied {min_mdp} to {temp_dir / 'minimization.mdp'}")
    else:
        print(f"Warning: {args.input_mdp_dir} not found")
    
    # Copy or create ions MDP
    ions_mdp_source = Path(args.ions_mdp_dir).resolve()
    if ions_mdp_source.exists():
        ions_mdp_dest = output_dir / "ions_mdp"
        if ions_mdp_dest.exists():
            shutil.rmtree(ions_mdp_dest)
        shutil.copytree(ions_mdp_source, ions_mdp_dest)
        print(f"Copied ions mdp directory to {ions_mdp_dest}")
        
        ions_mdp_file = ions_mdp_source / args.ions_file_mdp
        if ions_mdp_file.exists():
            shutil.copy2(ions_mdp_file, temp_dir / 'ions.mdp')
            print(f"Copied {ions_mdp_file} to {temp_dir / 'ions.mdp'}")
    else:
        print(f"Warning: {args.ions_mdp_dir} not found")
        # Create default ions.mdp
        with open(temp_dir / 'ions.mdp', 'w') as f:
            f.write("""; ions.mdp - used as input into grompp to generate ions.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme	= Verlet    ; Buffered neighbor searching 
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
rlist           = 1.2       ; Cut-off for making neighbor list (short range forces)
coulombtype     = cutoff    ; Treatment of long range electrostatic interactions
rcoulomb        = 1.2       ; Short-range electrostatic cut-off
rvdw            = 1.2       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
""")
        print(f"Created default ions.mdp in {temp_dir}")
    
    return True

def get_molecule_extent(gro_file):
    """
    Get minimum and maximum coordinates of molecule from GRO file
    
    Args:
        gro_file (Path): Path to GRO file
        
    Returns:
        tuple: (min_coord, max_coord) as tuples of (x,y,z) or (None, None) if error
    """
    with open(gro_file, 'r') as f:
        lines = f.readlines()
    
    if len(lines) < 3:
        return None, None
    
    num_atoms = int(lines[1].strip())
    coords = []
    
    for line in lines[2:2+num_atoms]:
        if len(line) >= 44:
            try:
                x = float(line[20:28])
                y = float(line[28:36])
                z = float(line[36:44])
                coords.append((x, y, z))
            except ValueError:
                pass
    
    if not coords:
        return None, None
    
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    
    min_coord = (min(xs), min(ys), min(zs))
    max_coord = (max(xs), max(ys), max(zs))
    
    return min_coord, max_coord

def center_and_box(temp_dir, args):
    """
    Center the molecule in the simulation box and define box size
    
    This function performs two critical operations:
    1. Centers the molecule at the origin using gmx editconf -c
    2. Defines the box size based on user specifications
    
    Args:
        temp_dir (Path): Temporary directory path
        args: Command line arguments
        
    Returns:
        bool: True if successful, False otherwise
    """
    gro_file = temp_dir / args.input_gro
    
    if not gro_file.exists():
        print(f"Error: {gro_file} not found")
        return False
    
    # Step 1: Center the molecule
    print(f"\nCentering molecule in the box...")
    cmd_center = ['gmx', 'editconf', '-f', args.input_gro, '-c', '-o', 'centered.gro']
    print(f"  {' '.join(cmd_center)}")
    
    result = subprocess.run(cmd_center, cwd=temp_dir, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error centering molecule: {result.stderr}")
        return False
    
    print("Molecule centered successfully")
    
    # Step 2: Define box size based on user preference
    centered_gro = temp_dir / 'centered.gro'
    
    # Priority: 1. Fixed box size, 2. Distance from atom, 3. Default padding
    if args.box_size:
        # Use fixed box size
        print(f"\nUsing fixed box size: {args.box_size} nm")
        cmd_box = ['gmx', 'editconf', '-f', 'centered.gro', '-o', 'boxed.gro', 
                   '-box', str(args.box_size), str(args.box_size), str(args.box_size)]
    
    elif args.use_distance_from_atom:
        # Calculate box based on distance from atoms
        min_coord, max_coord = get_molecule_extent(centered_gro)
        
        if min_coord is None:
            print("Error: Could not calculate molecule extent")
            return False
        
        # Calculate required box size: distance from farthest atom to edge
        # For a centered molecule, the farthest distance from origin is max(|min|, |max|)
        max_distance = max(abs(min_coord[0]), abs(max_coord[0]),
                          abs(min_coord[1]), abs(max_coord[1]),
                          abs(min_coord[2]), abs(max_coord[2]))
        
        box_size = 2 * (max_distance + args.distance_from_atom)
        
        print(f"\nUsing distance from atom: {args.distance_from_atom} nm")
        print(f"Farthest atom distance from center: {max_distance:.3f} nm")
        print(f"Calculated box size: {box_size:.3f} nm (cubic)")
        
        cmd_box = ['gmx', 'editconf', '-f', 'centered.gro', '-o', 'boxed.gro',
                   '-box', str(box_size), str(box_size), str(box_size)]
    
    else:
        # Default: use old padding method (backward compatibility)
        width, height, depth = get_molecule_dimensions(centered_gro)
        box_x = width + 2 * 1.5  # default 1.5 nm padding
        box_y = height + 2 * 1.5
        box_z = depth + 2 * 1.5
        print(f"\nUsing default padding (1.5 nm)")
        print(f"Molecule dimensions: {width:.3f} x {height:.3f} x {depth:.3f} nm")
        print(f"Box size: {box_x:.3f} x {box_y:.3f} x {box_z:.3f} nm")
        cmd_box = ['gmx', 'editconf', '-f', 'centered.gro', '-o', 'boxed.gro',
                   '-box', str(box_x), str(box_y), str(box_z)]
    
    print(f"\nRunning editconf for box definition:")
    print(f"  {' '.join(cmd_box)}")
    
    result = subprocess.run(cmd_box, cwd=temp_dir, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error running editconf: {result.stderr}")
        return False
    
    print("Box definition completed successfully")
    return True

def get_molecule_dimensions(gro_file):
    """
    Get dimensions of molecule from GRO file (backward compatibility)
    
    Args:
        gro_file (Path): Path to GRO file
        
    Returns:
        tuple: (width, height, depth) in nm
    """
    with open(gro_file, 'r') as f:
        lines = f.readlines()
    
    if len(lines) < 3:
        return 0.0, 0.0, 0.0
    
    num_atoms = int(lines[1].strip())
    coords = []
    
    for line in lines[2:2+num_atoms]:
        if len(line) >= 44:
            try:
                x = float(line[20:28])
                y = float(line[28:36])
                z = float(line[36:44])
                coords.append((x, y, z))
            except ValueError:
                pass
    
    if not coords:
        return 0.0, 0.0, 0.0
    
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    
    width = max(xs) - min(xs)
    height = max(ys) - min(ys)
    depth = max(zs) - min(zs)
    
    return width, height, depth

def fix_topology_includes(temp_dir, args):
    """
    Fix include paths in topology file to point to ff_files directory
    
    Args:
        temp_dir (Path): Temporary directory path
        args: Command line arguments with force field file names
        
    Returns:
        bool: True if successful, False otherwise
    """
    topol_file = temp_dir / 'topol_cg.top'
    
    if not topol_file.exists():
        print(f"Error: {topol_file} not found")
        return False
    
    with open(topol_file, 'r') as f:
        content = f.read()
    
    # Update include paths to point to ff_files directory with flexible file names
    # Replace any include of martini files with the correct ff_files path
    content = re.sub(r'#include\s+"([^/]*\.itp)"', r'#include "ff_files/\1"', content)
    content = re.sub(r'#include\s+"([^"]*martini[^"]*\.itp)"', r'#include "ff_files/\1"', content)
    content = re.sub(r'#include\s+"\.\./ff_files/', '#include "ff_files/', content)
    
    # Ensure cg.itp is included from current directory
    content = re.sub(r'#include\s+"ff_files/cg\.itp"', '#include "cg.itp"', content)
    content = re.sub(r'#include\s+"\.\./cg\.itp"', '#include "cg.itp"', content)
    
    with open(topol_file, 'w') as f:
        f.write(content)
    
    print(f"Fixed include paths in {topol_file}")
    return True

def solvate_box(temp_dir, args):
    """
    Solvate the box with water beads using Martini parameters
    
    This function uses gmx solvate with Martini-specific parameters:
    - -radius: Size of CG beads (0.21 nm for Martini)
    - -maxsol: Optional maximum number of solvent molecules
    
    Args:
        temp_dir (Path): Temporary directory path
        args: Command line arguments
        
    Returns:
        int: Number of water molecules added, or None if error
    """
    # Fix topology includes
    fix_topology_includes(temp_dir, args)
    
    # Create a clean topology for solvation (only solute)
    topol_file = temp_dir / 'topol_cg.top'
    
    # Backup original topology
    shutil.copy2(topol_file, temp_dir / 'topol_cg.top.backup')
    
    # Find molecule name from topology
    with open(topol_file, 'r') as f:
        content = f.read()
    
    molecule_name = "molecule"
    lines = content.split('\n')
    in_molecules = False
    
    for line in lines:
        if '[ molecules ]' in line:
            in_molecules = True
            continue
        if in_molecules and line.strip() and not line.strip().startswith(';'):
            parts = line.split()
            if len(parts) >= 1:
                molecule_name = parts[0]
                break
    
    print(f"Found molecule name: {molecule_name}")
    
    # Create temp topology for solvation using flexible force field file names
    temp_topol = temp_dir / 'topol_solvate.top'
    with open(temp_topol, 'w') as f:
        f.write(f'#include "ff_files/{args.ff}"\n')
        f.write(f'#include "cg.itp"\n')
        f.write(f'\n[ system ]\n')
        f.write(f'; System name\n')
        f.write(f'Solvation system\n')
        f.write(f'\n[ molecules ]\n')
        f.write(f'; Compound        #mols\n')
        f.write(f'{molecule_name}               1\n')
    
    # Build solvation command with Martini parameters
    cmd = ['gmx', 'solvate', '-cp', 'boxed.gro', '-cs', 'water.gro',
           '-o', 'solvated.gro', '-p', 'topol_solvate.top',
           '-radius', str(args.solvent_radius)]
    
    # Add maxsol if specified (useful for controlling system size)
    if args.max_solvent:
        cmd.extend(['-maxsol', str(args.max_solvent)])
        print(f"Limiting to maximum {args.max_solvent} solvent molecules")
    
    print(f"\nRunning solvation with Martini parameters:")
    print(f"  Solvent radius: {args.solvent_radius} nm")
    print(f"  {' '.join(cmd)}")
    
    result = subprocess.run(cmd, cwd=temp_dir, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error running solvate: {result.stderr}")
        return None
    
    print("Solvation completed successfully")
    
    # Parse water count from solvated topology
    with open(temp_topol, 'r') as f:
        solvated_content = f.read()
    
    w_count = 0
    in_molecules = False
    for line in solvated_content.split('\n'):
        if '[ molecules ]' in line:
            in_molecules = True
            continue
        if in_molecules and line.strip() and not line.strip().startswith(';'):
            parts = line.split()
            if len(parts) >= 2 and parts[0] == 'W':
                w_count = int(parts[1])
                break
    
    print(f"Added {w_count} water molecules")
    
    # Rename solvated gro file for next steps
    shutil.move(temp_dir / 'solvated.gro', temp_dir / 'solvated_temp.gro')
    
    return w_count

def add_ions_gromacs(temp_dir, args, w_count):
    """
    Add ions using GROMACS genion for neutralization and salt concentration
    
    This function:
    1. Creates a temporary topology with solute and water
    2. Runs grompp to prepare the tpr file
    3. Uses genion to replace water molecules with NA+ and CL- ions
    
    Args:
        temp_dir (Path): Temporary directory path
        args: Command line arguments
        w_count (int): Number of water molecules before ion addition
        
    Returns:
        tuple: (new_w_count, na_count, cl_count) or False if error
    """
    # Read original topology to find molecule name
    topol_file = temp_dir / 'topol_cg.top'
    
    with open(topol_file, 'r') as f:
        content = f.read()
    
    # Find molecule name
    in_molecules = False
    molecule_name = "molecule"
    for line in content.split('\n'):
        if '[ molecules ]' in line:
            in_molecules = True
            continue
        if in_molecules and line.strip() and not line.strip().startswith(';'):
            parts = line.split()
            if len(parts) >= 1:
                molecule_name = parts[0]
                break
    
    # Create a temporary topology with solute and water using flexible file names
    temp_topol = temp_dir / 'topol_ions.top'
    with open(temp_topol, 'w') as f:
        f.write('; Topology for ion addition\n')
        f.write(f'#include "ff_files/{args.ff}"\n')
        f.write(f'#include "ff_files/{args.ions}"\n')
        f.write(f'#include "ff_files/{args.solvent}"\n')
        f.write('#include "cg.itp"\n\n')
        
        f.write('[ system ]\n')
        f.write('; System name\n')
        f.write('CG system with water\n\n')
        
        f.write('[ molecules ]\n')
        f.write('; Compound        #mols\n')
        f.write(f'{molecule_name}               1\n')
        f.write(f'W              {w_count}\n')
    
    # Run grompp to prepare tpr for ion addition
    ions_mdp = temp_dir / 'ions.mdp'
    
    cmd_grompp = ['gmx', 'grompp', '-f', 'ions.mdp', '-c', 'solvated_temp.gro',
                  '-p', 'topol_ions.top', '-o', 'ions.tpr', '-maxwarn', '2']
    
    print(f"\nRunning grompp for ion addition:")
    print(f"  {' '.join(cmd_grompp)}")
    
    # Set GMXLIB environment variable for force field location
    env = os.environ.copy()
    ff_path = temp_dir / 'ff_files'
    if ff_path.exists():
        env['GMXLIB'] = str(ff_path)
        print(f"Set GMXLIB={env['GMXLIB']}")
    
    result = subprocess.run(cmd_grompp, cwd=temp_dir, env=env, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error running grompp: {result.stderr}")
        print(f"stdout: {result.stdout}")
        return False
    
    # Run genion to add ions
    # Input: Select water group (typically group 'W') for replacement
    genion_input = "W\n"
    
    cmd_genion = ['gmx', 'genion', '-s', 'ions.tpr', '-o', 'solv_ions_CG.gro',
                  '-p', 'topol_ions.top', '-pname', 'NA', '-nname', 'CL',
                  '-neutral', '-conc', str(args.salt)]
    
    print(f"\nRunning genion to add {args.salt} M salt:")
    print(f"  {' '.join(cmd_genion)}")
    
    result = subprocess.run(cmd_genion, cwd=temp_dir, input=genion_input,
                           capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error running genion: {result.stderr}")
        print(f"stdout: {result.stdout}")
        return False
    
    print("Ions added successfully")
    
    # Parse ion counts from output topology
    na_count = 0
    cl_count = 0
    new_w_count = 0
    
    with open(temp_dir / 'topol_ions.top', 'r') as f:
        in_molecules = False
        for line in f:
            if '[ molecules ]' in line:
                in_molecules = True
                continue
            if in_molecules and line.strip() and not line.strip().startswith(';'):
                parts = line.split()
                if len(parts) >= 2:
                    if parts[0] == 'W':
                        new_w_count = int(parts[1])
                    elif parts[0] == 'NA':
                        na_count = int(parts[1])
                    elif parts[0] == 'CL':
                        cl_count = int(parts[1])
    
    # Copy the final topology to topol_cg.top
    shutil.copy2(temp_dir / 'topol_ions.top', temp_dir / 'topol_cg.top')
    
    return new_w_count, na_count, cl_count

def update_topol_for_simulation(temp_dir, w_count, na_count, cl_count):
    """
    Final update of topology file for simulation (report only)
    
    Args:
        temp_dir (Path): Temporary directory path
        w_count (int): Final water molecule count
        na_count (int): Sodium ion count
        cl_count (int): Chloride ion count
        
    Returns:
        bool: True if successful
    """
    topol_file = temp_dir / 'topol_cg.top'
    
    if not topol_file.exists():
        print(f"Error: {topol_file} not found")
        return False
    
    print(f"\nFinal topology counts:")
    print(f"  Water molecules: {w_count}")
    print(f"  NA+ ions: {na_count}")
    print(f"  CL- ions: {cl_count}")
    
    return True

def run_grompp_minimization(args, temp_dir, output_dir):
    """
    Run grompp to prepare energy minimization input file
    
    Args:
        args: Command line arguments
        temp_dir (Path): Temporary directory path
        output_dir (Path): Output directory path
        
    Returns:
        bool: True if successful, False otherwise
    """
    mdp_file = temp_dir / 'minimization.mdp'
    topol_file = temp_dir / 'topol_cg.top'
    gro_file = temp_dir / 'solv_ions_CG.gro'
    
    # Check or create minimization MDP file
    if not mdp_file.exists():
        mdp_file = output_dir / "mdp" / args.input_name_file_mdp
        if not mdp_file.exists():
            print(f"Warning: {mdp_file} not found, creating default")
            with open(temp_dir / 'minimization.mdp', 'w') as f:
                f.write("""integrator               = steep
nsteps                   = 10000

nstxout                  = 0
nstvout                  = 0
nstfout                  = 0

cutoff-scheme            = Verlet
nstlist                  = 20
nsttcouple               = 20
nstpcouple               = 20
rlist                    = 1.35
verlet-buffer-tolerance  = -1
ns_type                  = grid
pbc                      = xyz

coulombtype              = reaction-field
rcoulomb                 = 1.1
epsilon_r                = 15	; 2.5 (with polarizable water)
epsilon_rf               = 0
vdw_type                 = cutoff
vdw-modifier             = Potential-shift-verlet
rvdw                     = 1.1

constraints              = none
constraint_algorithm     = Lincs
lincs_order              = 8
lincs_warnangle          = 90
lincs_iter               = 2
""")
            mdp_file = temp_dir / 'minimization.mdp'
    
    # Check required files
    if not topol_file.exists():
        print(f"Warning: {topol_file} not found, skipping grompp")
        return False
    
    if not gro_file.exists():
        print(f"Warning: {gro_file} not found, skipping grompp")
        return False
    
    # Set GMXLIB environment variable
    env = os.environ.copy()
    ff_path = temp_dir / 'ff_files'
    if ff_path.exists():
        env['GMXLIB'] = str(ff_path)
        print(f"Set GMXLIB={env['GMXLIB']}")
    
    cmd = [
        'gmx', 'grompp',
        '-f', str(mdp_file),
        '-p', 'topol_cg.top',
        '-c', 'solv_ions_CG.gro',
        '-o', 'em.tpr',
        '-maxwarn', '2'
    ]
    
    print(f"\nRunning grompp for energy minimization in {temp_dir}:")
    print(f"  {' '.join(cmd)}")
    
    result = subprocess.run(cmd, cwd=temp_dir, env=env, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error running grompp: {result.stderr}")
        print(f"stdout: {result.stdout}")
        return False
    
    print("grompp completed successfully")
    return True

def copy_results_to_output(temp_dir, output_dir):
    """
    Copy final results to output directory
    
    Args:
        temp_dir (Path): Temporary directory path
        output_dir (Path): Output directory path
    """
    files_to_copy = [
        'solv_ions_CG.gro',  # Solvated and neutralized system
        'topol_cg.top',       # Updated topology with correct ion counts
        'em.tpr',             # Energy minimization input
        'cg.gro',             # Original coordinates
        'cg.itp',             # Molecule topology
        'boxed.gro'           # System with box defined
    ]
    
    for filename in files_to_copy:
        src = temp_dir / filename
        if src.exists():
            dst = output_dir / filename
            shutil.copy2(src, dst)
            print(f"Copied {src} to {dst}")
    
    print(f"\n{'='*60}")
    print("FINAL OUTPUT DIRECTORY STRUCTURE:")
    print(f"{'='*60}")
    print(f"{output_dir}/")
    print(f"  ├── ff_files/          (Martini force field parameters)")
    print(f"  ├── mdp/               (Simulation parameter files)")
    print(f"  ├── ions_mdp/          (MDP files for ion addition)")
    print(f"  ├── solv_ions_CG.gro   (Solvated and neutralized system)")
    print(f"  ├── topol_cg.top       (Updated topology with correct ion counts)")
    print(f"  ├── em.tpr             (Energy minimization input)")
    print(f"  ├── cg.gro             (Original coordinates)")
    print(f"  ├── cg.itp             (Molecule topology)")
    print(f"  └── boxed.gro          (System with box defined)")

def cleanup_temp(temp_dir):
    """
    Clean up temporary directory
    
    Args:
        temp_dir (Path): Temporary directory path
    """
    try:
        shutil.rmtree(temp_dir)
        print(f"\nCleaned up temporary directory: {temp_dir}")
    except Exception as e:
        print(f"Warning: Could not clean up {temp_dir}: {e}")

def main():
    """Main pipeline function orchestrating all steps"""
    args = parse_arguments()
    
    print("=" * 60)
    print("CG System Preparation Pipeline (GROMACS-based with Martini 3)")
    print("=" * 60)
    print(f"Input arguments:")
    print(f"  input_ref_dir: {args.input_ref_dir}")
    print(f"  output_dir: {args.output_dir}")
    print(f"  water_dir: {args.water_dir}")
    print(f"  water_file: {args.water_file_gro}")
    print(f"  salt: {args.salt} M")
    print(f"  Force field files:")
    print(f"    FF: {args.ff}")
    print(f"    Ions: {args.ions}")
    print(f"    Solvent: {args.solvent}")
    if args.use_distance_from_atom:
        print(f"  Distance from atom: {args.distance_from_atom} nm")
    elif args.box_size:
        print(f"  Fixed box size: {args.box_size} nm")
    else:
        print(f"  Using default padding (1.5 nm)")
    print(f"  Solvent radius: {args.solvent_radius} nm")
    if args.max_solvent:
        print(f"  Max solvent molecules: {args.max_solvent}")
    print("=" * 60)
    
    # Setup directories
    output_dir, temp_dir = setup_directories(args.output_dir)
    print(f"\nOutput directory: {output_dir}")
    print(f"Temporary work directory: {temp_dir}")
    
    # Step 1: Copy force field and MDP files
    print("\n" + "=" * 60)
    print("STEP 1: Copying force field and MDP files")
    print("=" * 60)
    if not copy_ff_and_mdp_to_output(args, output_dir, temp_dir):
        print("Failed to copy ff_files. Exiting.")
        cleanup_temp(temp_dir)
        return
    
    # Step 2: Copy input files to working directory
    print("\n" + "=" * 60)
    print("STEP 2: Copying input files to working directory")
    print("=" * 60)
    if not copy_files_to_temp(args, temp_dir):
        print("Failed to copy files. Exiting.")
        cleanup_temp(temp_dir)
        return
    
    # Step 3: Center molecule and define simulation box
    print("\n" + "=" * 60)
    print("STEP 3: Centering molecule and defining simulation box")
    print("=" * 60)
    if not center_and_box(temp_dir, args):
        print("Failed to center molecule and define box. Exiting.")
        cleanup_temp(temp_dir)
        return
    
    # Step 4: Solvate system with water beads
    print("\n" + "=" * 60)
    print("STEP 4: Solvating system with water beads")
    print("=" * 60)
    w_count = solvate_box(temp_dir, args)
    if w_count is None:
        print("Failed to solvate system. Exiting.")
        cleanup_temp(temp_dir)
        return
    
    # Step 5: Add ions for neutralization and salt concentration
    print("\n" + "=" * 60)
    print("STEP 5: Adding ions to neutralize system")
    print("=" * 60)
    result = add_ions_gromacs(temp_dir, args, w_count)
    if not result:
        print("Failed to add ions. Exiting.")
        cleanup_temp(temp_dir)
        return
    
    w_count, na_count, cl_count = result
    
    # Step 6: Finalize topology file
    print("\n" + "=" * 60)
    print("STEP 6: Finalizing topology file")
    print("=" * 60)
    update_topol_for_simulation(temp_dir, w_count, na_count, cl_count)
    
    # Step 7: Prepare energy minimization with grompp
    print("\n" + "=" * 60)
    print("STEP 7: Preparing energy minimization (grompp)")
    print("=" * 60)
    run_grompp_minimization(args, temp_dir, output_dir)
    
    # Step 8: Copy final results to output directory
    print("\n" + "=" * 60)
    print("STEP 8: Copying final results")
    print("=" * 60)
    copy_results_to_output(temp_dir, output_dir)
    
    # Cleanup temporary directory
    cleanup_temp(temp_dir)
    
    print("\n" + "=" * 60)
    print("PIPELINE COMPLETED SUCCESSFULLY!")
    print("=" * 60)
    print(f"\nAll files are ready in: {output_dir}")
    print("\nTo run energy minimization:")
    print(f"  cd {output_dir}")
    print(f"  gmx mdrun -deffnm em -v")
    print("\nTo view the system:")
    print(f"  gmx view -f solv_ions_CG.gro")
    print("\nFor production MD simulation:")
    print(f"  gmx grompp -f production.mdp -c solv_ions_CG.gro -p topol_cg.top -o md.tpr")
    print(f"  gmx mdrun -deffnm md -v")
    print("=" * 60)

if __name__ == "__main__":
    main()
