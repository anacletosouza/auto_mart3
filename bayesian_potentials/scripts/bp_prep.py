#!/usr/bin/env python3
"""
Pipeline for CG system preparation and charge neutralization
"""

import argparse
import os
import shutil
import subprocess
import re
from pathlib import Path

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='CG system preparation pipeline')
    
    parser.add_argument('--input_ref_dir', 
                        default='/home/anacleto/projects/github/AA_TO_PARAMETRIZATION_CARB/MDRUN/A2/results_A2/GMX',
                        help='Directory containing input files')
    
    parser.add_argument('--input_gro', default='cg.gro',
                        help='Input GRO filename')
    
    parser.add_argument('--input_itp', default='cg.itp',
                        help='Input ITP filename')
    
    parser.add_argument('--input_topol', default='topol_cg.top',
                        help='Input topology filename')
    
    parser.add_argument('--input_ff_dir', default='ff_files/',
                        help='Force field directory')
    
    parser.add_argument('--output_dir', 
                        default='/home/anacleto/projects/github/AA_TO_PARAMETRIZATION_CARB/MDRUN/A2/results_A2/OPTIMIZATION_OF_POTENTIALS/ITER_0',
                        help='Output directory')
    
    parser.add_argument('--input_mdp_dir', default='../../../../mdp/',
                        help='MDP files directory')
    
    parser.add_argument('--input_name_file_mdp', default='minimization.mdp',
                        help='MDP filename')
    
    parser.add_argument('--pbc', default='cubic',
                        help='PBC type (cubic, rectangular, etc.)')
    
    parser.add_argument('--sol', default='W',
                        help='Solvent type')
    
    parser.add_argument('--salt', type=float, default=0.15,
                        help='Salt concentration (M)')
    
    return parser.parse_args()

def setup_directories(output_dir):
    """Create output directory and temp directory"""
    output_path = Path(output_dir).resolve()
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Create temp directory for manipulations
    temp_dir = output_path / "temp_work"
    temp_dir.mkdir(exist_ok=True)
    
    return output_path, temp_dir

def copy_files_to_temp(args, temp_dir):
    """Copy necessary files to temp directory"""
    ref_source = Path(args.input_ref_dir).resolve()
    
    if not ref_source.exists():
        print(f"Error: {ref_source} not found")
        return False
    
    # Copy specific files from input_ref_dir
    files_to_copy = [args.input_gro, args.input_itp, args.input_topol]
    for filename in files_to_copy:
        src = ref_source / filename
        if src.exists():
            shutil.copy2(src, temp_dir / filename)
            print(f"Copied {src} to {temp_dir / filename}")
        else:
            print(f"Warning: {src} not found")
    
    return True

def copy_ff_and_mdp_to_output(args, output_dir, temp_dir):
    """Copy ff_files and mdp directories to output directory and temp directory"""
    # Copy force field directory
    ff_source = Path(args.input_ff_dir).resolve() if args.input_ff_dir else None
    if ff_source and ff_source.exists():
        # Copy to output_dir
        ff_dest = output_dir / "ff_files"
        if ff_dest.exists():
            shutil.rmtree(ff_dest)
        shutil.copytree(ff_source, ff_dest)
        print(f"\nCopied force field directory to {ff_dest}")
        
        # Also copy to temp_dir for grompp
        ff_temp_dest = temp_dir / "ff_files"
        if ff_temp_dest.exists():
            shutil.rmtree(ff_temp_dest)
        shutil.copytree(ff_source, ff_temp_dest)
        print(f"Copied force field directory to {ff_temp_dest}")
    else:
        print(f"\nWarning: {args.input_ff_dir} not found")
    
    # Copy mdp directory
    mdp_source = Path(args.input_mdp_dir).resolve()
    if mdp_source.exists():
        mdp_dest = output_dir / "mdp"
        if mdp_dest.exists():
            shutil.rmtree(mdp_dest)
        shutil.copytree(mdp_source, mdp_dest)
        print(f"Copied mdp directory to {mdp_dest}")
    else:
        print(f"Warning: {args.input_mdp_dir} not found")
    
    return True

def run_insane(args, temp_dir):
    """Run insane command in temp directory"""
    gro_file = temp_dir / args.input_gro
    
    if not gro_file.exists():
        print(f"Error: {gro_file} not found")
        return None
    
    cmd = [
        'insane',
        '-f', args.input_gro,
        '-o', 'solv_ions_CG.gro',
        '-p', 'system.top',
        '-pbc', args.pbc,
        '-sol', args.sol,
        '-salt', str(args.salt)
    ]
    
    print(f"\nRunning insane in {temp_dir}:")
    print(f"  {' '.join(cmd)}")
    
    result = subprocess.run(cmd, cwd=temp_dir, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error running insane: {result.stderr}")
        return None
    
    print("insane completed successfully")
    return result

def parse_system_top(temp_dir):
    """Parse system.top to get ion counts"""
    system_top = temp_dir / 'system.top'
    if not system_top.exists():
        print(f"Error: {system_top} not found")
        return None, None, None
    
    w_count = None
    na_count = None
    cl_count = None
    
    with open(system_top, 'r') as f:
        for line in f:
            if line.startswith('W') and not line.startswith(';'):
                parts = line.split()
                if len(parts) >= 2:
                    w_count = int(parts[1])
            elif line.startswith('NA+') and not line.startswith(';'):
                parts = line.split()
                if len(parts) >= 2:
                    na_count = int(parts[1])
            elif line.startswith('CL-') and not line.startswith(';'):
                parts = line.split()
                if len(parts) >= 2:
                    cl_count = int(parts[1])
    
    return w_count, na_count, cl_count

def fix_gro_format(line):
    """Fix GRO format - ensure proper column alignment"""
    # GRO format: 
    # Columns 0-4: residue number (5 digits, right-aligned)
    # Columns 5-9: residue name (5 chars, left-aligned)
    # Columns 10-14: atom name (5 chars, left-aligned)
    # Columns 15-19: atom number (5 digits, right-aligned)
    # Columns 20-27: X coordinate (8.3f)
    # Columns 28-35: Y coordinate (8.3f)
    # Columns 36-43: Z coordinate (8.3f)
    
    # If line is too short, return as is
    if len(line) < 44:
        return line
    
    # Parse the line components
    try:
        # Extract residue number (first 5 chars)
        res_num = line[0:5].strip()
        if not res_num:
            res_num = "0"
        
        # Extract residue name (next 5 chars)
        res_name = line[5:10].strip()
        
        # Extract atom name (next 5 chars)
        atom_name = line[10:15].strip()
        
        # Extract atom number (next 5 chars)
        atom_num = line[15:20].strip()
        if not atom_num:
            atom_num = "0"
        
        # Extract coordinates
        x = line[20:28].strip()
        y = line[28:36].strip()
        z = line[36:44].strip()
        
        # Standardize ion names
        if res_name in ['NA+', 'Na+', 'Na', 'NA']:
            res_name = 'NA'
        if res_name in ['CL-', 'Cl-', 'Cl', 'CL']:
            res_name = 'CL'
        
        if atom_name in ['NA+', 'Na+', 'Na', 'NA']:
            atom_name = 'NA'
        if atom_name in ['CL-', 'Cl-', 'Cl', 'CL']:
            atom_name = 'CL'
        
        # Rebuild line with proper formatting
        new_line = f"{int(res_num):5d}{res_name:5s}{atom_name:5s}{int(atom_num):5d}{float(x):8.3f}{float(y):8.3f}{float(z):8.3f}"
        
        # Ensure line has exactly the right length
        if len(new_line) < 44:
            new_line = new_line.ljust(44)
        
        return new_line
    
    except (ValueError, IndexError):
        # If parsing fails, return original line
        return line

def standardize_gro_file(temp_dir):
    """Completely standardize GRO file with proper formatting"""
    gro_file = temp_dir / 'solv_ions_CG.gro'
    
    if not gro_file.exists():
        print(f"Error: {gro_file} not found")
        return False
    
    # Create backup
    backup_file = temp_dir / 'solv_ions_CG.gro.backup'
    shutil.copy2(gro_file, backup_file)
    
    with open(gro_file, 'r') as f:
        lines = f.readlines()
    
    if len(lines) < 3:
        print(f"Error: {gro_file} has invalid format")
        return False
    
    header = lines[0].rstrip('\n')
    num_atoms = int(lines[1].strip())
    atom_lines = lines[2:2+num_atoms]
    box_line = lines[2+num_atoms].rstrip('\n') if len(lines) > 2+num_atoms else ''
    
    # Fix each atom line
    fixed_lines = []
    modified = 0
    
    for i, line in enumerate(atom_lines):
        original = line.rstrip('\n')
        fixed = fix_gro_format(original)
        fixed_lines.append(fixed + '\n')
        if original != fixed:
            modified += 1
            print(f"  Fixed line {i+1}: {original[:20]}... -> {fixed[:20]}...")
    
    # Write fixed gro file
    with open(gro_file, 'w') as f:
        f.write(f"{header}\n")
        f.write(f"{num_atoms:5d}\n")
        f.writelines(fixed_lines)
        if box_line:
            f.write(f"{box_line}\n")
    
    print(f"\nFixed {modified} lines in {gro_file}")
    return True

def parse_charges_from_itp(temp_dir):
    """Parse charges from cg.itp file"""
    itp_file = temp_dir / 'cg.itp'
    
    if not itp_file.exists():
        print(f"Error: {itp_file} not found")
        return 0.0
    
    total_charge = 0.0
    in_atoms = False
    
    with open(itp_file, 'r') as f:
        for line in f:
            if '[ atoms ]' in line:
                in_atoms = True
                continue
            if in_atoms and '[' in line and 'atoms' not in line:
                break
            if in_atoms and line.strip() and not line.strip().startswith(';'):
                parts = line.split()
                if len(parts) >= 7:
                    try:
                        charge = float(parts[6])
                        total_charge += charge
                    except ValueError:
                        pass
    
    print(f"Total charge from ITP: {total_charge:.3f}")
    return total_charge

def neutralize_system(temp_dir, w_count, na_count, cl_count, total_charge):
    """Add or remove ions to neutralize the system"""
    gro_file = temp_dir / 'solv_ions_CG.gro'
    
    if not gro_file.exists():
        print(f"Error: {gro_file} not found")
        return False, 0
    
    # Calculate charge deficit
    charge_deficit = -total_charge
    
    if abs(charge_deficit) < 0.01:
        print("System is already neutral")
        return True, 0
    
    # We'll add NA+ ions (charge +1) to neutralize negative charge
    na_to_add = int(round(abs(charge_deficit)))
    
    print(f"Need to add {na_to_add} NA+ ions to neutralize (current charge: {total_charge:.3f})")
    
    # Read gro file
    with open(gro_file, 'r') as f:
        lines = f.readlines()
    
    header = lines[0].rstrip('\n')
    num_atoms = int(lines[1].strip())
    atom_lines = lines[2:2+num_atoms]
    box_line = lines[2+num_atoms].rstrip('\n') if len(lines) > 2+num_atoms else ''
    
    # Find water molecules (residue name W)
    water_indices = []
    for i, line in enumerate(atom_lines):
        if len(line) >= 10:
            residue_name = line[5:10].strip()
            if residue_name == 'W':
                water_indices.append(i)
    
    print(f"Found {len(water_indices)} water molecules")
    
    if na_to_add > len(water_indices):
        print(f"Warning: Need {na_to_add} NA+ but only {len(water_indices)} water molecules available")
        na_to_add = len(water_indices)
    
    # Replace last N water molecules with NA
    water_to_replace = water_indices[-na_to_add:]
    
    for idx in water_to_replace:
        line = atom_lines[idx]
        if len(line) >= 10:
            # Parse the line
            res_num = line[0:5].strip()
            atom_num = line[15:20].strip()
            x = line[20:28].strip()
            y = line[28:36].strip()
            z = line[36:44].strip()
            
            # Create new line with NA instead of W
            try:
                new_line = f"{int(res_num):5d}NA   NA   {int(atom_num):5d}{float(x):8.3f}{float(y):8.3f}{float(z):8.3f}"
                atom_lines[idx] = new_line + '\n'
                print(f"  Replaced water molecule at index {idx+1} with NA")
            except (ValueError, IndexError):
                print(f"  Warning: Could not parse line {idx+1}")
    
    # Write modified gro file
    with open(gro_file, 'w') as f:
        f.write(f"{header}\n")
        f.write(f"{num_atoms:5d}\n")
        f.writelines(atom_lines)
        if box_line:
            f.write(f"{box_line}\n")
    
    print(f"Added {na_to_add} NA+ ions to neutralize system")
    return True, na_to_add

def update_topol_cg(temp_dir, w_count, na_count, cl_count, na_to_add=0):
    """Update topol_cg.top with correct ion counts"""
    topol_file = temp_dir / 'topol_cg.top'
    
    if not topol_file.exists():
        print(f"Error: {topol_file} not found")
        return False
    
    # Calculate new counts
    new_w_count = w_count - na_to_add
    new_na_count = na_count + na_to_add
    new_cl_count = cl_count
    
    print(f"\nUpdating topology counts:")
    print(f"  W: {w_count} -> {new_w_count}")
    print(f"  NA: {na_count} -> {new_na_count}")
    print(f"  CL: {cl_count} -> {new_cl_count}")
    
    # Read original content
    with open(topol_file, 'r') as f:
        content = f.read()
    
    # Replace ion names
    content = content.replace('NA+', 'NA')
    content = content.replace('CL-', 'CL')
    
    # Update include paths
    content = re.sub(r'#include\s+"[^"]*martini_v3\.0\.0\.itp"', '#include "ff_files/martini_v3.0.0.itp"', content)
    content = re.sub(r'#include\s+"[^"]*martini_v3\.0\.0_ions_v1\.itp"', '#include "ff_files/martini_v3.0.0_ions_v1.itp"', content)
    content = re.sub(r'#include\s+"[^"]*martini_v3\.0\.0_solvents_v1\.itp"', '#include "ff_files/martini_v3.0.0_solvents_v1.itp"', content)
    content = re.sub(r'#include\s+"[^"]*cg\.itp"', '#include "cg.itp"', content)
    
    # Find and update [ molecules ] section
    lines = content.split('\n')
    new_lines = []
    in_molecules = False
    molecules_section_found = False
    
    for line in lines:
        if '[ molecules ]' in line:
            in_molecules = True
            molecules_section_found = True
            new_lines.append(line)
            new_lines.append('; Compound        #mols')
            new_lines.append(f'molecule               1')
            new_lines.append(f'W              {new_w_count}')
            new_lines.append(f'NA               {new_na_count}')
            new_lines.append(f'CL               {new_cl_count}')
            continue
        
        if in_molecules:
            # Skip old molecule lines
            if line.strip() and not line.strip().startswith(';') and not '[' in line:
                continue
            if '[' in line and line.strip() != '[ molecules ]':
                in_molecules = False
                new_lines.append(line)
        else:
            new_lines.append(line)
    
    if not molecules_section_found:
        # Append molecules section at the end
        new_lines.append('')
        new_lines.append('[ molecules ]')
        new_lines.append('; Compound        #mols')
        new_lines.append(f'molecule               1')
        new_lines.append(f'W              {new_w_count}')
        new_lines.append(f'NA               {new_na_count}')
        new_lines.append(f'CL               {new_cl_count}')
    
    # Write updated topology
    with open(topol_file, 'w') as f:
        f.write('\n'.join(new_lines))
    
    print(f"\nUpdated {topol_file}")
    return True

def run_grompp(args, temp_dir, output_dir):
    """Run grompp for energy minimization"""
    mdp_file = output_dir / "mdp" / args.input_name_file_mdp
    topol_file = temp_dir / 'topol_cg.top'
    gro_file = temp_dir / 'solv_ions_CG.gro'
    
    if not mdp_file.exists():
        print(f"Warning: {mdp_file} not found, skipping grompp")
        return False
    
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
        '-maxwarn', '1'
    ]
    
    print(f"\nRunning grompp in {temp_dir}:")
    print(f"  {' '.join(cmd)}")
    
    result = subprocess.run(cmd, cwd=temp_dir, env=env, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error running grompp: {result.stderr}")
        return False
    
    print("grompp completed successfully")
    return True

def copy_results_to_output(temp_dir, output_dir):
    """Copy final results to output directory"""
    files_to_copy = [
        'solv_ions_CG.gro',
        'topol_cg.top',
        'em.tpr',
        'system.top',
        'cg.gro',
        'cg.itp'
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
    print(f"  ├── solv_ions_CG.gro   (Solvated and neutralized system)")
    print(f"  ├── topol_cg.top       (Updated topology with correct ion counts)")
    print(f"  ├── em.tpr             (Energy minimization input)")
    print(f"  ├── cg.gro             (Original coordinates)")
    print(f"  ├── cg.itp             (Molecule topology)")
    print(f"  └── system.top         (Insane output for reference)")

def cleanup_temp(temp_dir):
    """Clean up temporary directory"""
    try:
        shutil.rmtree(temp_dir)
        print(f"\nCleaned up temporary directory: {temp_dir}")
    except Exception as e:
        print(f"Warning: Could not clean up {temp_dir}: {e}")

def main():
    """Main pipeline function"""
    args = parse_arguments()
    
    print("=" * 60)
    print("CG System Preparation Pipeline")
    print("=" * 60)
    print(f"Input arguments:")
    print(f"  input_ref_dir: {args.input_ref_dir}")
    print(f"  output_dir: {args.output_dir}")
    print(f"  pbc: {args.pbc}")
    print(f"  sol: {args.sol}")
    print(f"  salt: {args.salt}")
    print("=" * 60)
    
    # Setup directories
    output_dir, temp_dir = setup_directories(args.output_dir)
    print(f"\nOutput directory: {output_dir}")
    print(f"Temporary work directory: {temp_dir}")
    
    # Copy ff_files and mdp to output
    print("\n" + "=" * 60)
    print("STEP 1: Copying force field and MDP files")
    print("=" * 60)
    copy_ff_and_mdp_to_output(args, output_dir, temp_dir)
    
    # Copy files to temp
    print("\n" + "=" * 60)
    print("STEP 2: Copying input files to working directory")
    print("=" * 60)
    if not copy_files_to_temp(args, temp_dir):
        print("Failed to copy files. Exiting.")
        cleanup_temp(temp_dir)
        return
    
    # Run insane
    print("\n" + "=" * 60)
    print("STEP 3: Running insane (solvation and ionization)")
    print("=" * 60)
    if run_insane(args, temp_dir) is None:
        print("Failed to run insane. Exiting.")
        cleanup_temp(temp_dir)
        return
    
    # Parse system.top
    print("\n" + "=" * 60)
    print("STEP 4: Parsing system topology")
    print("=" * 60)
    w_count, na_count, cl_count = parse_system_top(temp_dir)
    if w_count is None:
        print("Failed to parse system.top. Exiting.")
        cleanup_temp(temp_dir)
        return
    
    print(f"\nInitial counts from system.top:")
    print(f"  Water molecules (W): {w_count}")
    print(f"  Sodium ions (NA+): {na_count}")
    print(f"  Chloride ions (CL-): {cl_count}")
    
    # Standardize GRO file
    print("\n" + "=" * 60)
    print("STEP 5: Standardizing ion names in GRO file")
    print("=" * 60)
    standardize_gro_file(temp_dir)
    
    # Parse charges from ITP
    print("\n" + "=" * 60)
    print("STEP 6: Calculating molecular charge from ITP")
    print("=" * 60)
    total_charge = parse_charges_from_itp(temp_dir)
    
    # Neutralize system if needed
    print("\n" + "=" * 60)
    print("STEP 7: Checking and adjusting system neutrality")
    print("=" * 60)
    if abs(total_charge) > 0.01:
        print(f"System has charge {total_charge:.3f}. Adding counter-ions...")
        success, na_to_add = neutralize_system(temp_dir, w_count, na_count, cl_count, total_charge)
        if success:
            update_topol_cg(temp_dir, w_count, na_count, cl_count, na_to_add)
    else:
        print("System is already neutral.")
        update_topol_cg(temp_dir, w_count, na_count, cl_count, 0)
    
    # Run grompp
    print("\n" + "=" * 60)
    print("STEP 8: Preparing energy minimization (grompp)")
    print("=" * 60)
    run_grompp(args, temp_dir, output_dir)
    
    # Copy results
    print("\n" + "=" * 60)
    print("STEP 9: Copying final results")
    print("=" * 60)
    copy_results_to_output(temp_dir, output_dir)
    
    # Cleanup
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
    print("=" * 60)

if __name__ == "__main__":
    main()
