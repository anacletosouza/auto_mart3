#!/usr/bin/env python3
"""
Usage example:

For original ITP only (no modifications) - DEFAULT

python3 5-bonds_general.py --input_aa_itp carb.itp --input_aa_gro carb.gro --input_cg_ndx cg.ndx --input_cg_gro cg.gro --output_bonds_ndx cg_bonds.ndx --output_cg_renamed cg_renamed.gro --report conection_report.txt

For cycles (fixing 3 internal members based on the distance)

python3 5-bonds_general.py --input_aa_itp carb.itp --input_aa_gro carb.gro --input_cg_ndx cg.ndx --input_cg_gro cg.gro --cycle_restr "fix=3,mode=cycle" --output_bonds_ndx cg_bonds.ndx --output_cg_renamed cg_renamed.gro --report conection_report.txt

For linear (all residues become linear chains)

python3 5-bonds_general.py --input_aa_itp carb.itp --input_aa_gro carb.gro --input_cg_ndx cg.ndx --input_cg_gro cg.gro --cycle_restr "mode=linear" --output_bonds_ndx cg_bonds.ndx --output_cg_renamed cg_renamed.gro --report conection_report.txt
"""

import argparse
import numpy as np
from collections import defaultdict
from itertools import combinations, permutations

# =========================
# ITP PARSER
# =========================
def parse_itp(file_path):
    """Parse GROMACS topology file (.itp) to extract atoms and bonds"""
    atoms = {}
    bonds = {}
    with open(file_path, 'r') as f:
        lines = f.readlines()
    section = None
    atom_count = 0
    bond_count = 0
    for line in lines:
        line = line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith("["):
            if "atoms" in line:
                section = "atoms"
            elif "bonds" in line:
                section = "bonds"
            else:
                section = None
            continue
        if section == "atoms":
            parts = line.split()
            atom_count += 1
            atoms[f"atom_{atom_count}"] = tuple(parts)
        elif section == "bonds":
            parts = line.split()
            if len(parts) >= 2:
                bond_count += 1
                bonds[f"bond_{bond_count}"] = (int(parts[0]), int(parts[1]))
    return {"atoms": atoms, "bonds": bonds}


# =========================
# GRO PARSER
# =========================
def parse_gro(file_path):
    """Parse GROMACS coordinate file (.gro) to extract atom positions"""
    atoms = {}
    with open(file_path, 'r') as f:
        lines = f.readlines()
    atom_lines = lines[2:-1]
    for i, line in enumerate(atom_lines):
        resnr   = int(line[0:5].strip())
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        atomnr  = int(line[15:20].strip())
        x = float(line[20:28].strip())
        y = float(line[28:36].strip())
        z = float(line[36:44].strip())
        atoms[f"atom_{i+1}"] = (resnr, resname, atomname, atomnr, x, y, z)
    return atoms


# =========================
# NDX PARSER
# =========================
def parse_ndx(file_path):
    """Parse GROMACS index file (.ndx) to extract bead definitions"""
    dic_map_ndx = {}
    bead_count = 0
    current_group = None
    key = None
    with open(file_path, 'r') as f:
        lines = f.readlines()
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            current_group = line.strip("[]").strip()
            bead_count += 1
            key = f"bead_{bead_count}"
            dic_map_ndx[key] = {"name": current_group, "atoms": [], "residue": None}
            continue
        if current_group is not None and key is not None:
            nums = [int(x) for x in line.split()]
            dic_map_ndx[key]["atoms"].extend(nums)
    return dic_map_ndx


# =========================
# CG GRO PARSER
# =========================
def parse_cg_gro(file_path):
    """Parse CG GROMACS coordinate file to extract bead positions"""
    dic_cg_gro = {}
    with open(file_path, 'r') as f:
        lines = f.readlines()
    bead_lines = lines[2:-1]
    
    residue_mapping = {}
    new_resnr = 1
    
    for i, line in enumerate(bead_lines):
        bead_id = f"bead_{i+1}"
        resnr   = int(line[0:5].strip())
        resname = line[5:10].strip()
        beadname = line[10:15].strip()
        beadnr  = int(line[15:20].strip())
        x = float(line[20:28].strip())
        y = float(line[28:36].strip())
        z = float(line[36:44].strip())
        
        if resnr not in residue_mapping:
            residue_mapping[resnr] = new_resnr
            new_resnr += 1
        
        new_resnr_value = residue_mapping[resnr]
        dic_cg_gro[bead_id] = (new_resnr_value, resname, beadname, beadnr, x, y, z)
    
    return dic_cg_gro


# =========================
# BUILD ATOM -> BEAD MAP
# =========================
def build_atom_to_bead_map(dic_map_ndx):
    """Create mapping from atom index to bead index"""
    atom_to_bead = {}
    for bead_key, bead_data in dic_map_ndx.items():
        for atom_idx in bead_data["atoms"]:
            atom_to_bead[atom_idx] = bead_key
    return atom_to_bead


# =========================
# VERIFY INDEX CONSISTENCY
# =========================
def verify_index_consistency(dic_aa_itp, dic_map_ndx):
    """Verify that all AA indices from ITP are covered in NDX mapping"""
    print("\n" + "="*60)
    print("INDEX CONSISTENCY VERIFICATION")
    print("="*60)
    
    aa_atom_indices = set()
    for bond_key, (a1, a2) in dic_aa_itp["bonds"].items():
        aa_atom_indices.add(a1)
        aa_atom_indices.add(a2)
    
    mapped_atom_indices = set()
    for bead_key, bead_data in dic_map_ndx.items():
        for atom_idx in bead_data["atoms"]:
            mapped_atom_indices.add(atom_idx)
    
    missing_indices = aa_atom_indices - mapped_atom_indices
    extra_indices = mapped_atom_indices - aa_atom_indices
    
    print(f"\n STATISTICS:")
    print(f"  Total AA atom indices in ITP: {len(aa_atom_indices)}")
    print(f"  Total atom indices mapped in NDX: {len(mapped_atom_indices)}")
    print(f"  Coverage: {len(mapped_atom_indices & aa_atom_indices)}/{len(aa_atom_indices)} ({100*len(mapped_atom_indices & aa_atom_indices)/len(aa_atom_indices):.1f}%)")
    
    if missing_indices:
        print(f"\n  WARNING: {len(missing_indices)} AA atom indices NOT mapped to any CG bead")
        print(f"  Missing indices: {sorted(missing_indices)}")
    else:
        print(f"\n SUCCESS: All AA atom indices are mapped to CG beads!")
    
    print("="*60 + "\n")
    return missing_indices, extra_indices


# =========================
# FIND BEST CYCLE BASED ON ITP AND DISTANCE
# =========================
def find_best_cycle_from_itp(beads, bead_positions, existing_connections, cycle_size):
    """
    Find the best cycle of specified size based on ITP connections and minimum distance
    Returns: tuple of beads forming the cycle in order, or None if not possible
    """
    if len(beads) < cycle_size:
        return None
    
    # Collect all possible edges from existing connections within these beads
    possible_edges = []
    for conn in existing_connections:
        if conn[0] in beads and conn[1] in beads:
            dist = np.linalg.norm(np.array(bead_positions[conn[0]]) - np.array(bead_positions[conn[1]]))
            possible_edges.append((conn[0], conn[1], dist))
    
    # Build adjacency list
    adjacency = defaultdict(set)
    for b1, b2, _ in possible_edges:
        adjacency[b1].add(b2)
        adjacency[b2].add(b1)
    
    # Find all cycles of size 'cycle_size' using existing connections
    cycles_found = []
    
    def find_cycles(current_cycle, start_bead, depth):
        if depth == cycle_size:
            # Check if last bead connects back to start
            if start_bead in adjacency[current_cycle[-1]]:
                cycles_found.append(current_cycle.copy())
            return
        
        last_bead = current_cycle[-1]
        for neighbor in adjacency[last_bead]:
            if neighbor not in current_cycle:
                # For cycles > 3, avoid revisiting
                if depth < cycle_size - 1 or neighbor == start_bead:
                    current_cycle.append(neighbor)
                    find_cycles(current_cycle, start_bead, depth + 1)
                    current_cycle.pop()
    
    # Try to find cycles starting from each bead
    for bead in beads:
        find_cycles([bead], bead, 1)
    
    if cycles_found:
        # Calculate perimeter for each cycle found
        best_cycle = None
        min_perimeter = float('inf')
        
        for cycle in cycles_found:
            perimeter = 0
            for i in range(len(cycle)):
                b1 = cycle[i]
                b2 = cycle[(i+1) % len(cycle)]
                dist = np.linalg.norm(np.array(bead_positions[b1]) - np.array(bead_positions[b2]))
                perimeter += dist
            if perimeter < min_perimeter:
                min_perimeter = perimeter
                best_cycle = cycle
        
        return tuple(best_cycle)
    
    # If no cycle found in existing connections, create one using minimum distances
    return find_minimum_distance_cycle(beads, bead_positions, cycle_size)


def find_minimum_distance_cycle(beads, bead_positions, cycle_size):
    """
    Find cycle with minimum perimeter using distance-based selection
    """
    if len(beads) < cycle_size:
        return None
    
    best_cycle = None
    min_perimeter = float('inf')
    
    # Try all combinations of beads to form a cycle
    for cycle_candidates in combinations(beads, cycle_size):
        # Try all permutations to find a valid cycle
        for perm in permutations(cycle_candidates):
            perimeter = 0
            for i in range(len(perm)):
                b1 = perm[i]
                b2 = perm[(i+1) % len(perm)]
                dist = np.linalg.norm(np.array(bead_positions[b1]) - np.array(bead_positions[b2]))
                perimeter += dist
            
            if perimeter < min_perimeter:
                min_perimeter = perimeter
                best_cycle = perm
    
    return best_cycle


# =========================
# BUILD LINEAR CHAIN (MINIMUM SPANNING TREE)
# =========================
def build_linear_chain(beads, bead_positions):
    """
    Build a linear chain (path) connecting all beads using minimum distances
    This creates a Minimum Spanning Tree which will be a path (linear)
    """
    if len(beads) <= 1:
        return set()
    
    if len(beads) == 2:
        return {tuple(sorted(beads))}
    
    bead_list = list(beads)
    edges = set()
    
    # Calculate all distances between beads
    distances = []
    for i, b1 in enumerate(bead_list):
        for j, b2 in enumerate(bead_list[i+1:], i+1):
            dist = np.linalg.norm(np.array(bead_positions[b1]) - np.array(bead_positions[b2]))
            distances.append((dist, b1, b2))
    
    # Sort by distance
    distances.sort(key=lambda x: x[0])
    
    # Kruskal's algorithm to build MST
    parent = {bead: bead for bead in bead_list}
    
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    
    def union(x, y):
        px, py = find(x), find(y)
        if px == py:
            return False
        parent[px] = py
        return True
    
    # Build MST
    mst_edges = []
    for dist, b1, b2 in distances:
        if union(b1, b2):
            mst_edges.append((b1, b2))
            if len(mst_edges) == len(bead_list) - 1:
                break
    
    # Convert to sorted edges
    for b1, b2 in mst_edges:
        edges.add(tuple(sorted([b1, b2])))
    
    return edges


# =========================
# CONNECT BRANCH BEADS TO CYCLE
# =========================
def connect_branches_to_cycle(cycle_beads, branch_beads, bead_positions, existing_connections):
    """
    Connect branch beads to the closest point on the cycle
    Prioritize existing connections from ITP
    Returns: set of edges connecting branches to cycle
    """
    edges = set()
    
    for branch in branch_beads:
        # First check if there are existing connections from this branch to cycle beads
        existing_cycle_connections = []
        for cycle_bead in cycle_beads:
            if tuple(sorted([branch, cycle_bead])) in existing_connections:
                dist = np.linalg.norm(np.array(bead_positions[branch]) - np.array(bead_positions[cycle_bead]))
                existing_cycle_connections.append((cycle_bead, dist))
        
        if existing_cycle_connections:
            # Use existing connection with minimum distance
            closest_cycle_bead = min(existing_cycle_connections, key=lambda x: x[1])[0]
            edges.add(tuple(sorted([branch, closest_cycle_bead])))
        else:
            # Find closest bead in the cycle based on distance
            min_dist = float('inf')
            closest_cycle_bead = None
            
            for cycle_bead in cycle_beads:
                dist = np.linalg.norm(np.array(bead_positions[branch]) - np.array(bead_positions[cycle_bead]))
                if dist < min_dist:
                    min_dist = dist
                    closest_cycle_bead = cycle_bead
            
            if closest_cycle_bead:
                edges.add(tuple(sorted([branch, closest_cycle_bead])))
    
    return edges


# =========================
# BUILD INTERNAL CONNECTIONS
# =========================
def build_internal_connections(beads, bead_positions, aa_connections, cycle_size=3, cycle_mode='none'):
    """
    Build internal connections based on physico-chemical principles:
    - If cycle_mode='cycle': Find optimal cycle of specified size (default 3)
    - Connect all other beads to the cycle at closest points
    - If cycle_mode='linear': Create linear chain (MST) for ALL residues
    - If cycle_mode='none' or not specified: Preserve all ITP connections
    """
    if len(beads) <= 1:
        return set()
    
    if len(beads) == 2:
        return {tuple(sorted(beads))}
    
    edges = set()
    
    # Collect existing connections within these beads
    existing_connections = set()
    for conn in aa_connections:
        if conn[0] in beads and conn[1] in beads:
            existing_connections.add(tuple(sorted([conn[0], conn[1]])))
    
    # Mode 'none' - preserve all ITP connections
    if cycle_mode == 'none':
        for conn in existing_connections:
            edges.add(conn)
        return edges
    
    # Mode 'linear' - create linear chain for ALL residues (ignore cycle_size)
    if cycle_mode == 'linear':
        return build_linear_chain(beads, bead_positions)
    
    # Mode 'cycle' - apply cycle logic
    if cycle_mode == 'cycle':
        # If number of beads is less than or equal to cycle_size, preserve ITP connections
        if len(beads) <= cycle_size:
            # Preserve all existing ITP connections
            for conn in existing_connections:
                edges.add(conn)
            
            # For exactly 3 beads, ensure they form a triangle if ITP doesn't specify
            if len(beads) == 3 and len(edges) < 3:
                # Add missing edges to complete triangle based on minimum distance
                bead_list = list(beads)
                for i in range(3):
                    for j in range(i+1, 3):
                        pair = tuple(sorted([bead_list[i], bead_list[j]]))
                        if pair not in edges:
                            edges.add(pair)
        else:
            # Find best cycle based on ITP and distance
            best_cycle = find_best_cycle_from_itp(beads, bead_positions, existing_connections, cycle_size)
            
            if best_cycle:
                # Add cycle edges
                for i in range(len(best_cycle)):
                    b1 = best_cycle[i]
                    b2 = best_cycle[(i+1) % len(best_cycle)]
                    edges.add(tuple(sorted([b1, b2])))
                
                # Identify branch beads (not in cycle)
                cycle_set = set(best_cycle)
                branch_beads = [b for b in beads if b not in cycle_set]
                
                # Connect branches to cycle, prioritizing ITP connections
                branch_edges = connect_branches_to_cycle(best_cycle, branch_beads, bead_positions, existing_connections)
                edges.update(branch_edges)
                
                print(f"    Created cycle of size {cycle_size} with {len(branch_beads)} branch beads")
            else:
                # Fallback: preserve ITP connections
                for conn in existing_connections:
                    edges.add(conn)
                
                # If still no edges, use linear chain
                if len(edges) == 0:
                    edges = build_linear_chain(beads, bead_positions)
    
    return edges


# =========================
# FIND EXTERNAL CONNECTIONS (Based on AA topology only)
# =========================
def find_external_connections(dic_joined, aa_connections, cycle_size=3, cycle_mode='none'):
    """
    Find external connections between residues based on AA topology
    """
    # Get bead positions
    bead_positions = {}
    for bead_key, bead_data in dic_joined.items():
        bead_positions[bead_key] = (bead_data["x"], bead_data["y"], bead_data["z"])
    
    # Group by residue
    residue_beads = defaultdict(list)
    for bead_key, bead_data in dic_joined.items():
        residue_beads[bead_data["resnr"]].append(bead_key)
    
    # Track all connections
    all_connections = set()
    
    # First, add all AA-based connections
    for conn in aa_connections:
        all_connections.add(conn)
    
    # Process internal connections for each residue
    print("\n Building internal connections for residues:")
    for resnr, beads in sorted(residue_beads.items()):
        if len(beads) >= 2:
            if cycle_mode == 'linear':
                print(f"  Residue {resnr}: {len(beads)} beads - building linear chain")
            elif cycle_mode == 'cycle':
                print(f"  Residue {resnr}: {len(beads)} beads, cycle_size={cycle_size}")
            else:
                print(f"  Residue {resnr}: {len(beads)} beads - preserving ITP connections")
            
            # Build new internal connections
            new_internal = build_internal_connections(beads, bead_positions, all_connections, 
                                                      cycle_size, cycle_mode)
            
            # Remove old internal connections
            to_remove = []
            for conn in all_connections:
                if conn[0] in beads and conn[1] in beads:
                    to_remove.append(conn)
            
            for conn in to_remove:
                all_connections.discard(conn)
            
            # Add new internal connections
            for conn in new_internal:
                all_connections.add(conn)
            
            print(f"    Created {len(new_internal)} internal bonds")
    
    return all_connections


# =========================
# JOIN CG GRO + NDX
# =========================
def join_cg_data(dic_map_ndx, dic_cg_gro):
    """Join CG index data with CG coordinate data"""
    dic_joined = {}
    for bead_key in dic_map_ndx.keys():
        if bead_key not in dic_cg_gro:
            continue
        ndx_data = dic_map_ndx[bead_key]
        gro_data = dic_cg_gro[bead_key]
        dic_joined[bead_key] = {
            "bead_name": ndx_data["name"],
            "aa_atoms": ndx_data["atoms"],
            "resnr": gro_data[0],
            "resname": gro_data[1],
            "cg_name": gro_data[2],
            "cg_index": gro_data[3],
            "x": gro_data[4],
            "y": gro_data[5],
            "z": gro_data[6],
        }
    return dic_joined


# =========================
# COUNT INTERNAL / EXTERNAL BONDS
# =========================
def count_internal_external(dic_joined, bead_connections):
    """Count internal (within residue) and external (between residues) bonds"""
    internal = 0
    external = 0
    internal_pairs = []
    external_pairs = []
    
    for b1, b2 in bead_connections:
        if b1 == b2:
            continue
            
        r1 = dic_joined[b1]["resnr"]
        r2 = dic_joined[b2]["resnr"]
        if r1 == r2:
            internal += 1
            internal_pairs.append((b1, b2))
        else:
            external += 1
            external_pairs.append((b1, b2))
    
    return internal, external, internal_pairs, external_pairs


# =========================
# FILTER VALID BONDS
# =========================
def filter_valid_bonds(connections):
    """Remove self-connections and duplicate bonds"""
    valid_bonds = set()
    for b1, b2 in connections:
        if b1 != b2:
            valid_bonds.add(tuple(sorted([b1, b2])))
    return valid_bonds


# =========================
# WRITE BONDS NDX FILE
# =========================
def write_bonds_ndx(filename, dic_joined, internal_pairs, external_pairs):
    """Write bonds.ndx file with internal and external bonds"""
    bead_to_index = {}
    for i, bead_key in enumerate(dic_joined.keys(), start=1):
        bead_to_index[bead_key] = i
    
    index_to_info = {}
    for bead_key, bead_data in dic_joined.items():
        idx = bead_to_index[bead_key]
        index_to_info[idx] = {
            "resnr": bead_data["resnr"],
            "bead_name": bead_data["bead_name"],
            "cg_name": bead_data["cg_name"]
        }
    
    def calculate_distance(bead1_key, bead2_key, dic_joined):
        x1, y1, z1 = dic_joined[bead1_key]["x"], dic_joined[bead1_key]["y"], dic_joined[bead1_key]["z"]
        x2, y2, z2 = dic_joined[bead2_key]["x"], dic_joined[bead2_key]["y"], dic_joined[bead2_key]["z"]
        dist_nm = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
        dist_angstrom = dist_nm * 10.0
        return dist_angstrom
    
    def is_valid_bond(bead1_key, bead2_key):
        return bead1_key != bead2_key
    
    def get_sort_key(bond_pair):
        bead1_key, bead2_key = bond_pair
        idx1 = bead_to_index[bead1_key]
        idx2 = bead_to_index[bead2_key]
        res1 = index_to_info[idx1]["resnr"]
        res2 = index_to_info[idx2]["resnr"]
        return (min(res1, res2), min(idx1, idx2))
    
    filtered_internal = [(b1, b2) for b1, b2 in internal_pairs if is_valid_bond(b1, b2)]
    filtered_external = [(b1, b2) for b1, b2 in external_pairs if is_valid_bond(b1, b2)]
    
    sorted_internal = sorted(filtered_internal, key=get_sort_key)
    sorted_external = sorted(filtered_external, key=get_sort_key)
    
    with open(filename, 'w') as f:
        f.write("[ bonds ]\n")
        f.write(";;;;;;; internal bonds\n\n")
        
        for bead1_key, bead2_key in sorted_internal:
            idx1 = bead_to_index[bead1_key]
            idx2 = bead_to_index[bead2_key]
            dist = calculate_distance(bead1_key, bead2_key, dic_joined)
            res1 = index_to_info[idx1]["resnr"]
            res2 = index_to_info[idx2]["resnr"]
            cg_name1 = index_to_info[idx1]["cg_name"]
            cg_name2 = index_to_info[idx2]["cg_name"]
            f.write(f"{idx1:6d}{idx2:6d} ; {res1}_{cg_name1}---{res2}_{cg_name2} dist={dist:.3f}\n")
        
        f.write("\n;;;;;;; external bonds\n\n")
        for bead1_key, bead2_key in sorted_external:
            idx1 = bead_to_index[bead1_key]
            idx2 = bead_to_index[bead2_key]
            dist = calculate_distance(bead1_key, bead2_key, dic_joined)
            res1 = index_to_info[idx1]["resnr"]
            res2 = index_to_info[idx2]["resnr"]
            cg_name1 = index_to_info[idx1]["cg_name"]
            cg_name2 = index_to_info[idx2]["cg_name"]
            f.write(f"{idx1:6d}{idx2:6d} ; {res1}_{cg_name1}---{res2}_{cg_name2} dist={dist:.3f}\n")


# =========================
# WRITE DETAILED REPORT
# =========================
def write_report(filename, bead_connections, dic_joined, dic_aa_itp, dic_map_ndx, dic_cg_gro,
                internal, external, internal_pairs, external_pairs):
    """Write detailed connectivity report with ITP and CG information"""
    
    # Build atom to bead mapping
    atom_to_bead = build_atom_to_bead_map(dic_map_ndx)
    
    # Build bead to atoms mapping with detailed info
    bead_to_atoms = {}
    for bead_key, bead_data in dic_map_ndx.items():
        bead_to_atoms[bead_key] = bead_data["atoms"]
    
    # Get AA atom names from ITP
    aa_atom_names = {}
    for atom_key, atom_data in dic_aa_itp["atoms"].items():
        # Parse atom index from the tuple
        if len(atom_data) >= 6:
            try:
                atom_idx = int(atom_data[0])
                atom_name = atom_data[4]  # Atom name is typically at position 4
                aa_atom_names[atom_idx] = atom_name
            except (ValueError, IndexError):
                pass
    
    with open(filename, 'w') as f:
        f.write("=== Bead Connectivity Report ===\n\n")
        
        # Group beads by residue
        residue_beads = defaultdict(list)
        for bead_key, bead_data in dic_joined.items():
            residue_beads[bead_data["resnr"]].append(bead_key)
        
        # Get all ITP connections grouped by residue
        itp_connections_by_residue = defaultdict(list)
        aa_connections, _ = find_bead_connections_aa(dic_aa_itp, dic_map_ndx)
        
        # Organize connections by residue
        for conn in aa_connections:
            b1, b2 = conn
            if b1 in dic_joined and b2 in dic_joined:
                res1 = dic_joined[b1]["resnr"]
                res2 = dic_joined[b2]["resnr"]
                if res1 == res2:
                    itp_connections_by_residue[res1].append(conn)
        
        f.write("=" * 100 + "\n")
        f.write("DETAILED RESIDUE CONNECTIVITY ANALYSIS\n")
        f.write("=" * 100 + "\n\n")
        
        for resnr, beads in sorted(residue_beads.items()):
            f.write(f"\n{'='*80}\n")
            f.write(f"RESIDUE {resnr}\n")
            f.write(f"{'='*80}\n\n")
            
            # Header for the table
            f.write(f"{'----- ' + str(resnr) + ' -----':<20} ||----- group of AA model-----||----- CG model-----\n")
            f.write(f"{'-'*80}\n")
            
            # Get ITP connections for this residue
            itp_conns = itp_connections_by_residue.get(resnr, [])
            
            # Get final connections for this residue
            final_conns = [conn for conn in bead_connections 
                          if conn[0] in beads and conn[1] in beads]
            
            # Process each connection
            all_processed_pairs = set()
            
            # First, list all ITP connections
            if itp_conns:
                f.write("\nITP (AA-level) connections:\n")
                for b1, b2 in itp_conns:
                    pair_key = tuple(sorted([b1, b2]))
                    all_processed_pairs.add(pair_key)
                    
                    # Get bead names
                    bead1_name = dic_joined[b1]["bead_name"]
                    bead2_name = dic_joined[b2]["bead_name"]
                    
                    # Get CG names
                    cg1_name = dic_joined[b1]["cg_name"]
                    cg2_name = dic_joined[b2]["cg_name"]
                    
                    # Get AA atom indices for each bead
                    atoms1 = bead_to_atoms.get(b1, [])
                    atoms2 = bead_to_atoms.get(b2, [])
                    
                    # Get CG indices
                    cg_idx1 = dic_joined[b1]["cg_index"]
                    cg_idx2 = dic_joined[b2]["cg_index"]
                    
                    # Format the output line
                    line = f"  bead_{bead1_name} - bead_{bead2_name}"
                    line += f" ==== bead_{bead1_name}[{','.join(map(str, atoms1))}] - bead_{bead2_name}[{','.join(map(str, atoms2))}]"
                    line += f" ==== {cg1_name}[{cg_idx1}] - {cg2_name}[{cg_idx2}]"
                    f.write(line + "\n")
            
            # Then, list final connections (including new ones from cycle formation)
            if final_conns:
                new_conns = [conn for conn in final_conns if tuple(sorted(conn)) not in all_processed_pairs]
                removed_conns = [conn for conn in itp_conns if tuple(sorted(conn)) not in [tuple(sorted(c)) for c in final_conns]]
                
                if new_conns:
                    f.write("\nFinal CG connections (after optimization):\n")
                    for b1, b2 in new_conns:
                        # Get bead names
                        bead1_name = dic_joined[b1]["bead_name"]
                        bead2_name = dic_joined[b2]["bead_name"]
                        
                        # Get CG names
                        cg1_name = dic_joined[b1]["cg_name"]
                        cg2_name = dic_joined[b2]["cg_name"]
                        
                        # Get AA atom indices for each bead
                        atoms1 = bead_to_atoms.get(b1, [])
                        atoms2 = bead_to_atoms.get(b2, [])
                        
                        # Get CG indices
                        cg_idx1 = dic_joined[b1]["cg_index"]
                        cg_idx2 = dic_joined[b2]["cg_index"]
                        
                        # Format the output line
                        line = f"  bead_{bead1_name} - bead_{bead2_name}"
                        line += f" ==== bead_{bead1_name}[{','.join(map(str, atoms1))}] - bead_{bead2_name}[{','.join(map(str, atoms2))}]"
                        line += f" ==== {cg1_name}[{cg_idx1}] - {cg2_name}[{cg_idx2}]"
                        
                        # Indicate if this is a cycle bond or branch connection
                        if len(final_conns) != len(itp_conns):
                            line += " [NEW - from optimization]"
                        
                        f.write(line + "\n")
                
                if removed_conns:
                    f.write("\nRemoved ITP connections (not in final model):\n")
                    for b1, b2 in removed_conns:
                        bead1_name = dic_joined[b1]["bead_name"]
                        bead2_name = dic_joined[b2]["bead_name"]
                        f.write(f"  bead_{bead1_name} - bead_{bead2_name} [REMOVED]\n")
            
            # Summary for this residue
            f.write(f"\nSummary: {len(beads)} beads, {len(final_conns)} internal bonds")
            if len(final_conns) != len(itp_conns):
                f.write(f" (original ITP had {len(itp_conns)} bonds)")
            f.write("\n")
        
        # Overall statistics
        f.write(f"\n{'='*80}\n")
        f.write("OVERALL STATISTICS\n")
        f.write(f"{'='*80}\n\n")
        f.write(f"Total unique bead connections: {len(bead_connections)}\n")
        f.write(f"Internal bonds (within residue): {internal}\n")
        f.write(f"External bonds (between residues): {external}\n")
        
        if external_pairs:
            f.write("\n=== External Connections ===\n")
            for b1, b2 in external_pairs:
                bead1_name = dic_joined[b1]["bead_name"]
                bead2_name = dic_joined[b2]["bead_name"]
                res1 = dic_joined[b1]["resnr"]
                res2 = dic_joined[b2]["resnr"]
                f.write(f"Residue {res1}_{bead1_name} ----- Residue {res2}_{bead2_name}\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("END OF REPORT\n")
        f.write("="*80 + "\n")


# =========================
# WRITE GRO
# =========================
def write_gro(filename, atoms, box):
    """Write GROMACS coordinate file (.gro)"""
    with open(filename, 'w') as f:
        f.write("Structure written by Python\n")
        f.write(f"{len(atoms)}\n")
        for a in atoms:
            line = f"{a[7]:5d}{a[1]:<5s}{a[2]:>5s}{a[3]:5d}{a[4]:8.3f}{a[5]:8.3f}{a[6]:8.3f}\n"
            f.write(line)
        f.write(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}\n")


# =========================
# PREPARE ATOMS LIST FOR WRITE_GRO
# =========================
def prepare_atoms_list(dic_joined):
    """Prepare atoms list for GRO file writing"""
    atoms_list = []
    for i, bead_key in enumerate(dic_joined.keys(), start=1):
        d = dic_joined[bead_key]
        atoms_list.append((
            d["resnr"], d["bead_name"], d["cg_name"], i,
            d["x"], d["y"], d["z"], i
        ))
    return atoms_list


# =========================
# FIND BEAD CONNECTIONS (AA-based)
# =========================
def find_bead_connections_aa(dic_aa_itp, dic_map_ndx):
    """Identify bead connections based on AA topology"""
    bonds = dic_aa_itp["bonds"]
    atom_to_bead = build_atom_to_bead_map(dic_map_ndx)
    bead_connections = set()
    connection_counts = defaultdict(int)
    
    for bond in bonds.values():
        a1, a2 = bond
        bead1 = atom_to_bead.get(a1)
        bead2 = atom_to_bead.get(a2)
        if bead1 is None or bead2 is None:
            continue
        if bead1 != bead2:
            pair = tuple(sorted([bead1, bead2]))
            bead_connections.add(pair)
            connection_counts[pair] += 1
    
    return bead_connections, connection_counts


# =========================
# MAIN
# =========================
def main():
    parser = argparse.ArgumentParser(description='Generate CG bonds from AA topology with intelligent cycle detection')
    parser.add_argument("--input_aa_itp", required=True, help="AA topology file (.itp)")
    parser.add_argument("--input_aa_gro", required=True, help="AA coordinate file (.gro)")
    parser.add_argument("--input_cg_ndx", required=True, help="CG index file (.ndx)")
    parser.add_argument("--input_cg_gro", required=True, help="CG coordinate file (.gro)")
    parser.add_argument("--output_cg_renamed", default=None, help="Output CG GRO file (optional)")
    parser.add_argument("--output_bonds_ndx", default="bonds.ndx", help="Output bonds index file")
    parser.add_argument("--box", nargs=3, type=float, default=[0.0,0.0,0.0], help="Box dimensions")
    parser.add_argument("--report", default="connectivity_report.txt", help="Output report file")
    parser.add_argument("--cycle_restr", default="none", 
                       help="Cycle restriction: fix=N (default=3), mode=cycle/linear, or none (default)")
    parser.add_argument("--skip_verification", action="store_true", help="Skip index verification")
    args = parser.parse_args()

    # Parse cycle_restr argument
    cycle_size = 3  # Default (only used for cycle mode)
    cycle_mode = 'none'  # Default is 'none' (preserve original ITP)
    
    if args.cycle_restr and args.cycle_restr.lower() != "none":
        parts = args.cycle_restr.split(',')
        for part in parts:
            if 'fix=' in part:
                try:
                    cycle_size = int(part.split('=')[1].strip('[]'))
                except:
                    pass
            elif 'mode=' in part:
                cycle_mode = part.split('=')[1].strip().lower()
        
        # If mode is specified without fix, use default cycle_size for cycle mode
        if cycle_mode == 'cycle' and 'fix=' not in args.cycle_restr:
            cycle_size = 3
            
        print(f"\n  CONFIGURATION:")
        print(f"  Mode: {cycle_mode}")
        if cycle_mode == 'cycle':
            print(f"  Cycle size: {cycle_size}")
        elif cycle_mode == 'linear':
            print(f"  Building linear chains for all residues")
    else:
        print(f"\n  CONFIGURATION: Using original ITP connections only (mode=none)")
    
    # Parse all inputs
    print("\n Reading input files...")
    dic_aa_itp = parse_itp(args.input_aa_itp)
    print(f"  ✓ AA ITP: {len(dic_aa_itp['atoms'])} atoms, {len(dic_aa_itp['bonds'])} bonds")
    
    dic_aa_gro = parse_gro(args.input_aa_gro)
    print(f"  ✓ AA GRO: {len(dic_aa_gro)} atoms")
    
    dic_map_ndx = parse_ndx(args.input_cg_ndx)
    print(f"  ✓ CG NDX: {len(dic_map_ndx)} beads")
    
    dic_cg_gro = parse_cg_gro(args.input_cg_gro)
    print(f"  ✓ CG GRO: {len(dic_cg_gro)} beads")
    
    # Verify index consistency
    if not args.skip_verification:
        verify_index_consistency(dic_aa_itp, dic_map_ndx)
    
    # Join CG
    dic_joined = join_cg_data(dic_map_ndx, dic_cg_gro)
    print(f"\n✓ Joined CG data: {len(dic_joined)} beads successfully matched")
    
    # Find AA-based bead connections
    aa_connections, _ = find_bead_connections_aa(dic_aa_itp, dic_map_ndx)
    print(f"✓ AA-based connections: {len(aa_connections)} unique bead pairs")
    
    # Find all connections with cycle detection or use original only
    if cycle_mode == 'none':
        all_connections = aa_connections.copy()
        print("\n✓ Using original ITP connections only (no modifications)")
    else:
        all_connections = find_external_connections(dic_joined, aa_connections, 
                                                    cycle_size, cycle_mode)
    
    # Filter out invalid bonds
    all_connections = filter_valid_bonds(all_connections)
    print(f"\n✓ Final connections after filtering: {len(all_connections)}")
    
    # Count internal and external
    internal, external, internal_pairs, external_pairs = count_internal_external(
        dic_joined, all_connections
    )
    
    print(f"\n FINAL STATISTICS:")
    print(f"  Internal bonds: {internal}")
    print(f"  External bonds: {external}")
    print(f"  Total bonds: {len(all_connections)}")
    
    # Write output files
    write_bonds_ndx(args.output_bonds_ndx, dic_joined, internal_pairs, external_pairs)
    print(f"\n✓ Bonds NDX written to {args.output_bonds_ndx}")
    
    write_report(args.report, all_connections, dic_joined, dic_aa_itp, dic_map_ndx, dic_cg_gro,
                internal, external, internal_pairs, external_pairs)
    print(f"✓ Report written to {args.report}")
    
    if args.output_cg_renamed:
        atoms_list = prepare_atoms_list(dic_joined)
        write_gro(args.output_cg_renamed, atoms_list, args.box)
        print(f"✓ CG GRO written to {args.output_cg_renamed}")
    
    print("\n Done!\n")


if __name__ == "__main__":
    main()
