#!/usr/bin/env python3
import argparse
import numpy as np
from scipy.spatial import cKDTree
from collections import defaultdict

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
    bead_lines = lines[2:-1]  # skip header and box
    for i, line in enumerate(bead_lines):
        bead_id = f"bead_{i+1}"
        resnr   = int(line[0:5].strip())
        resname = line[5:10].strip()
        beadname = line[10:15].strip()
        beadnr  = int(line[15:20].strip())
        x = float(line[20:28].strip())
        y = float(line[28:36].strip())
        z = float(line[36:44].strip())
        dic_cg_gro[bead_id] = (resnr, resname, beadname, beadnr, x, y, z)
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
# GROUP BEADS BY RESIDUE
# =========================
def group_beads_by_residue(dic_joined):
    """Group beads by residue number"""
    residue_beads = defaultdict(list)
    for bead_key, bead_data in dic_joined.items():
        resnr = bead_data["resnr"]
        residue_beads[resnr].append(bead_key)
    return residue_beads


# =========================
# FIND BEST TRIANGLE FROM BEADS
# =========================
def find_best_triangle(beads, bead_positions, aa_edges):
    """
    Find the best triangle from a set of beads (minimum perimeter)
    Returns: tuple of 3 beads forming the triangle
    """
    if len(beads) < 3:
        return None
    
    best_triangle = None
    min_perimeter = float('inf')
    
    for i in range(len(beads)):
        for j in range(i+1, len(beads)):
            for k in range(j+1, len(beads)):
                b1, b2, b3 = beads[i], beads[j], beads[k]
                d12 = np.linalg.norm(np.array(bead_positions[b1]) - np.array(bead_positions[b2]))
                d13 = np.linalg.norm(np.array(bead_positions[b1]) - np.array(bead_positions[b3]))
                d23 = np.linalg.norm(np.array(bead_positions[b2]) - np.array(bead_positions[b3]))
                perimeter = d12 + d13 + d23
                
                if perimeter < min_perimeter:
                    min_perimeter = perimeter
                    best_triangle = (b1, b2, b3)
    
    return best_triangle


# =========================
# FIND INTERNAL CONNECTIONS WITHIN RESIDUE
# =========================
def find_internal_connections(beads, bead_positions, aa_connections):
    """
    Find internal connections within a residue.
    For 4 beads: triangle (3 edges) + 1 branch connection (1 edge) = 4 edges max
    For 3 beads: triangle (3 edges)
    For 2 beads: 1 edge
    """
    if len(beads) == 1:
        return set()
    
    if len(beads) == 2:
        # Single bond
        return {tuple(sorted(beads))}
    
    # Get AA-based connections within this residue
    bead_set = set(beads)
    aa_edges = set()
    for b1, b2 in aa_connections:
        if b1 in bead_set and b2 in bead_set:
            aa_edges.add(tuple(sorted([b1, b2])))
    
    # For 3 beads, we need exactly 3 edges (triangle)
    if len(beads) == 3:
        triangle_edges = set()
        for i in range(3):
            for j in range(i+1, 3):
                triangle_edges.add(tuple(sorted([beads[i], beads[j]])))
        return triangle_edges
    
    # For 4 beads: need exactly 4 edges (triangle + branch)
    if len(beads) == 4:
        # Find the best triangle (minimum perimeter)
        best_triangle = find_best_triangle(beads, bead_positions, aa_edges)
        
        if best_triangle is None:
            return set()
        
        triangle_set = set(best_triangle)
        branch_beads = [b for b in beads if b not in triangle_set]
        
        # Calculate all possible distances to find the optimal structure
        edges = set()
        b1, b2, b3 = best_triangle
        
        # Add all triangle edges
        edges.add(tuple(sorted([b1, b2])))
        edges.add(tuple(sorted([b1, b3])))
        edges.add(tuple(sorted([b2, b3])))
        
        # Connect branch bead to the closest triangle vertex
        for branch in branch_beads:
            # Calculate distances to triangle vertices
            dist_to_b1 = np.linalg.norm(np.array(bead_positions[branch]) - np.array(bead_positions[b1]))
            dist_to_b2 = np.linalg.norm(np.array(bead_positions[branch]) - np.array(bead_positions[b2]))
            dist_to_b3 = np.linalg.norm(np.array(bead_positions[branch]) - np.array(bead_positions[b3]))
            
            # Find closest vertex
            distances = [(dist_to_b1, b1), (dist_to_b2, b2), (dist_to_b3, b3)]
            distances.sort(key=lambda x: x[0])
            
            # Add connection to closest vertex
            closest = distances[0][1]
            edges.add(tuple(sorted([branch, closest])))
        
        # Ensure we have exactly 4 edges (remove the longest if we have more)
        if len(edges) > 4:
            # Convert to list of edges with distances
            edge_list = []
            for edge in edges:
                bead_a, bead_b = edge
                dist = np.linalg.norm(np.array(bead_positions[bead_a]) - 
                                     np.array(bead_positions[bead_b]))
                edge_list.append((dist, edge))
            
            # Sort by distance (shorter is better)
            edge_list.sort(key=lambda x: x[0])
            
            # Keep only the 4 shortest edges
            edges = {edge for dist, edge in edge_list[:4]}
        
        return edges
    
    # For >4 beads, similar logic but keep structure compact
    # Find best triangle first
    best_triangle = find_best_triangle(beads, bead_positions, aa_edges)
    
    if best_triangle is None:
        return set()
    
    triangle_set = set(best_triangle)
    branch_beads = [b for b in beads if b not in triangle_set]
    
    edges = set()
    b1, b2, b3 = best_triangle
    
    # Add triangle edges
    edges.add(tuple(sorted([b1, b2])))
    edges.add(tuple(sorted([b1, b3])))
    edges.add(tuple(sorted([b2, b3])))
    
    # Connect each branch to its closest triangle vertex
    for branch in branch_beads:
        dist_to_b1 = np.linalg.norm(np.array(bead_positions[branch]) - np.array(bead_positions[b1]))
        dist_to_b2 = np.linalg.norm(np.array(bead_positions[branch]) - np.array(bead_positions[b2]))
        dist_to_b3 = np.linalg.norm(np.array(bead_positions[branch]) - np.array(bead_positions[b3]))
        
        distances = [(dist_to_b1, b1), (dist_to_b2, b2), (dist_to_b3, b3)]
        distances.sort(key=lambda x: x[0])
        
        closest = distances[0][1]
        edges.add(tuple(sorted([branch, closest])))
    
    return edges


# =========================
# FIND EXTERNAL CONNECTIONS (between residues)
# =========================
def find_external_connections(dic_joined, aa_connections, max_connections_per_bead=4):
    """
    Find external connections between residues using kNN and AA priorities
    """
    # Get bead positions
    bead_positions = {}
    bead_residues = {}
    for bead_key, bead_data in dic_joined.items():
        bead_positions[bead_key] = (bead_data["x"], bead_data["y"], bead_data["z"])
        bead_residues[bead_key] = bead_data["resnr"]
    
    # Group by residue
    residue_beads = group_beads_by_residue(dic_joined)
    
    # Track connections and their counts
    all_connections = set(aa_connections)
    connection_counts = defaultdict(int)
    for conn in aa_connections:
        connection_counts[conn] += 1
    
    # Track connection counts per bead
    bead_degree = defaultdict(int)
    for b1, b2 in all_connections:
        bead_degree[b1] += 1
        bead_degree[b2] += 1
    
    # Process each residue to ensure proper internal connections (triangulation)
    for resnr, beads in residue_beads.items():
        if len(beads) >= 3:
            # Get internal connections for this residue
            internal_edges = find_internal_connections(beads, bead_positions, all_connections)
            
            # First, remove existing internal connections for this residue
            to_remove = []
            for conn in all_connections:
                if conn[0] in beads and conn[1] in beads:
                    to_remove.append(conn)
            
            for conn in to_remove:
                all_connections.discard(conn)
                # Decrease degree counts
                bead_degree[conn[0]] -= 1
                bead_degree[conn[1]] -= 1
            
            # Add the new internal edges
            for edge in internal_edges:
                if edge not in all_connections:
                    all_connections.add(edge)
                    bead_degree[edge[0]] += 1
                    bead_degree[edge[1]] += 1
    
    # Now find external connections between residues using kNN
    # Create list of all beads with their positions
    bead_list = list(bead_positions.keys())
    positions = np.array([bead_positions[b] for b in bead_list])
    
    # Build kd-tree
    tree = cKDTree(positions)
    
    # For each bead, find potential external connections
    potential_connections = []
    for i, bead1 in enumerate(bead_list):
        # Find k nearest neighbors (k=5 to have options)
        distances, indices = tree.query(positions[i], k=min(6, len(bead_list)))
        
        for j, idx in enumerate(indices[1:], 1):  # skip self
            bead2 = bead_list[idx]
            dist = distances[j]
            
            # Skip if same residue
            if bead_residues[bead1] == bead_residues[bead2]:
                continue
            
            # Skip if connection already exists
            if (bead1, bead2) in all_connections or (bead2, bead1) in all_connections:
                continue
            
            # Check max connections per bead
            if bead_degree[bead1] >= max_connections_per_bead or \
               bead_degree[bead2] >= max_connections_per_bead:
                continue
            
            potential_connections.append((dist, bead1, bead2))
    
    # Sort by distance and add closest connections first
    potential_connections.sort(key=lambda x: x[0])
    
    # We need exactly 8 external connections total
    external_connections_needed = 8
    current_external = len([c for c in all_connections 
                           if bead_residues[c[0]] != bead_residues[c[1]]])
    
    for dist, bead1, bead2 in potential_connections:
        if current_external >= external_connections_needed:
            break
        
        # Check again if connection still valid
        if (bead1, bead2) in all_connections or (bead2, bead1) in all_connections:
            continue
        
        if bead_degree[bead1] >= max_connections_per_bead or \
           bead_degree[bead2] >= max_connections_per_bead:
            continue
        
        all_connections.add((bead1, bead2))
        bead_degree[bead1] += 1
        bead_degree[bead2] += 1
        
        if bead_residues[bead1] != bead_residues[bead2]:
            current_external += 1
    
    return all_connections, connection_counts


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
# WRITE BONDS NDX FILE
# =========================
def write_bonds_ndx(filename, dic_joined, internal_pairs, external_pairs):
    """
    Write bonds.ndx file with internal and external bonds
    Format: bead_index1 bead_index2 ; residue1_bead---residue2_bead dist=value (in Angstroms)
    """
    # First, create a mapping from bead_key to bead index (1-based as in cg.gro)
    bead_to_index = {}
    for i, bead_key in enumerate(dic_joined.keys(), start=1):
        bead_to_index[bead_key] = i
    
    # Also store bead names and residues for each index
    index_to_info = {}
    for bead_key, bead_data in dic_joined.items():
        idx = bead_to_index[bead_key]
        index_to_info[idx] = {
            "resnr": bead_data["resnr"],
            "bead_name": bead_data["bead_name"],
            "cg_name": bead_data["cg_name"]
        }
    
    # Calculate distances for each bond (convert nm to Angstroms: multiply by 10)
    def calculate_distance(bead1_key, bead2_key, dic_joined):
        """Calculate Euclidean distance between two beads in Angstroms"""
        x1, y1, z1 = dic_joined[bead1_key]["x"], dic_joined[bead1_key]["y"], dic_joined[bead1_key]["z"]
        x2, y2, z2 = dic_joined[bead2_key]["x"], dic_joined[bead2_key]["y"], dic_joined[bead2_key]["z"]
        # Convert from nm to Angstroms (1 nm = 10 Angstroms)
        dist_nm = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
        dist_angstrom = dist_nm * 10.0
        return dist_angstrom
    
    # Sort internal bonds by residue number and then by index
    def get_sort_key(bond_pair):
        bead1_key, bead2_key = bond_pair
        idx1 = bead_to_index[bead1_key]
        idx2 = bead_to_index[bead2_key]
        res1 = index_to_info[idx1]["resnr"]
        res2 = index_to_info[idx2]["resnr"]
        # Sort primarily by residue number, then by index
        return (min(res1, res2), min(idx1, idx2))
    
    sorted_internal = sorted(internal_pairs, key=get_sort_key)
    sorted_external = sorted(external_pairs, key=get_sort_key)
    
    with open(filename, 'w') as f:
        f.write("[ bonds ]\n")
        
        # Write internal bonds
        f.write(";;;;;;; internal bonds\n\n")
        for bead1_key, bead2_key in sorted_internal:
            idx1 = bead_to_index[bead1_key]
            idx2 = bead_to_index[bead2_key]
            dist = calculate_distance(bead1_key, bead2_key, dic_joined)
            
            # Get residue and bead info
            res1 = index_to_info[idx1]["resnr"]
            res2 = index_to_info[idx2]["resnr"]
            cg_name1 = index_to_info[idx1]["cg_name"]
            cg_name2 = index_to_info[idx2]["cg_name"]
            
            f.write(f"{idx1:6d}{idx2:6d} ; {res1}_{cg_name1}---{res2}_{cg_name2} dist={dist:.3f}\n")
        
        # Write external bonds
        f.write("\n;;;;;;; external bonds\n\n")
        for bead1_key, bead2_key in sorted_external:
            idx1 = bead_to_index[bead1_key]
            idx2 = bead_to_index[bead2_key]
            dist = calculate_distance(bead1_key, bead2_key, dic_joined)
            
            # Get residue and bead info
            res1 = index_to_info[idx1]["resnr"]
            res2 = index_to_info[idx2]["resnr"]
            cg_name1 = index_to_info[idx1]["cg_name"]
            cg_name2 = index_to_info[idx2]["cg_name"]
            
            f.write(f"{idx1:6d}{idx2:6d} ; {res1}_{cg_name1}---{res2}_{cg_name2} dist={dist:.3f}\n")


# =========================
# WRITE REPORT
# =========================
def write_report(filename, bead_connections, connection_counts, dic_joined, 
                internal, external, internal_pairs, external_pairs):
    """Write detailed connectivity report"""
    with open(filename, 'w') as f:
        f.write("=== Bead Connectivity Report ===\n\n")
        f.write("Connections:\n")
        
        # Sort connections for consistent output
        sorted_conns = sorted(bead_connections)
        for b1, b2 in sorted_conns:
            name1 = dic_joined[b1]["bead_name"]
            name2 = dic_joined[b2]["bead_name"]
            f.write(f"{b1}({name1}) ----- {b2}({name2})\n")
        
        f.write("\n=== Statistics ===\n\n")
        f.write(f"Total unique bead connections: {len(bead_connections)}\n")
        f.write(f"Internal bonds (within residue): {internal}\n")
        f.write(f"External bonds (between residues): {external}\n")
        f.write(f"Total AA bonds contributing: {sum(connection_counts.values())}\n")
        
        f.write("\nDetailed counts per connection:\n")
        for conn in sorted_conns:
            count = connection_counts.get(conn, 1)
            name1 = dic_joined[conn[0]]["bead_name"]
            name2 = dic_joined[conn[1]]["bead_name"]
            f.write(f"{conn[0]}({name1}) - {conn[1]}({name2}): {count} AA bonds\n")
        
        f.write("\n=== Internal Connections ===\n")
        for b1, b2 in internal_pairs:
            name1 = dic_joined[b1]["bead_name"]
            name2 = dic_joined[b2]["bead_name"]
            f.write(f"{b1}({name1}) - {b2}({name2})\n")
        
        f.write("\n=== External Connections ===\n")
        for b1, b2 in external_pairs:
            name1 = dic_joined[b1]["bead_name"]
            name2 = dic_joined[b2]["bead_name"]
            f.write(f"{b1}({name1}) - {b2}({name2})\n")


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
            d["resnr"],    # 0
            d["bead_name"],# 1
            d["cg_name"],  # 2
            i,             # 3: bead index
            d["x"],        # 4
            d["y"],        # 5
            d["z"],        # 6
            i              # 7: renum index for GRO
        ))
    return atoms_list


# =========================
# MAIN
# =========================
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_aa_itp", required=True, help="AA topology file (.itp)")
    parser.add_argument("--input_aa_gro", required=True, help="AA coordinate file (.gro)")
    parser.add_argument("--input_cg_ndx", required=True, help="CG index file (.ndx)")
    parser.add_argument("--input_cg_gro", required=True, help="CG coordinate file (.gro)")
    parser.add_argument("--output_cg_renamed", default=None, help="Output CG GRO file (optional)")
    parser.add_argument("--output_bonds_ndx", default="bonds.ndx", help="Output bonds index file")
    parser.add_argument("--box", nargs=3, type=float, default=[0.0,0.0,0.0], help="Box dimensions")
    parser.add_argument("--report", default="connectivity_report.txt", help="Output report file")
    args = parser.parse_args()

    # Parse all inputs
    dic_aa_itp = parse_itp(args.input_aa_itp)
    dic_aa_gro = parse_gro(args.input_aa_gro)
    dic_map_ndx = parse_ndx(args.input_cg_ndx)
    dic_cg_gro = parse_cg_gro(args.input_cg_gro)

    # Join CG
    dic_joined = join_cg_data(dic_map_ndx, dic_cg_gro)

    # Find AA-based bead connections
    aa_connections, aa_connection_counts = find_bead_connections_aa(dic_aa_itp, dic_map_ndx)
    
    # Find all connections (AA + geometric)
    all_connections, connection_counts = find_external_connections(dic_joined, aa_connections)
    
    # Update counts for AA connections
    for conn in aa_connections:
        connection_counts[conn] = aa_connection_counts[conn]
    
    # Count internal and external
    internal, external, internal_pairs, external_pairs = count_internal_external(
        dic_joined, all_connections
    )
    
    print(f"Internal bonds: {internal}")
    print(f"External bonds: {external}")
    print(f"Total bonds: {len(all_connections)}")
    
    # Write bonds.ndx file
    write_bonds_ndx(args.output_bonds_ndx, dic_joined, internal_pairs, external_pairs)
    print(f"Bonds NDX written to {args.output_bonds_ndx}")
    
    # Write report
    write_report(args.report, all_connections, connection_counts, dic_joined,
                internal, external, internal_pairs, external_pairs)
    print(f"Report written to {args.report}")

    # Prepare atoms list and write GRO (only if output_cg_renamed is provided)
    if args.output_cg_renamed:
        atoms_list = prepare_atoms_list(dic_joined)
        write_gro(args.output_cg_renamed, atoms_list, args.box)
        print(f"CG GRO written to {args.output_cg_renamed}")
    else:
        print("No CG GRO output requested (--output_cg_renamed not provided)")


if __name__ == "__main__":
    main()
