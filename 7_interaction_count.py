import os
import csv
from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch  # For handling PDB and CIF files
from concurrent.futures import ProcessPoolExecutor #for running simultaneously

def classify_chain(chain): #classify RNA, DNA or Protein chains on the first res only (hybrid chains ignored)
    rna_residues = {'A', 'U', 'G', 'C'}
    dna_residues = {'DA', 'DT', 'DG', 'DC'}
    protein_residues = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'}
    for residue in chain: #check the first res and break
        if residue.resname in rna_residues:
            return "RNA"
        elif residue.resname in dna_residues:
            return "DNA"
        elif residue.resname in protein_residues:
            return "Protein"
        break
    return None  #if its not identified

def process_structure(file_path, interaction_threshold=10.0): # IMP: threshold change????
    #PDB/CIF as input. analyze ribose interactions. 
    #tally interactions at each ribose pos
    #RNA/DNA-Protein combinations only.
    print(f"Analyzing file: {file_path}")
    if file_path.endswith(".pdb"):
        parser = PDBParser(QUIET=True)
    elif file_path.endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    else:
        print(f"NOT PDB/CIF, {file_path}")
        return None, "NOT PDB/CIF"
    try:
        structure = parser.get_structure('Structure', file_path)
    except Exception as e:
        return None, f"Parsing error: {e}"

    interaction_counts = {'C1\'': 0, 'C2\'': 0, 'C3\'': 0, 'C4\'': 0, 'O4\'': 0} 
    #how these positions are classified in PDB/CIF file (rmsd_bond.auth_atom_id)

    molecule_types = set()  # RNA, DNA, Protein. No dup. easy comparing
    all_atoms = [] # list of all atoms for spatial indexing
    for model in structure:  # hierarchy of divisions in PDB/CIF parser
        for chain in model:  # chains: RNA, DNA, PRO
            for residue in chain:  # AA, N, DN
                for atom in residue:
                    all_atoms.append(atom)
    neighbor_search = NeighborSearch(all_atoms) #K dimensions trees for spatial indexing (makes grids to handle coordinates efficiently)
    #(O(n\logn) instead of O(n2) with brute force calc)

    for model in structure:
        for chain in model:
            chain_type = classify_chain(chain) # Classify chain from first res
            molecule_types.add(chain_type)
            print(f"{chain_type} chain")
            if chain_type in {"RNA", "DNA"}:  # Only process RNA/DNA chains and consider interactions with proteins
                for residue in chain:
                    if residue.resname in ['DA', 'DT', 'DG', 'DC', 'A', 'U', 'G', 'C']: #only N/dN
                        for atom in residue:
                            if atom.name in interaction_counts:  # (rmsd_bond.auth_atom_id)= ribose atoms
                                close_atoms = neighbor_search.search(atom.coord, interaction_threshold)
                                for other_atom in close_atoms:
                                    # Check if the other atom belongs to a protein residue
                                    parent_residue = other_atom.get_parent()
                                    if parent_residue.resname in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']:
                                        interaction_counts[atom.name] += 1
    return interaction_counts, molecule_types

def process_file(file_path):
    return file_path, process_structure(file_path)

def save_ribose_results_to_csv(results, output_file):#csv output 
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['PDB File', 'Molecule Type', 'C1\'', 'C2\'', 'C3\'', 'C4\'', 'O4\'']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()  # Write header row
        for pdb_file, data in results.items():
            row = {'PDB File': pdb_file, 'Molecule Type': data['type']}
            row.update(data['counts'])  # Add counts to the row
            writer.writerow(row)

if __name__ == "__main__":
    try:
        input_directory = "RCSB_pdb_files"
        print(f"Processing from {input_directory}")
        # all PDB/CIF in folder
        files_to_process = []
        for f in os.listdir(input_directory):
            if f.endswith(".pdb") or f.endswith(".cif"):
                full_path = os.path.join(input_directory, f)
                files_to_process.append(full_path)
        with ProcessPoolExecutor() as executor: #concurrent or parallel processing
            results = list(executor.map(process_file, files_to_process))
        ribose_results = {} #combining result
        for file_path, result in results:
            molecule_counts = result[0]  # Extracting
            molecule_types = result[1] 
            if molecule_types: #only if not empty
                filtered_types = [] #filtering none molecule_type
                for m_type in molecule_types:
                    if m_type:
                        filtered_types.append(m_type)
                ribose_results[file_path] = {'type': ", ".join(filtered_types), 'counts': molecule_counts}
        # Save results to CSV
        save_ribose_results_to_csv(ribose_results, "ribose_interactions.csv")
        print("saved to ribose_interactions.csv")
    except KeyboardInterrupt:
        print("Interrupted by keyboard")
        executor.shutdown(wait=False)
        raise