import os
import csv
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch
from concurrent.futures import ProcessPoolExecutor

charged_residues = { #charged atom of charged residue present in AA and Nucleotides
    'ASP': ['OD1', 'OD2'],  # -COO⁻
    'GLU': ['OE1', 'OE2'],  # -COO⁻
    'LYS': ['NZ'],          # -NH₃⁺
    'ARG': ['NH1', 'NH2'], 
    'HIS': ['NE2', 'ND1']
}
ribose_dipoles = {'C1\'', 'C2\'', 'C3\'', 'C4\'', 'O4\''}  
rna_residues = {'A', 'U', 'G', 'C'}
dna_residues = {'DA', 'DT', 'DG', 'DC'}
#20 Standard AAs
protein_residues = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'}

def calculate_angle(atom1, atom2, atom3): #Calculating angle between three atoms (degrees)
    vector1 = atom1.coord - atom2.coord
    vector2 = atom3.coord - atom2.coord
    cosine_angle = np.dot(vector1, vector2) / (np.linalg.norm(vector1) * np.linalg.norm(vector2))
    angle = np.degrees(np.arccos(cosine_angle))
    return angle

def classify_chain(chain, num_residues=3):#Classifies a chain based on the first `num_residues` residues. to account for non-standard first residues
    molecule_types = set()
    for i, residue in enumerate(chain):
        if residue.resname in rna_residues:
            molecule_types.add("RNA")
        elif residue.resname in dna_residues:
            molecule_types.add("DNA")
        elif residue.resname in protein_residues:
            molecule_types.add("Protein")
        if i + 1 == num_residues:
            break
    return ", ".join(molecule_types) if molecule_types else "Unknown" #unknown to check later in case theres any

def is_hbond(donor, hydrogen, acceptor): #Checks if an interaction is a HBond based on distance and angle on 4 strength group (cutoff subject to change)
    donor_acceptor_dist = np.linalg.norm(donor.coord - acceptor.coord)
    hydrogen_acceptor_dist = np.linalg.norm(hydrogen.coord - acceptor.coord)
    angle = calculate_angle(donor, hydrogen, acceptor)
    
    if donor_acceptor_dist > 3.9 or hydrogen_acceptor_dist > 3.0 or angle < 120:
        return None
    if donor_acceptor_dist <= 2.8 and hydrogen_acceptor_dist <= 2.0:
        return "Very Strong H-bond"
    elif donor_acceptor_dist <= 3.2 and hydrogen_acceptor_dist <= 2.2:
        return "Strong H-bond"
    elif donor_acceptor_dist <= 3.5 and hydrogen_acceptor_dist <= 2.5:
        return "Moderate H-bond"
    elif donor_acceptor_dist <= 3.9 and hydrogen_acceptor_dist <= 3.0:
        return "Weak H-bond"
    return None

def is_ion_dipole(charged_atom, dipole_atom): #Checks if an interaction is an ion-dipole interaction and classifies strength. (cutoff subject to change)
    distance = np.linalg.norm(charged_atom.coord - dipole_atom.coord)
    if distance <= 4.0:
        return "Strong Ion-Dipole"
    elif distance <= 5.0:
        return "Weak Ion-Dipole"
    return None

def classify_interaction(donor=None, hydrogen=None, acceptor=None, charged_atom=None, dipole_atom=None): #Classifies all interaction types (currently Hbonds and Ion-Dipole)
    interactions = []
    if donor and hydrogen and acceptor:
        hbond_strength = is_hbond(donor, hydrogen, acceptor)
        if hbond_strength:
            interactions.append(hbond_strength)
    if charged_atom and dipole_atom:
        ion_dipole_strength = is_ion_dipole(charged_atom, dipole_atom)
        if ion_dipole_strength:
            interactions.append(ion_dipole_strength)
    return interactions  # List of interaction types

def save_hbonds_to_csv(hbond_data, output_file):
    fieldnames = ['donor', 'donor_residue', 'donor_resi', 'donor_chain',
                  'acceptor', 'acceptor_residue', 'acceptor_resi', 'acceptor_chain',
                  'distance_DA', 'distance_HA', 'angle', 'strength']
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for entry in hbond_data:
            writer.writerow(entry)

def save_ion_dipole_to_csv(ion_dipole_data, output_file):
    fieldnames = ['PDB_ID', 'Chain', 'Charged Atom', 'Charged Residue', 
                  'Charged Residue Index', 'Dipole Atom', 'Dipole Residue Index', 
                  'Distance', 'Interaction Type']
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for entry in ion_dipole_data:
            writer.writerow(entry)

def save_results_to_csv(results, output_file): # for number of each interaction types
    fieldnames = ['PDB File', 'Molecule Type'] + list(results[next(iter(results))]['counts'].keys())
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for pdb_file, data in results.items():
            row = {'PDB File': pdb_file, 'Molecule Type': data['type']}
            row.update(data['counts'])
            writer.writerow(row)

def process_structure(file_path): #main func
    print(f"Analyzing file: {file_path}")
    if file_path.endswith(".pdb"):
        parser = PDBParser(QUIET=True)
    elif file_path.endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    else:
        return {}, set(), []  # empty structures if file is not valid
    try:
        structure = parser.get_structure('Structure', file_path)
    except Exception as e:
        print(f"Error parsing {file_path}: {e}")
        return {}, set(), []  # Return empty results in case of a parsing error
    
    interaction_counts = {
        "Very Strong H-bond": 0, "Strong H-bond": 0, "Moderate H-bond": 0,
        "Weak H-bond": 0, "Strong Ion-Dipole": 0, "Weak Ion-Dipole": 0
    }
    interaction_data = []
    molecule_types = set()
    all_atoms = []
    atom_to_chain = {}

    for model in structure:
        for chain in model:
            chain_type = classify_chain(chain, num_residues=3)
            molecule_types.add(chain_type)
            for residue in chain:
                for atom in residue:
                    all_atoms.append(atom)
                    atom_to_chain[atom] = chain_type

    neighbor_search = NeighborSearch(all_atoms)
    
    for atom in all_atoms:
        if atom.name == "O4'":
            close_atoms = [a for a in neighbor_search.search(atom.coord, 3.9) if a.element in {'N', 'O'}]
            for donor in close_atoms:
                if donor == atom:
                    continue
                interacting_chain_type = atom_to_chain.get(donor, None)
                if interacting_chain_type in {"RNA", "DNA", "Protein"}:
                    for hydrogen in donor.parent:
                        if hydrogen.element == 'H':
                            interaction_types = classify_interaction(donor, hydrogen, atom)
                            for interaction in interaction_types:
                                interaction_counts[interaction] += 1
                                interaction_data.append({
                                    'donor': donor.name,
                                    'donor_residue': donor.get_parent().resname,
                                    'donor_resi': donor.get_parent().id[1],
                                    'donor_chain': donor.get_parent().get_parent().id,
                                    'acceptor': atom.name,
                                    'acceptor_residue': atom.get_parent().resname,
                                    'acceptor_resi': atom.get_parent().id[1],
                                    'acceptor_chain': atom.get_parent().get_parent().id,
                                    'distance_DA': np.linalg.norm(donor.coord - atom.coord),
                                    'distance_HA': np.linalg.norm(hydrogen.coord - atom.coord),
                                    'angle': calculate_angle(donor, hydrogen, atom),
                                    'strength': interaction
                                })
        # Ion-Dipole Interaction Detection
        resname = atom.parent.resname
        if (resname in charged_residues and atom.name in charged_residues[resname]) or \
        ((resname in rna_residues or resname in dna_residues) and atom.name in {'OP1', 'OP2'}):
            close_atoms = [a for a in neighbor_search.search(atom.coord, 6.0) if a.name in ribose_dipoles]
            for dipole_atom in close_atoms:
                interaction = classify_interaction(None, None, None, atom, dipole_atom)
                for ion_dipole in interaction:
                    if atom.element=='H':
                        print(f"DEBUG: Element={atom.element}, Charged Atom={atom.name}, Residue={resname}, Residue Index={atom.get_parent().id[1]}")
                    interaction_counts[ion_dipole] += 1
                    interaction_data.append({
                        'PDB_ID': file_path.split('/')[-1],  # Add PDB filename for multiple structures
                        'Chain': atom.get_parent().get_parent().id,  # Get chain ID
                        'Charged Atom': atom.name,
                        'Charged Residue': resname,
                        'Charged Residue Index': atom.get_parent().id[1],  # Get residue number
                        'Dipole Atom': dipole_atom.name,
                        'Dipole Residue Index': dipole_atom.get_parent().id[1],  # Get residue number
                        'Distance': np.linalg.norm(atom.coord - dipole_atom.coord),
                        'Interaction Type': ion_dipole
                    })
    print(f"Finished analyzing: {file_path}")
    print(f"  Found {len(molecule_types)} molecule types: {molecule_types}")
    print(f"  Interaction counts: {interaction_counts}")

    return interaction_counts, molecule_types, interaction_data

if __name__ == "__main__":
    try:
        input_directory = "PDB_test"
        files_to_process = [os.path.join(input_directory, f) for f in os.listdir(input_directory) if f.endswith(".pdb") or f.endswith(".cif")]
        ribose_results = {}
        hbond_results = []
        ion_dipole_results = []
        
        with ProcessPoolExecutor() as executor:
            results = list(executor.map(process_structure, files_to_process))
        
        for file_path, result in zip(files_to_process, results):
            molecule_counts, molecule_types, interactions = result
            ribose_results[file_path] = {'type': ", ".join(molecule_types), 'counts': molecule_counts}
            for interaction in interactions:
                if 'donor' in interaction:
                    hbond_results.append(interaction)
                elif 'Charged Atom' in interaction:
                    ion_dipole_results.append(interaction)
        save_results_to_csv(ribose_results, "ribose_interactions.csv")
        save_hbonds_to_csv(hbond_results, "hydrogen_bonds.csv")
        save_ion_dipole_to_csv(ion_dipole_results, "ion_dipole.csv")

        print("Saved to ribose_interactions.csv, hydrogen_bonds.csv, and ion_dipole.csv")
    except KeyboardInterrupt:
        print("Interrupted by keyboard")


