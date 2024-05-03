from Bio.PDB import PDBParser

class AminoAcid:
    def __init__(self, residue_name, chain_id, residue_number, atoms):
        self.residue_name = residue_name
        self.chain_id = chain_id
        self.residue_number = residue_number
        self.atoms = atoms

def distance_squared(aa1, aa2):
    min_dist_sq = float('inf')  # Initialiser la distance minimale entre 2 atomes 
    for atom1 in aa1.atoms:
        for atom2 in aa2.atoms:
            dist_sq = (atom1[0] - atom2[0])**2 + \
                      (atom1[1] - atom2[1])**2 + \
                      (atom1[2] - atom2[2])**2
            if dist_sq < min_dist_sq:
                min_dist_sq = dist_sq  # prendre toujours la valeur la plus petite 
    return min_dist_sq

def detect_interface(aa_list, threshold_distance):
    interface_aa = []
    for i in range(len(aa_list)):
        for j in range(i + 1, len(aa_list)):
            aa_a = aa_list[i]
            aa_b = aa_list[j]
            if aa_a.chain_id != aa_b.chain_id:  # Vérifier s'ils appartiennent à des chaînes différentes
                dist_sq = distance_squared(aa_a, aa_b)
                if dist_sq <= threshold_distance**2:
                    interface_aa.append((aa_a, aa_b))
    return interface_aa

def write_interface_to_file(interface_aa, filename):
    interface_dict = {}  # Dictionnaire pour regrouper les résidus par chaîne
    for aa_pair in interface_aa:
        aa_a, aa_b = aa_pair
        # Ajouter les résidus à la liste correspondante dans le dictionnaire
        if aa_a.chain_id not in interface_dict:
            interface_dict[aa_a.chain_id] = []
        interface_dict[aa_a.chain_id].append((aa_a.residue_name, aa_a.residue_number))
        if aa_b.chain_id not in interface_dict:
            interface_dict[aa_b.chain_id] = []
        interface_dict[aa_b.chain_id].append((aa_b.residue_name, aa_b.residue_number))

    # Trier les résidus par numéro de résidu
    for chain_id, residues in interface_dict.items():
        interface_dict[chain_id] = sorted(set(residues))

    # Écrire les informations dans le fichier
    with open(filename, 'w') as f:
        for chain_id, residues in interface_dict.items():
            f.write(f"Chain {chain_id}:\n")
            for residue_name, residue_number in residues:
                f.write(f"{residue_name}{residue_number}\n")
            f.write("\n")

def parse_pdb(filename):
    parser = PDBParser()
    structure = parser.get_structure("protein", filename)
    amino_acids = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id()[0] == ' ' and residue.get_resname() != 'HOH':  # enlever les molécules d'eau si elles se trouvent dans le fichier 
                    res_name = residue.get_resname()
                    chain_id = chain.get_id()
                    res_number = residue.get_id()[1]
                    atoms = [atom.get_coord() for atom in residue.get_atoms()]
                    amino_acids.append(AminoAcid(res_name, chain_id, res_number, atoms))
    return amino_acids

# Programme principal
pdb_filename = input("Veuillez entrer le nom du fichier PDB : ")
# Saisie utilisateur pour la distance
threshold_distance = float(input("Veuillez entrer la distance en Angströms pour détecter l'interface : "))
all_aa = parse_pdb(pdb_filename)
interface_aa = detect_interface(all_aa, threshold_distance)
write_interface_to_file(interface_aa, f"interface_residues_{threshold_distance}.txt")
