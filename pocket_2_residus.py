import os
from Bio.PDB import PDBParser

# Fonction pour extraire les résidus d'un fichier PDB
def extract_residues_from_pdb(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)
    residues = set()
    for model in structure:
        for chain in model:
            for residue in chain:
                residues.add(residue.get_resname() + str(residue.get_id()[1]))
    return residues

# Charger les résidus de l'interface depuis le fichier
with open("interface.txt", "r") as interface_file:
    interface_residues = set(interface_file.read().split())

# Répertoire contenant les fichiers PDB
pdb_directory = "/home/sdv/m1isdd/ksebouai/Bureau/STAGE/fpocket_results/pockets_pdb_files"

# Liste pour stocker les noms des fichiers contenant des résidus de l'interface
files_with_interface_residues = []

# Vérifier quels fichiers PDB contiennent des résidus de l'interface
for filename in os.listdir(pdb_directory):
    if filename.endswith(".pdb"):
        pdb_file = os.path.join(pdb_directory, filename)
        pdb_residues = extract_residues_from_pdb(pdb_file)
        intersection = pdb_residues.intersection(interface_residues)
        # Vérifier si l'intersection contient au moins deux résidus
        if len(intersection) >= 1:
            files_with_interface_residues.append(filename)

# Enregistrer les noms des fichiers dans un fichier texte
output_file = "fichiers_interface_1res.txt"
with open(output_file, "w") as f:
    for filename in files_with_interface_residues:
        f.write(filename + "\n")

print("Les noms des fichiers contenant au moins deux résidus de l'interface ont été enregistrés dans", output_file)

