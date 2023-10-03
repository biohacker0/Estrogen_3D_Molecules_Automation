from rdkit import Chem
from rdkit.Chem import AllChem
import json
import os

# This scrit converts all those smile format isomers that we stored in the json file to 3d molecules and saves all of them in a folder

#WE Create a directory to save the 3D structure files in SDF format
output_dir = "output_3d_structures_sdf"
os.makedirs(output_dir, exist_ok=True)

#This Load the JSON file containing isomer data
with open("images/estradiol_isomers.json", "r") as json_file:
    isomer_data = json.load(json_file)

#Then we Generate 3D structures and save them to SDF files
for i, isomer_info in enumerate(isomer_data):
    isomer_smiles = isomer_info["SMILES"]
    
    # Convert SMILES to an RDKit molecule
    molecule = Chem.MolFromSmiles(isomer_smiles)
    
    if molecule is not None:
        # Add explicit hydrogen atoms
        molecule = Chem.AddHs(molecule)
        
        # Generate 3D coordinates
        AllChem.EmbedMolecule(molecule, AllChem.ETKDG())
        
        # Define the output file path
        output_file = os.path.join(output_dir, f"isomer_{i + 1}.sdf")
        
        # Save the 3D structure to an SDF file
        Chem.MolToMolFile(molecule, output_file)

print(f"3D structures saved in the '{output_dir}' directory in SDF format.")
