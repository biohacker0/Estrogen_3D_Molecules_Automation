from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import json
import os
# this program is to generate sterioisomers of the compund and its image and save them in a images folder
# Define the SMILES representation of estradiol
estradiol_smiles = "CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O"

# Convert the SMILES string to an RDKit molecule object
estradiol_molecule = Chem.MolFromSmiles(estradiol_smiles)

# Generate and save isomers
isomers = list(AllChem.EnumerateStereoisomers(estradiol_molecule))
num_isomers = len(isomers)

# Create a directory to save images
os.makedirs("images", exist_ok=True)

isomer_data = []

# Generate images and save isomers to JSON
for i, isomer in enumerate(isomers):
    isomer_smiles = Chem.MolToSmiles(isomer)
    
    # Generate image of the isomer and save it in the "images" folder
    img = Draw.MolToImage(isomer, size=(300, 300))
    img.save(f"images/isomer_{i + 1}.png")
    
    # Collect isomer data
    isomer_data.append({"SMILES": isomer_smiles, "Image": f"isomer_{i + 1}.png"})

# Save isomer data to a JSON file in the "images" folder
with open("images/estradiol_isomers.json", "w") as json_file:
    json.dump(isomer_data, json_file, indent=4)

print(f"Number of isomers generated: {num_isomers}")
print("Isomers saved to images/estradiol_isomers.json and images/ folder.")
