from rdkit import Chem
import json

#This might look stupid , but yes this script if to verify if the isomers we generated are correct and valid.
# Load the JSON file containing isomer data
with open("images/estradiol_isomers.json", "r") as json_file:
    isomer_data = json.load(json_file)

# Initialize a list to store verification results
verification_results = []

# Verify each isomer's SMILES representation
for i, isomer_info in enumerate(isomer_data):
    isomer_smiles = isomer_info["SMILES"]
    
    # Attempt to convert the SMILES back to an RDKit molecule
    molecule = Chem.MolFromSmiles(isomer_smiles)
    
    # Check if the conversion was successful
    if molecule is not None:
        verification_results.append({"Isomer": f"Isomer_{i + 1}", "Verification": "Correct"})
    else:
        verification_results.append({"Isomer": f"Isomer_{i + 1}", "Verification": "Incorrect"})

# Print verification results
for result in verification_results:
    print(f"{result['Isomer']}: {result['Verification']}")
