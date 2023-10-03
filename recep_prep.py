import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

#this script is to automate the receptor preparation step, we have functions to all those things that you do in GUI, if you dont wnat something to be done to your file, just remove that part of code and modify script, also if you want to add something,then too you can add some of code from your side.

# Step 1: Remove Ligand (ESTRADIOL) and Save as New PDB (optional, you can remove this code if this is not of your use case)
def remove_ligand(input_pdb, output_pdb):
    mol_supplier = Chem.SDMolSupplier(input_pdb)
    ligand_atoms = set()
    for mol in mol_supplier:
        if mol is not None:
            for atom in mol.GetAtoms():
                if atom.GetPDBResidueInfo().GetName() == "EST":
                    ligand_atoms.add(atom.GetPDBResidueInfo().GetIdx())

    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        for line in infile:
            if not line.startswith('HETATM') or int(line[6:11]) not in ligand_atoms:
                outfile.write(line)

# Step 2: Add Missing Hydrogen Atoms
def add_hydrogens(input_pdb, output_pdb):
    mol = Chem.MolFromPDBFile(input_pdb)
    mol = Chem.AddHs(mol)
    Chem.MolToPDBFile(mol, output_pdb)

def remove_water(input_pdb, output_pdb):
    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        for line in infile:
            if not line.startswith('HETATM') or 'HOH' not in line[17:20]:
                outfile.write(line)

# Step 3: Convert Receptor PDB to PDBQT Using Open Babel
def convert_to_pdbqt(input_pdb, output_pdbqt):
    obabel_command = 'obabel'
    input_format = '-ipdb'
    output_format = '-opdbqt'
    
    # Run Open Babel subprocess for format conversion
    subprocess.run([obabel_command, input_pdb, input_format, '-O', output_pdbqt, output_format])

if __name__ == "__main__":
    # Define input and output file paths
    input_pdb_file = r'WRITE_THE_INPUT_PATH_OF_PDB_FILE_OR_RECEPTOR'  # Replace with your receptor PDB file path
    output_pdb_file = "apo_receptor.pdb"  # Output apo receptor PDB file name
    output_pdbqt_file = "apo_receptor.pdbqt"  # Output apo receptor PDBQT file name

    # Call the remove_water function before other steps
    remove_water(output_pdb_file, output_pdb_file)

    # Step 1: Remove Ligand (ESTRADIOL)
    remove_ligand(input_pdb_file, output_pdb_file)

    # Step 2: Add Missing Hydrogen Atoms
    add_hydrogens(output_pdb_file, output_pdb_file)


    # Step 3: Convert to PDBQT Using Open Babel
    convert_to_pdbqt(output_pdb_file, output_pdbqt_file)

    print("Receptor preparation completed.")
