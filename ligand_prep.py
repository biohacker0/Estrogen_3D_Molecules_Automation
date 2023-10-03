import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

#This script is to prepare the ligand, look at the fucntion names to see what things its doing, if you dont want a partucular thing, just modify the script to scipt that or add something to do something you want in your ligand prep.

#This Define the path to the input SDF file containing isomer molecules
input_sdf_file = '_WRITE_THE_INPUT_FILE_PATH'

#This Define the path to the output SDF file for prepared ligands
output_sdf_file = '_WRITE_THE_OUTPUT_FILE_PATH'

#This Create an SDMolSupplier to read molecules from the input SDF file
suppl = Chem.SDMolSupplier(input_sdf_file)

#WE Create a list to store molecule objects
prepared_ligands = []

#Then we Iterate through the molecules and perform ligand preparation
for molecule in suppl:
    if molecule is not None:
        # Standardize the molecule
        molecule = Chem.AddHs(molecule)
        Chem.SanitizeMol(molecule)

        # Perform 3D geometry optimization using ETKDG
        AllChem.EmbedMolecule(molecule, AllChem.ETKDG())
        AllChem.MMFFOptimizeMolecule(molecule)

        # Compute Gasteiger charges
        AllChem.ComputeGasteigerCharges(molecule)

        # Handle stereochemistry
        Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)

        # Append the prepared ligand to the list
        prepared_ligands.append(molecule)

#After that Write the prepared ligands to the output SDF file
w = Chem.SDWriter(output_sdf_file)
for mol in prepared_ligands:
    w.write(mol)
w.close()

#In end Use Open Babel to convert the SDF file to PDBQT format
output_pdbqt_file = 'WRITE_YOUR_PATH_HERE where you want your output to be with extension'
subprocess.run(["obabel", "-isdf", output_sdf_file, "-opdbqt", "-O", output_pdbqt_file])

print(f'Ligand preparation completed. Prepared ligands saved to {output_pdbqt_file}')
