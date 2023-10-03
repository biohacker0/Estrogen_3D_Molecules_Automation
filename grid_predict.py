import os
import subprocess

#This program is to automatically make a grid config file for you which covers the whole receptor, I haven't figured out a algorithmic way to get the perfect coordinate, so lets go the brute force way and cover the whole thing

#This Define the paths to the receptor and output configuration file
original_pdb = r'WRITE_THE_PATH_HERE'  # Replace with your original receptor PDB file
receptor_pdbqt = 'receptor.pdbqt'  # Output receptor PDBQT file name
output_config = 'config.txt'  # Output configuration file name

#WE Use Open Babel to convert the original receptor to PDBQT format
obabel_command = 'obabel'
input_format = '-ipdb'
output_format = '-opdbqt'

#this Run Open Babel to convert the original receptor to PDBQT format
subprocess.run([obabel_command, original_pdb, input_format, '-O', receptor_pdbqt, output_format])

#To Calculate receptor dimensions based on the original PDB file
with open(original_pdb, 'r') as pdb_file:
    min_x = min_y = min_z = float('inf')
    max_x = max_y = max_z = float('-inf')
    
    for line in pdb_file:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            
            min_x = min(min_x, x)
            min_y = min(min_y, y)
            min_z = min(min_z, z)
            max_x = max(max_x, x)
            max_y = max(max_y, y)
            max_z = max(max_z, z)

#this is to Calculate the center and size of the grid box
center_x = (min_x + max_x) / 2.0
center_y = (min_y + max_y) / 2.0
center_z = (min_z + max_z) / 2.0
size_x = max_x - min_x
size_y = max_y - min_y
size_z = max_z - min_z

# Write the AutoDock Vina configuration file
with open(output_config, 'w') as config_file:
    config_file.write(f"receptor = {receptor_pdbqt}\n")
    config_file.write(f"ligand = \n")  # hey You can specify the ligand path here 
    config_file.write(f"out = output.pdbqt\n")
    config_file.write(f"center_x = {center_x}\n")
    config_file.write(f"center_y = {center_y}\n")
    config_file.write(f"center_z = {center_z}\n")
    config_file.write(f"size_x = {size_x}\n")
    config_file.write(f"size_y = {size_y}\n")
    config_file.write(f"size_z = {size_z}\n")
    config_file.write(f"exhaustiveness = 8\n")
    config_file.write(f"num_modes = 9\n")
    config_file.write(f"energy_range = 3\n")

print(f'AutoDock Vina configuration file generated: {output_config}')
