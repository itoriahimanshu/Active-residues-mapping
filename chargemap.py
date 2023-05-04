   import nglview as nv

# Load the PDB file
pdb_file = "5fr9.pdb"
view = nv.show_file(pdb_file)

# Select the target amino acid by its ID number and hide everything else
target_id = 128  # Replace with the ID number of the target amino acid
selection = f"resid {target_id}"
view.add_representation("ball+stick", selection=selection)
view.add_representation("surface", selection=selection)

# Add electrostatic potential coloring
view.add_surface(opacity=0.8, surfaceType="av", colorScheme="electrostatic", contour=50)

# Hide the PDB file
view.clear_representations()

# Show the view
view
