# BIC (Bilayer Initial Constructor)
This Python script facilitates the generation of lipid bilayers for molecular simulations through a step-by-step process:

Input Acquisition:

Users provide .pdb files containing molecular structures.
Specify the number of molecules for each molecule type in each leaflet, ensuring a total of 64 molecules in each leaflet.
## File Parsing:

The script reads and extracts relevant information from the .pdb files, focusing on lines containing 'ATOM' or 'HETATM' keywords, while ignoring extraneous content.

## Dimension Computation:

Calculate the dimensions of each molecule within the xyz space to ensure accurate bilayer construction.
## Direction Determination:

Identify the direction of each molecule, essential for proper bilayer orientation.
## Head Atom Selection:

Prompt users to select a reference atom for each molecule to facilitate alignment.
## Bilayer Generation:

Create bilayer coordinates using chosen reference atoms, ensuring optimized spacing between molecules in both leaflets.
## Coordinate Adjustment:

Fit the head atoms of each molecule randomly to scaled coordinates of the reference nitrogen atoms, with the rest of the atoms within each molecule translated accordingly.
## File Output:

Generate a new .pdb file containing the constructed bilayer, including additional header information such as creation time, molecule counts, and simulation details.
This procedure allows users to customize the composition and structure of lipid bilayers for molecular simulations, ensuring accuracy and efficiency in simulation setups.
