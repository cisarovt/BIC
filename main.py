#!/usr/bin/ipython
# -*- coding: utf-8 -*-
# author: Císařová Terezie
"""The Bilayer Initial Constructor (BIC) script constructs 128 molecule phospholipid bilayer from up to 3 molecule 
types. Classes from module classes and exceptions from module errors are imported and used.
 
This script requires that 'os','numpy', 'pandas', 'random', 'scipy', 'datetime', 
'glob' and 'abc' be installed within the Python environment you are running this script in.

Briefly, it takes provided *.pdb files, stripes the header and the footer and uses only the molecule information.
Such information is loaded into pandas dataframe, and the coordinates for bilayer component in both upper and 
lower leaflet are created. Then, the head template is rescaled according to the size of bilayer components to
prevent atom overlaping. The counts of molecules in upper and lower layer are also provided by user. After that, 
to each coordinate in head template is randomly assigned one bilayer component, depending on the upper/lower 
location it is chosen the direction of the bilayer component. Such component is the randomly rotated along the 
Z-axis passing through its head atom. The head atom is also provided by user. Finaly, each bilayer component 
molecule is translated on its new position given by the head templat. Whole bilayer is printed with new header 
and footer to a new *.pdb file.

User input
----------
Names of files with bilayer components.
Acronyms of bilayer components.
Counts of bilayer components in both upper and lower layer.
Head atom of each bilayer component.
Name of the new *.pdb file (if not provided, is created automaticaly).

"""
import os
import errors
import classes

print("This script will generate mixed phospholipid bilayer from given\
\nmolecules (phospholipids or other bilayer components) in given ratio\
\n(by given numbers of each molecule type for each leaflet).")

num_comp = int(input("How many types of molecules do you \
\nwish to incorporate into the bilayer (1 - 3)?:\n"))
if num_comp < 1 or num_comp > 3:
    raise errors.WrongMoleculeTypeNumber

## Data input
# initialisation of variables for loop
molecule_filename = [] #list with names of pdb files of the molecules
number_mol_lower = [] # list with numbers of molecules of given type in the lower layer
number_mol_upper = [] # list with numbers of molecules of given type in the upper layer
name_mol = [] #list with molecule names
num_mol_lower = 0 # initalisation of the overall number of molecules to use on lower layer,
#should be 64 in total
num_mol_upper = 0 # initalisation of the overall number of molecules to use on upper layer,\
#should be 64 in total

for i in range(num_comp):
    while True: # should give user another chance if the input is wrong, does not work
        file_input = str(input(f"Type the name of the pdb file of the {i+1}. component:\n"))
        if os.path.isfile(file_input) is True: # does the file exist?
            break
        else:
            raise errors.FileNotFoundError(file_input = file_input)
    molecule_filename.append(file_input)

    name_mol.append(str(input("Type the acronym of the molecule:\n")))

    while True:
        try:
            number_mol_upper.append(int(input(f"Type the number of molecules of the {i+1}.\
 component in upper layer (int, 0 - {64 - num_mol_upper}): \n")))
        except ValueError:
            raise ValueError
        else:
            break         
    while True:
        try:
            number_mol_lower.append(int(input(f"Type the number of molecules of the {i+1}.\
 component in lower layer (int, 0 - {64 - num_mol_lower}): \n")))
        except ValueError:
            raise ValueError
        else:
            break

    num_mol_lower = num_mol_lower+number_mol_lower[i]
    num_mol_upper = num_mol_upper+number_mol_upper[i]

## raises exception if the overall count is not 128
if num_mol_lower != 64:
    raise errors.WrongMoleculeNumber(number_mol_lower)
if num_mol_upper != 64:
    raise errors.WrongMoleculeNumber(number_mol_upper)

## loading and extracting data from .pdb file of the head template
nheads = classes.NHeads(molname="Nheads")

##initialisation of variables for loop
molecule = [0] * num_comp # list with zeros of length of objects
spans_X = [0] * num_comp
spans_Y = [0] * num_comp
spans_Z =[0] * num_comp
for i in range(num_comp):
    ## extracting the information from *.pdb files of bilayer components
    molecule[i]= classes.BilayerComp(pdb_filename=molecule_filename[i], molname=name_mol[i])
    ## creating list of atom names for the choice of head atom for each bilayer component
    atom_names = set(molecule[i].atom_name.squeeze())
    atom_names = [x for x in atom_names if not x.startswith(("H","C"))] # assuming head atom cannot be hydrogen or carbon
    print(f"All non-carbon and non-hydrogen atom names in {molecule[i].molecule_name} are:")
    for idx, ele in enumerate(atom_names):
        print(f"{idx}: {ele}")

    while True:
        try:
            head_atom = int(input(f"Type the number of the atom representing the head of the molecule {name_mol[i]}\n \
(The  head of the molecule is the atom, which is inserted instead of the nitrogen atom in the phospholipid bilayer.)\n: "))
        except ValueError:
            raise ValueError
        except NameError:
            raise NameError
        except IndexError:
            raise IndexError
        else:
            break
    head_atom = atom_names[head_atom]

    #molecule[i].find_dim(head_atom)

    molecule[i].lowr(head_atom) #creates coordinates for bilayer component in the lower sheet
    molecule[i].uppr(head_atom) #creates coordinates for bilayer component in the upper sheet
    ## extract molecule spans for the scaling of the template
    molecule[i].find_dim()
    spans_X[i] = molecule[i].span_X
    spans_Y[i] = molecule[i].span_Y
    spans_Z[i] = molecule[i].Zspan#.span_Z # taking distance headatom-tail, not max to min

max_span_X = max(spans_X)
max_span_Y = max(spans_Y)
max_span_Z = max(spans_Z)

## rescaling of the nitrogen atoms to prevent atom overlaping
nheads.scale(mol_maxspan_X=max_span_X,mol_maxspan_Y=max_span_Y,mol_maxspan_Z=max_span_Z) 

new_name = str(input("Type in the whole name of your bilayer .pdb:"))

## creation of new *.pdb file
bilayer = classes.NewPDB(pdb_filename=new_name, N_dims=[nheads.new_X,nheads.new_Y,nheads.new_Z],\
                         comp_list=molecule,number_mol_upper=number_mol_upper,number_mol_lower=number_mol_lower)
## print into fwf
bilayer.print_PDB()
  