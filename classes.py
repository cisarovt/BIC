#!/usr/bin/ipython
# -*- coding: utf-8 -*-
# author: Císařová Terezie
"""This module contains classes and their functions used to construct phospholipid bilayer with 128 molecules 
from up to 3 molecule types. Module abstract_classes is loaded and used for class creation. Imports 'os',
'numpy','random', 'scipy.spacial.transform','datetime' and 'abstract_classes' modules.

Classes
-------
NHeads(pdb_filename="N_heads.pdb", molname = "nheads")
    constructs class handling head template
BilayerComp(pdb_filename: str, molname: str)
    constructs class handing bilayer component data
NewPDB(pdb_filename: str, N_dims, comp_list, number_mol_upper, number_mol_lower, molname="bilayer")
    constructs class handling bilayer creation and file print
"""

import os
import numpy as np
import pandas as pd
import random as rdn
from scipy.spatial.transform import Rotation
from datetime import datetime
from abstract_classes import PDBfile

class NHeads(PDBfile):
    """Class extracting the coordinates of bilayer component 'heads' and scaling them for bilayer creation. 
    """
    object_type = "nheads"
    def __init__(self,pdb_filename="N_heads.pdb", molname = "nheads"):
        """NHeads class constructor. Extracts the coordinate information.

        Parameters
        ----------
        pdb_filename : str, optional
            PDB file with head coordinates, by default "N_heads.pdb".
        molname : str, optional
            Name of atom heads, by default "nheads".
        """
        super().__init__(pdb_filename, molname)

        self.filename = pdb_filename
        self.molecule_name = molname
        new_first_line = "_ATOM_ATNUM|ATNA@RESN$RESQ€|||__dimX____dimY____dimZ__OCCUPYTMPFAC||||||SGIDEL+-"
            #line to create header = column names for the pd dataframe

        with open(pdb_filename, 'r',encoding='utf-8') as original: data = original.read()
        with open("temporary.pdb", 'w',encoding='utf-8') as modified: modified.write(new_first_line+"\n" + data)
        self.data = pd.read_fwf(filepath_or_buffer="temporary.pdb",\
                                colspecs=[(30,38),(38,46),(46,54)],header = 0)

        os.remove("temporary.pdb")
        ## extracting the information from the *.pdb file
        # the XYZ coordinates of the template are the only intersting information from the *.pdb file
        self.X = self.data["__dimX__"]
        self.Y = self.data["__dimY__"]
        self.Z = self.data["__dimZ__"]
        
        ## calculating the maximum space provided for one layer - needed for further scaling
        self.max_span_Z = (max(self.Z) - min(self.Z))/2 
        self.avgZ_upper = np.average(self.Z[self.Z > 0])
        self.avgZ_lower = np.average(self.Z[self.Z < 0])
        #print(self.avgZ_lower, self.avgZ_upper)

    def scale(self, mol_maxspan_X, mol_maxspan_Y, mol_maxspan_Z):
        """Function scaling the head coordinates accordingly to the size of the components.
        It moves the XY zero plane to the middle of the bilayer.

        Parameters
        ----------
        mol_maxspan_X : float
            Maximum width of bilayer components in X-axis.
        mol_maxspan_Y : float
            Maximum width of bilayer components in Y-axis.
        mol_maxspan_Z : float
            Maximum length of bilayer components.
        """
        self.span_X = 6.965 #the smallest span between two heads in head template in X-axis
        self.span_Y = 5.265 #the smallest span between two heads in head template in Y-axis
        #xy_span_N_neighbours = 8.731062363767652

        comp_spanX = float(mol_maxspan_X)
        comp_spanY = float(mol_maxspan_Y)
        comp_spanZ = float(mol_maxspan_Z)

        diagonal = np.sqrt(comp_spanX**2 + comp_spanY**2) * 0.7 #0.9

        ## computing the scaling factors
        #scale_Z = self.max_span_Z/comp_spanZ/3 #3
        scale_X = self.span_X/diagonal
        scale_Y = self.span_Y/diagonal

        ##new coordinates
        #self.new_Z = self.Z / scale_Z \ #1.4, meli by bejt na sobe tesne ted
        self.new_Z = [value - self.avgZ_upper + comp_spanZ * 1  if value > 0 \
                    else value - self.avgZ_lower - comp_spanZ * 1  if value < 0 \
                    else value for value in self.Z]
        #print(self.new_Z)
        self.new_Z = pd.Series(self.new_Z)
        self.new_X = self.X / scale_X
        self.new_Y = self.Y / scale_Y

        ## shift of Z axis - XY zero plane is now in the middle of the bilayer
        subtr = max(self.new_Z) - (max(self.new_Z) - min(self.new_Z))/2
        self.new_Z = self.new_Z.subtract(subtr) # well its already in the middle but whatever

class BilayerComp(PDBfile):
    """Class managing information from *.prb files of bilayer components.
    """
    def __init__(self, pdb_filename: str, molname: str):
        """Extracts the information from *.pdb files of bilayer components. Skips the header and the footer. 
        Calculates the widths of the molecule in XY-axes and its length in Z-axis.

        Parameters
        ----------
        pdb_filename : str
            Name of PDB file of the bilayer component.
        molname : str
            Name (acronym) of the bilayer component.
        """
        super().__init__(pdb_filename, molname)
        self.filename = pdb_filename
        self.molecule_name = molname

        ## skipping header and footer of the *.pdb file
        with open(pdb_filename, 'r',encoding='utf-8') as original:  
            lines = original.readlines()
            line_idxs = []
            for row in lines:
                if row.find("ATOM") or row.find("HETATM") != -1:
                    line_idxs.append(lines.index(row) + 1)

            if len(line_idxs) == 0: # no header nor footer
                original.seek(0) 
                data = original.read()        
            elif line_idxs[0] != 1 and line_idxs[-1] > 50 : # only footer
                last_line_idx = min(line_idxs) -1
                original.seek(0)
                data = original.readlines()[:last_line_idx]
            elif line_idxs[0] == 1 and line_idxs[-1] < 30 : # only header
                first_line_idx = max(line_idxs)
                original.seek(0)
                data = original.readlines()[first_line_idx:]
            else:                                           # header and footer
                for i in range(len(line_idxs)):
                    if line_idxs[i] + 1 != line_idxs[i+1]:
                        span_after_idx = i
                        break

                first = line_idxs[:span_after_idx+1]
                last = line_idxs[span_after_idx+1:]
                first_line_idx = max(first)
                last_line_idx = min(last) - 1
                original.seek(0)
                data = original.readlines()[first_line_idx:last_line_idx]

            string = ""
            for item in data:   #transforming the datatype
                string = string + item
            data = string

        new_first_line = "_ATOM_ATNUM|ATNA@RESN$RESQ€|||__dimX____dimY____dimZ__OCCUPYTMPFAC||||||SGIDEL+-"
            #line to create header = column names for the pd dataframe
        with open("temporary.pdb", 'w',encoding='utf-8') as modified: modified.write(new_first_line+"\n" + data)
        self.data = pd.read_fwf(filepath_or_buffer="temporary.pdb",\
                                colspecs=[(0,6),(6,11),(12,16),(16,17),(17,21),(21,22),(22,26),(26,27),(30,38),(38,46),(46,54),(54,60),(60,66),(72,76),(76,78),(78,80)],header = 0)
        os.remove("temporary.pdb")

        ## extracting single pieces of information from *.pdb file
        self.ATOM = self.data["_ATOM_"]
        self.atom_num = self.data["ATNUM"]
        #here one space
        self.atom_name = self.data["ATNA"]
        self.alt_loc = self.data["@"]
        self.res_name = self.data["RESN"]
        self.chain_identif = self.data["$"]
        self.res_seq_num = self.data["RESQ"]
        self.code_insert_res = self.data["€"]
        #here 3 spaces
        self.X = self.data["__dimX__"]
        self.Y = self.data["__dimY__"]
        self.Z = self.data["__dimZ__"]
        self.occupancy = self.data["OCCUPY"]
        self.temp_fac= self.data["TMPFAC"]
        #here 6 spaces
        self.seg_id = self.data["SGID"]
        self.ele_symb = self.data["EL"]
        self.charge = self.data["+-"]
        self.molname = pd.Series(molname, index = range(len(self.data)))

        ## calculating molecule span in all dimensions for scaling 
        self.span_X = np.max(self.X) - np.min(self.X)
        self.span_Y = np.max(self.Y) - np.min(self.Y)
        self.span_Z = np.max(self.Z) - np.min(self.Z)

        print(f"{self.molecule_name} -- X:{self.span_X}, Y: {self.span_Y}, Z: {self.span_Z}")

    def find_dim(self):#, head_atom = "N"):
        #self.head_atom = head_atom
        #self.head_idx = self.atom_name.index[self.atom_name == str(head_atom)].tolist()[0] 
        self.Zspan = np.max([np.abs(np.max(self.Z) - self.Z.iloc[self.head_idx]), np.abs(np.min(self.Z) - self.Z.iloc[self.head_idx])])


    def rescale(self, head_atom = "N"):
        """Moving the head of the bilayer component to coordinates [0; 0; 0].
        Observing the direction of the bilayer component in provided *.pdb file.

        Parameters
        ----------
        head_atom : str, optional
            Name of the atom, which later replaces the heads in the head template, by default "N".
            This atom is moved to [0; 0; 0]

        Returns
        -------
        direction: str
            Direction of the bilayer component, either "upper" or "lower".
        """
        self.head_atom = head_atom
        self.head_idx = self.atom_name.index[self.atom_name == str(head_atom)].tolist()[0] #finding the index of the head atom
        head_X = self.X.iloc[self.head_idx]
        head_Y = self.Y.iloc[self.head_idx]
        head_Z = self.Z.iloc[self.head_idx]

        ## transporting the head atom to [0,0,0]
        self.new_Z = self.Z.subtract(head_Z)
        self.new_X = self.X.subtract(head_X)
        self.new_Y = self.Y.subtract(head_Y)

        ## observing the direction of the bilayer component
        mean = np.mean(self.new_Z)
        if mean > 0:
            direction = "lower"
        elif mean < 0:
            direction = "upper"
        return direction

    def uppr(self, head_atom ="N"):
        """Creates or copies the coordinates for the bilayer component for the upper leaflet of the bilayer.

        Parameters
        ----------
        head_atom : str, optional
            Name of the atom, which replaces the heads in the head template, by default "N".
        """
        ## gets the direction of the bilayer component in provided *.pdb file
        direction = self.rescale(head_atom=head_atom) 
 
        if direction == "upper":    
            self.upper_Z = self.new_Z
        else: ## if the direction is not upper, it creates the coordinates for the upper leaflet by mirroring the atoms 
            ## through the XY zero plane = inverting the Z-coordinates
            self.upper_Z = self.new_Z.multiply(-1)
        ## XY stays unchanged
        self.upper_X = self.new_X
        self.upper_Y = self.new_Y

    def lowr(self, head_atom = "N"):
        """Creates or copies the coordinates for the bilayer component for the lower leaflet of the bilayer.

        Parameters
        ----------
        head_atom : str, optional
            Name of the atom, which replaces the heads in the head template, by default "N".
        """
        direction = self.rescale(head_atom=head_atom)
        if direction == "lower":
            self.lower_Z = self.new_Z
        else: ## if the direction is not upper, it creates the coordinates for the lower leaflet by inverting the Z-coordinates 
             self.lower_Z = self.new_Z.multiply(-1)
        ## XY stays unchanged
        self.lower_Y = self.new_Y
        self.lower_X = self.new_X

    def random_rotate(self, direction):
        """Rotates randomly the bilayer component around the Z-axis passing through its head atom.

        Parameters
        ----------
        direction : str
            Direction of the bilayer component.

        Returns
        -------
        rotated_X, rotated_Y, rotated_Z : float
            Rotated coordinates.
        """
        theta = rdn.uniform(0,2*np.pi) # random angle
        ## creating rotation matrix
        r = Rotation.from_matrix([[np.cos(theta), -np.sin(theta), 0],
                   [np.sin(theta), np.cos(theta), 0],
                   [0, 0, 1]])
        
        if direction == "upper":
            ## obtaining correct coordinates in needed datatype
            v = [self.upper_X.tolist(),self.upper_Y.tolist(),self.upper_Z.tolist()] 
            v = np.array(v)
        if direction == "lower":
            ## obtaining correct coordinates in needed datatype
            v = [self.lower_X.tolist(),self.lower_Y.tolist(),self.lower_Z.tolist()]
            v = np.array(v)

        v = v.T #needs to be transposed

        rotated_X = []
        rotated_Y = []
        rotated_Z = []
        for i in range(len(self.upper_X)):
            ## aplying rotation matrix
            rotated_XYZ = r.apply(v[i]) 
            rotated_X.append(rotated_XYZ[0])
            rotated_Y.append(rotated_XYZ[1])
            rotated_Z.append(rotated_XYZ[2]) # Z should stay the same, but just to be safe
        return rotated_X, rotated_Y, rotated_Z

class NewPDB(PDBfile): 
    """Class managing the information of the new *.pdb file.
    """
    object_type = "newPDB"
    def __init__(self, pdb_filename: str, N_dims, comp_list, number_mol_upper, number_mol_lower, molname="bilayer"):
        """Creates automaticaly the name of the new *.pdb file if not provided. 
        Constructs the bilayer by copying information from bilayer components and rotating and translating 
        the coordinates of each single molecule.

        Parameters
        ----------
        pdb_filename : str
            Name of the new .*.pdb file of the phospholipid bilayer. 
            If empty, it will be created automaticaly, containing the counts of the phospholipids 
            in each leaflet and their acronyms.
        molname : str
            Name of the bilayer, by default "bilayer".
        N_dims : list
            Coordinates of the head template.
        comp_list : list
            List containing the constructed classes of all bilayer components.
        number_mol_upper : list
            List of molecule counts for the upper leaflet.
        number_mol_lower : list
            List of molecule counts for the lower leaflet.
        """
        super().__init__(pdb_filename, molname)
        self.components = comp_list
        self.num_mol_up = number_mol_upper
        self.num_mol_low = number_mol_lower
        ## constructing the name of the new *.pdb file if not provided.
        if pdb_filename == "":
            name = ""
            for idx, item in enumerate(comp_list):
                name = name + str(item.molecule_name) + "_"+ str(number_mol_lower[idx])+"low" + "_" + str(number_mol_upper[idx]) + "up" + "_"
            self.filename = name[:-1] + ".pdb"
        else:
            self.filename = pdb_filename
        print(self.filename, "\n")

        ## the dictionary creation with counts and types of molecules
        ## constructed for easier handling with the information
        mol_dic = {}
        molecule = [0]*len(comp_list)

        layer_upper = []
        layer_lower = []

        N_X, N_Y, N_Z = N_dims

        for idx, item in enumerate(comp_list):
            molecule[idx] = item 
            mol_dic[str(molecule[idx].molecule_name)] = molecule[idx]
            layer_upper = layer_upper + [str(molecule[idx].molecule_name)]*number_mol_upper[idx]
            layer_lower = layer_lower + [str(molecule[idx].molecule_name)]*number_mol_lower[idx]

        ## containers for information initialisation
        self.ATOM = pd.Series()
        self.atom_num = pd.Series()
        #here one space
        self.atom_name = pd.Series()
        self.alt_loc = pd.Series()
        self.res_name = pd.Series()
        self.chain_identif = pd.Series()
        self.res_seq_num = pd.Series()
        self.code_inset_res = pd.Series()
        #here 3 spaces
        self.X = pd.Series()
        self.Y = pd.Series()
        self.Z = pd.Series()
        self.occupancy = pd.Series()
        self.temp_fac= pd.Series()
        #here 6 spaces
        self.seg_id = pd.Series()
        self.ele_symb = pd.Series()
        self.charge = pd.Series()
        self.molname = pd.Series()

        ## assigning random molecule to each head in template
        for i in range(len(N_Z)): 
            if N_Z[i] < 0 : # it is the upper layer
                vyber = np.random.choice(layer_upper) # choosing randomly from list
                pos = layer_upper.index(vyber) 
                layer_upper.pop(pos) # removing chosen option from list

                mol_X, mol_Y, mol_Z = mol_dic[vyber].random_rotate("upper") #obtaining rotated coordinates of the phospholipid

            elif N_Z[i] > 0 : # its the lower layer
                vyber = np.random.choice(layer_lower) # choosing randomly from list
                pos = layer_lower.index(vyber)
                layer_lower.pop(pos) # removing chosen option from list

                mol_X, mol_Y, mol_Z = mol_dic[vyber].random_rotate("lower") #obtaining rotated coordinates of the phospholipid

            ## obtaining head coordinates, should be 0.0; 0.0; 0.0
            head_X = mol_X[mol_dic[vyber].head_idx] 
            head_Y = mol_Y[mol_dic[vyber].head_idx]
            head_Z = mol_Z[mol_dic[vyber].head_idx]

            ## translating the component head on the place of the template head,
            ## subtracting head coordinates just for sure,they should be 0.0; 0.0 ;0.0
            translate_X = N_X[i] - head_X 
            translate_Y = N_Y[i] - head_Y 
            translate_Z = N_Z[i] - head_Z  

            ## translating the rest of atoms accordingly
            newmol_X = [x - translate_X for x in mol_X] 
            newmol_Y = [x - translate_Y for x in mol_Y] 
            newmol_Z = [x - translate_Z for x in mol_Z] 

            ## collecting and creating the information of the bilayer
            self.ATOM =pd.concat([self.ATOM, mol_dic[vyber].ATOM], ignore_index=True)
            self.atom_num = pd.concat([self.atom_num, mol_dic[vyber].atom_num], ignore_index=True) # je jedno, co, jeste se bude prehazovat
            self.atom_name = pd.concat([self.atom_name, mol_dic[vyber].atom_name], ignore_index=True)
            self.alt_loc = pd.concat([self.alt_loc, mol_dic[vyber].alt_loc], ignore_index=True)
            self.res_name = pd.concat([self.res_name, mol_dic[vyber].res_name], ignore_index=True)
            self.chain_identif = pd.concat([self.chain_identif, mol_dic[vyber].chain_identif], ignore_index=True)
            self.res_seq_num = pd.concat([self.res_seq_num, pd.Series([i+1]*len(mol_dic[vyber].X))], ignore_index=True)# pdb cisluje od jedne
            self.code_inset_res = pd.concat([self.code_inset_res,mol_dic[vyber].code_insert_res],ignore_index=True)
            self.X = pd.concat([self.X, pd.Series(newmol_X)], ignore_index=True)
            self.Y = pd.concat([self.Y, pd.Series(newmol_Y)], ignore_index=True)
            self.Z = pd.concat([self.Z, pd.Series(newmol_Z)], ignore_index=True)
            self.occupancy = pd.concat([self.occupancy, mol_dic[vyber].occupancy], ignore_index=True)
            self.temp_fac= pd.concat([self.temp_fac, mol_dic[vyber].temp_fac], ignore_index=True)
            self.seg_id = pd.concat([self.seg_id, mol_dic[vyber].seg_id], ignore_index=True)
            self.ele_symb = pd.concat([self.ele_symb, mol_dic[vyber].ele_symb], ignore_index=True)
            self.charge = pd.concat([self.charge, mol_dic[vyber].charge], ignore_index=True)
            self.molname = pd.concat([self.molname, mol_dic[vyber].molname], ignore_index = True)

                                                   # alignment
        dataframe = {"_ATOM_": self.ATOM,          # left - in case of "HETATM"
                    "ATNUM": self.atom_num,        # right
                    "ATNA" : self.atom_name,       # left, but is complicated
                    "@" : self.alt_loc,            # -
                    "RESN" : self.res_name,        # right
                    "$": self.chain_identif,       # -
                    "RESQ" : self.res_seq_num,     # right
                    "€" : self.code_inset_res,     # -
                    "__dimX__" : self.X,           # right
                    "__dimY__" : self.Y,           # right
                    "__dimZ__" : self.Z,           # right
                    "OCCUPY" : self.occupancy,     # right
                    "TMPFAC" : self.temp_fac,      # right
                    "SGID" : self.seg_id,          # left
                    "EL" : self.ele_symb,          # right
                    "+-": self.charge,             # -
                    "molname": self.molname}       # not printed
        self.data = pd.DataFrame(dataframe)

        ## sorting that each bilyer component is only at one place
        self.data.sort_values(by=['molname','RESQ'], ignore_index=True, inplace=True)#self.data.sort_values(by=['RESN','RESQ'], ignore_index=True, inplace=True)
        self.data["ATNUM"] = pd.Series([i for i in range(1, len(self.data["ATNA"])+1)]) #correcting atom numbers
        self.data = self.data.fillna("") # getting rid of NaN

        ## correcting molecule numbers
        mol_lengths = []
        names = []
        for item in molecule:
            mol_lengths.append(len(item.X))
            #names.append(item.data["RESN"][0])
            names.append(item.molname[0])

        sorting = {"names": names,
                   "lens": mol_lengths,
                   "up": number_mol_upper,
                   "low": number_mol_lower} 
        sorting = pd.DataFrame(sorting)
        sorting.sort_values(by=["names"], ignore_index=True, inplace=True) # by default it is not known which molecule comes first, sorting by alphabet

        num_mol = []
        start = 1
        for ind in sorting.index: # is done in correct order 
            for j in range(start, start + sorting["up"][ind] + sorting["low"][ind]):
               num_mol =  num_mol + [j,]* sorting["lens"][ind]
            start = j + 1
        self.data["RESQ"] = pd.Series(num_mol)

    def print_PDB(self):
        """Prints the *.pdb file of the new bilayer. Containing new header with information
        about bilayer components in each leaflet, box and hydration vectors for further processing in GROMACS
        software. Then comes the *.pdb information. Footer containing "END" statement.
        """
        ## box and hydration parameters
        # box - adding 3 Ankstroms in each direction
        x_box_left = float(min(self.X) - 3)
        x_box_right = float(max(self.X) + 3)

        y_box_left = float(min(self.Y) - 3)
        y_box_right = float(max(self.Y) + 3)

        z_box_low = float(min(self.Z) - 3)
        z_box_up = float(max(self.Z) + 3)

        print("Actual system vectors: ")
        print("x:", '%.3f'%x_box_left, ";", '%.3f'%x_box_right)
        print("y:", '%.3f'%y_box_left, ";", '%.3f'%y_box_right)
        print("z:", ' %.3f'%z_box_low, ";", '%.3f'%z_box_up,'\n')

        # hydration
        #x_box = (abs(x_box_left)/10 + abs(x_box_right)/10) * 1.02
        #y_box = (abs(y_box_left)/10 + abs(y_box_right)/10) * 1.02
        #z_hydr =(420 / x_box / y_box ) * 2.5 #  TODO prepocist na volume
        #z_box = abs(z_box_low)/10 + abs(z_box_up)/10 + 2 * z_hydr

        # hydration
        x_box = (abs(x_box_left)/10 + abs(x_box_right)/10) * 1.02
        y_box = (abs(y_box_left)/10 + abs(y_box_right)/10) * 1.02
        z_hydr = 420 / x_box / y_box #  volume of 420 nm2 for hydration of 3356 water molecules in one leaflet should be sufficient
        z_box = abs(z_box_low)/10*1.01 + abs(z_box_up)/10*1.01 + 2 * z_hydr

        print("Input for cell creation and system hydration:")
        print("Box vector:", '%.3f'%(x_box),";", '%.3f'%(y_box), ";", '%.3f'%(z_box), "[nm]")
        print("Hydration vector: ",'%.3f'%(x_box), ";", '%.3f'%(y_box), ";", '%.3f'%(z_hydr), "[nm]")

        ## printing bilayer information in needed format (fwf)
        file = ""
        for i in range(len(self.data["__dimX__"])):
            line = '{0:<6}'.format(str(self.data["_ATOM_"][i])) + '{0:>5}'.format(str(self.data["ATNUM"][i])) + ' ' + '{0:<4}'.format(str(self.data["ATNA"][i])) + '{0:<1}'.format(str(self.data["@"][i]))\
                  + '{0:>4}'.format(str(self.data["RESN"][i])) + '{0:<1}'.format(str(self.data["$"][i])) + '{0:>4}'.format(str(self.data["RESQ"][i])) + '{0:<1}'.format(str(self.data["€"][i])) + '   ' + '{0:>8}'.format(str(round(self.data["__dimX__"][i],3)))\
                  + '{0:>8}'.format(str(round(self.data["__dimY__"][i],3))) + '{0:>8}'.format(str(round(self.data["__dimZ__"][i],3))) + '{0:>6}'.format(str(self.data["OCCUPY"][i])) +'{0:>6}'.format(str(self.data["TMPFAC"][i]))\
                  + '      ' + '{0:>4}'.format(str(self.data["SGID"][i])) + '{0:>2}'.format(str(self.data["EL"][i])) + '{0:>2}'.format(str(self.data["+-"][i] + "\n"))
            file = file + line

        # header and footer 
        string = ""
        upper = ""
        lower = ""
        for idx, item in enumerate(self.components):
            string = string + str(self.num_mol_low[idx] + self.num_mol_up[idx]) + " " + str(item.molecule_name) + " molecules, " 
            upper = upper + str(self.num_mol_up[idx]) + " " +  str(item.molecule_name) + ", "
            lower = lower + str(self.num_mol_low[idx]) + " " +  str(item.molecule_name) + ", "
        
        string = string[:-2] + " in 128 molecules phospholipid bilayer." 
        upper = upper[:-2] +  " molecules randomly distributed in the upper bilayer." 
        lower = lower[:-2] +  " molecules randomly distributed in the lower bilayer."

        now = datetime.now()
        date = datetime.today()

        self.head = f"REMARK    GENERATED BY BIC ON {now}\
        \nREMARK    INPUT GENERATION:\
        \nREMARK    {string}\
        \nREMARK    {upper}\
        \nREMARK    {lower}\
        \nREMARK    DATE:   {date}\
        \nREMARK    BOX VECTOR XYZ: {'%.3f'%(x_box)}; {'%.3f'%(y_box)}; {'%.3f'%(z_box)} [nm]\
        \nREMARK    HYDRATION VECTOR XYZ: {'%.3f'%(x_box)}, {'%.3f'%(y_box)}; {'%.3f'%(z_hydr)} [nm]\n"
        self.end = "END"

        ## print into *.pdb file
        with open(f"{self.filename}","w") as new_pdb:
            new_pdb.writelines(self.head)
            new_pdb.writelines(file)
            new_pdb.writelines(self.end)
