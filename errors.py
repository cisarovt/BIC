#!/usr/bin/ipython
# -*- coding: utf-8 -*-
# author: Císařová Terezie
"""This module contains error classes (parent class = Exception) 
for construction of phospholipid bilayer. Imports the 'glob' module.

Classes
-------
WrongMoleculeTypeNumber()
FileNotFoundError(file_input, *args: object)
WrongMoleculeNumber(value, message="Wrong number of molecules, should be 64 in total."))
"""
import glob

class WrongMoleculeTypeNumber(Exception):
    """Exception class for wrong number of molecule types input.
    Attributes
        ----------
        message : str
            Error message.
    """
    def __init__(self, *args: object):
        super().__init__(*args)
        self.message = "Wrong number of molecule types. Should be integer between 1 and 3."
    def __str__(self) -> str:
        return self.message

class FileNotFoundError(Exception):
    """Exception class for bad file path input.

    Attributes:
    ----------
    file_input: str
        File path refering to nonexisting file.
    """
    def __init__(self, file_input, *args: object):
        """Constructs all the necessary attributes for the FileNotFoundError exception object.

        Parameters
        ----------
             file_input: str
                File path refering to nonexisting file.

        """
        super().__init__(*args)

        file_input = str(file_input)
        separator1 = "\\"
        separator2 = "/"

        ending = "." + str(file_input.split(".",1)[1])

        if separator1 in file_input:
            path = file_input.split(separator1, 1)[0] + separator1
        elif separator2 in file_input:
            path = file_input.split(separator2, 1)[0] + separator2
        else:
            path = ".\\"

        list_of_files = glob.glob(path + "*" + ending)

        self.message = f"Given file {file_input} does not exist.\n \
            Maybe you intended to use one of those files in the directory {path}:\n \
                  {list_of_files}"
    def __str__(self) -> str:
        return self.message

class WrongMoleculeNumber(Exception):
    """ Class for raising error when the sum of given molecules is not 128.

    Attributes
    ----------
    value : list
        List of numbers of molecules.
    message : str
        Error message.
    """
    def __init__(self, value, message="Wrong number of molecules, should be 64 in total."):
        self.mol_nums = value
        self.message = message
        self.num_mol = sum(i for i in self.mol_nums)
        super().__init__(self, message)
    def __str__(self):
        return f"Given coumts of molecules: {self.mol_nums},\
              sums to {self.num_mol} -- {self.message}"
    