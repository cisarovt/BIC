#!/usr/bin/ipython
# -*- coding: utf-8 -*-
# author: Císařová Terezie
"""Module with abstract classes for class creation in module used for phospholipid bilayer construction.
Imports 'abc' module.
"""
import abc

class PDBfile(abc.ABC):
    """Abstract class for creation of class working with information gained from *.pdb files
    """
    @abc.abstractmethod
    def __init__(self, pdb_filename:str, molname:str):
        self.molname = molname
        self.pdbfile = pdb_filename
