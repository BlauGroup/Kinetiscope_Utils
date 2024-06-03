# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 13:40:39 2024

@author: jacob
"""
__version__ = '1.0.0'
# import unittest
from monty.serialization import loadfn
from monty.json import MSONable
# import copy
# import pickle
# import time
# from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash
# import sqlite3
# from pymatgen.core.structure import Molecule
# from pymatgen.analysis.graphs import MoleculeGraph
# from pymatgen.analysis.local_env import OpenBabelNN
# from mc_analysis import ReportGenerator
# import glob
# import pickle
# from monty.serialization import loadfn, dumpfn
# import os
# import copy
# import networkx as nx
# import time
# import csv
# import sqlite3
# import copy
# import operator

# class HiPRGen_Reaction(MSONable):
#     def __init__(self, reaction_dict, phase=None, tag=None):
#         self.reactants = [r for r in reaction_dict["reactants"] if r is not None]
#         self.products = [p for p in reaction_dict["products"] if p is not None]
#         self.phase = phase
#         self.reaction_hash = self.write_reaction_hash()
    
#     def write_reaction_hash(self): #mpculeids should be unique, so their sum also should be
#         sorted_reactant_mpculeids = sorted(self.reactants)
#         sorted_product_mpculeids = sorted(self.products)
#         reaction_hash = ''.join(sorted_reactant_mpculeids) + ''.join(sorted_product_mpculeids)
#         return reaction_hash
            
#     def as_dict(self):
#         # Convert reaction attributes to dictionary
#         return {
#             "@module": self.__class__.__module__,
#             "@class": self.__class__.__name__,
#             "reactants": self.reactants,
#             "products": self.products,
#             "phase": self.phase,
#             "reaction_hash": self.reaction_hash
#         }

#     def __str__(self):
#         # Customized string representation
#         return str({
#             "reactants": self.reactants,
#             "products": self.products,
#             "phase": self.phase,
#             "reaction_hash": self.reaction_hash
#         })
    
#     @classmethod
#     def from_dict(cls, d):
#         # Create a Reaction object from a dictionary
#         return cls(**d)
class HiPRGen_Reaction(MSONable):
    def __init__(self, reaction_dict, phase=None, tag=None):
        self.reactants = [r for r in reaction_dict["reactants"] if r is not None]
        self.products = [p for p in reaction_dict["products"] if p is not None]
        self.phase = phase
        self.reaction_hash = self.write_reaction_hash()
        self.tag = tag
    
    def write_reaction_hash(self): #after testing these are unique
        reaction_hash = ''.join(self.reactants) + "=>"  + ''.join(self.products)
        return reaction_hash
            
    def as_dict(self):
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "reactants": self.reactants,
            "products": self.products,
            "phase": self.phase,
            "reaction_hash": self.reaction_hash,
            "tag":self.tag
        }

    def __str__(self):
        return str({
            "reactants": self.reactants,
            "products": self.products,
            "phase": self.phase,
            "reaction_hash": self.reaction_hash,
            "tag":self.tag
        })
    
    @classmethod
    def from_dict(cls, d):
        return cls(d)
    
class Kinetiscope_Reaction(MSONable):
    def __init__(self, kinetiscope_name, rate_coefficient, order, tag):
        self.kinetiscope_name = kinetiscope_name
        self.rate_coefficient = rate_coefficient
        self.reaction_order = order
        self.tag = tag
        
    def as_dict(self):
        # Convert reaction attributes to dictionary
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "kinetiscope_name": self.kinetiscope_name,
            "rate_coefficient": self.rate_coefficient,
            "reaction_order": self.reaction_order,
            "tag": self.tag
        }
    
    @classmethod
    def from_dict(cls, d):
        # Create a Reaction object from a dictionary
        return cls(**d)
    
# class ChemicalReaction(Reaction):
#     def __init__(self, molecule_id_name, reaction_dict, kinetiscope_name, rate_coefficient, reaction_order, transition_state, classification):
#         super().__init__(molecule_id_name, reaction_dict, kinetiscope_name, rate_coefficient, reaction_order)
#         self.transition_state = transition_state
#         self.classification = classification

#     def as_dict(self):
#         # Convert chemical reaction attributes to dictionary
#         chemical_dict = super().as_dict()
#         chemical_dict.update({
#             "transition_state": self.transition_state,
#             "classification": self.classification
#         })
#         return chemical_dict

# class LightReaction(Reaction):
#     pass  # No additional attributes for light reactions

# class ElectronReaction(Reaction):
#     pass  # No additional attributes for electron reactions

#charge neutralization: reactants of opposite charges, products neutral
#radical reactions require a total spin of 2
#charge transfer: charge does not change
 

#ionization reactions have electron as a product, attachment have electrons as a reactant
#think I can find starting species with their indicies
#write light and electron reactions for starting species
#write LEE reactions for starting species
#output everything to a kinetiscope-readable file