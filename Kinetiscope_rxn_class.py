# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 13:40:39 2024

@author: jacob
"""

# import unittest
# from monty.serialization import loadfn, dumpfn
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

class HiPRGen_Reaction(MSONable):
    def __init__(self, mpcule_id_name):
        self.mpcule_id_name = mpcule_id_name
        self.reaction_dict = self.generate_reaction_dict(mpcule_id_name)
        self.reactants = self.reaction_dict["reactants"]
        self.products = self.reaction_dict["products"]
        
    @staticmethod
    def generate_reaction_dict(mpcule_id_name):
        """
        For a given reaction, we often need to access the reactants and products.
        This just uses the given mpculeid_name of the reaction to associate
        the reactants and products of a reaction with a Reaction object.

        Parameters
        ----------
        mpcule_id_name : str
            the name of a chemical reaction, written as mpculeid + mpculeid ->
            mpculeid + mpculeid. These will always involve 1-2 reactants and/or
            products.
            
        Returns
        -------
        dict
            a dictionary of the form: {"reactants": reactant_list, "products":
                                       product_list}

        """
        both_sides = mpcule_id_name.split('->')
        reactants = both_sides[0].split('+') if '+' in both_sides[0] else [both_sides[0]]
        products = both_sides[1].split('+') if '+' in both_sides[1] else [both_sides[1]]
        return {"reactants":reactants, "products":products}
            
    def as_dict(self):
        # Convert reaction attributes to dictionary
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "mpculecule_id_name": self.mpculecule_id_name,
            "reaction_dict": self.reaction_dict,
            "reactants": self.reactants,
            "products": self.products
        }

    @classmethod
    def from_dict(cls, d):
        # Create a Reaction object from a dictionary
        return cls(**d)
    
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