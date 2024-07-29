# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 13:40:39 2024

@author: jacob
"""

__version__ = '1.0.0'

from monty.json import MSONable

class HiPRGen_Reaction(MSONable):
    def __init__(self, reaction_dict, phase=None, tag=None):
        self.reactants = sorted([r for r in reaction_dict["reactants"] if r is not None])
        self.products = sorted([p for p in reaction_dict["products"] if p is not None])
        self.phase = phase
        # self.reaction_hash = self.write_reaction_hash()
        self.tag = tag
        self.name = self.generate_reaction_string()

    def generate_reaction_string(self):
        reactants_str = " + ".join(self.reactants)
        products_str = " + ".join(self.products)
        return f"{reactants_str} => {products_str}"
    
    # def write_reaction_hash(self): #after testing these are unique
    #     reaction_hash = ''.join(self.reactants) + "=>"  + ''.join(self.products)
    #     return reaction_hash
            
    def as_dict(self):
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "reactants": self.reactants,
            "products": self.products,
            "phase": self.phase,
            # "reaction_hash": self.reaction_hash,
            "tag": self.tag,
            "name": self.name
        }

    def __str__(self):
        return self.name
    
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