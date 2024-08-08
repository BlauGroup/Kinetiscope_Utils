# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 13:40:39 2024

@author: jacob
"""

__version__ = '1.0.0'

from monty.json import MSONable

class HiPRGen_Reaction(MSONable):
    def __init__(self, reaction_dict, classification_list=None, phase=None, tag=None):
        self.reactants = sorted([r for r in reaction_dict["reactants"] if r is not None])
        self.products = sorted([p for p in reaction_dict["products"] if p is not None])
        self.phase = phase
        self.classification_list = classification_list
        self.tag = tag
        self.name = self.generate_reaction_string()

    def generate_reaction_string(self):
        reactants_str = " + ".join(self.reactants)
        products_str = " + ".join(self.products)
        return f"{reactants_str} => {products_str}"
            
    def as_dict(self):
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "reactants": self.reactants,
            "products": self.products,
            "phase": self.phase,
            "classification_list": self.classification_list,
            "tag": self.tag,
            "name": self.name
        }

    def __str__(self):
        return self.name
    
    @classmethod
    def from_dict(cls, d):
        return cls(
            reaction_dict={"reactants": d["reactants"], "products": d["products"]},
            classification_list=d.get("classification_list"),
            phase=d.get("phase"),
            tag=d.get("tag")
        )
    
class Kinetiscope_Reaction(MSONable):
    def __init__(self, HiPRGen_rxn, kinetiscope_name, rate_coefficient, order, marker_species):
        self.HiPRGen_rxn = HiPRGen_rxn
        self.kinetiscope_name = kinetiscope_name
        self.rate_coefficient = rate_coefficient
        self.reaction_order = order
        self.marker_species = marker_species
        
    def as_dict(self):
        # Convert reaction attributes to dictionary
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "HiPRGen_rxn": self.HiPRGen_rxn,
            "kinetiscope_name": self.kinetiscope_name,
            "rate_coefficient": self.rate_coefficient,
            "reaction_order": self.reaction_order,
            "marker_species": self.marker_species
        }

    @classmethod
    def from_dict(cls, d):
        # Create a Reaction object from a dictionary
        return cls(
            HiPRGen_rxn=d.get("HiPRGen_rxn"),
            kinetiscope_name=d.get("kinetiscope_name"),
            rate_coefficient=d.get("rate_coefficient"),
            order=d.get("reaction_order"),
            marker_species=d.get("marker_species")
        )
    
    def __str__(self):
        return (f"Kinetiscope_Reaction(HiPRGen_rxn={self.HiPRGen_rxn}, "
                f"kinetiscope_name={self.kinetiscope_name}, "
                f"rate_coefficient={self.rate_coefficient}, "
                f"reaction_order={self.reaction_order}, "
                f"marker_species={self.marker_species})")