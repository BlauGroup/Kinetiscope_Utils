# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 13:40:39 2024

@author: jacob
"""

__version__ = '1.0.2'

from monty.json import MSONable
import networkx as nx
from pymatgen.core.structure import Molecule


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


class RecombMolEntry(MSONable):
    def __init__(self, graph, molecule, entry_id):
        self.graph = graph
        self.molecule = molecule
        self.entry_id = entry_id

    def as_dict(self):
        """
        Serialize the object to a dictionary for MSONable compatibility.
        """
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "graph": nx.node_link_data(self.graph),  # Converts the NetworkX graph to serializable data
            "molecule": self.molecule.as_dict(),  # Assumes pymatgen Molecule object has as_dict
            "entry_id": self.entry_id
        }

    @classmethod
    def from_dict(cls, d):
        """
        Deserialize an object from a dictionary.
        """
        graph = nx.node_link_graph(d["graph"])  # Reconstruct the NetworkX graph
        molecule = Molecule.from_dict(d["molecule"])  # Assumes pymatgen Molecule has from_dict
        entry_id = d["entry_id"]
        return cls(graph, molecule, entry_id)

    def __str__(self):
        """
        Provides a readable string representation of the object.
        """
        return (
            f"RecombMolEntry(\n"
            f"  entry_id: {self.entry_id},\n"
            f"  molecule: {self.molecule.composition},\n"
            f"  graph nodes: {len(self.graph.nodes)},\n"
            f"  graph edges: {len(self.graph.edges)}\n"
            f")"
        )

class Kinetiscope_Reaction(MSONable):  # TODO rename this to better represent what it is
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
    
class Kinetiscope_simulation_reaction(MSONable):
    def __init__(self, kinetiscope_name, kinetiscope_index, selection_freq):
        self.kinetiscope_name = kinetiscope_name
        self.kinetiscope_index = kinetiscope_index
        self.selection_freq = selection_freq
        self.reactants = self.find_reactants()
        self.products = self.find_products()
    
    def find_reactants(self):
        return self.find_species(0)

    def find_products(self):
        return self.find_species(1)
    
    def find_species(self, index):
        unsorted_species = self.find_both_sides()[index]
        unsorted_species_list = unsorted_species.split(" + ")
        sorted_species = sorted(species.strip() for species in unsorted_species_list)
        return sorted_species
    
    def find_both_sides(self):
        return self.kinetiscope_name.split( "=> ")
    
    def as_dict(self):
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "kinetiscope_name": self.kinetiscope_name,
            "kinetiscope_index": self.kinetiscope_index,
            "reactants": self.reactants,
            "products": self.products,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            kinetiscope_name=d["kinetiscope_name"],
            kinetiscope_index=d["kinetiscope_index"]
        )
    
    def __str__(self):
        return (f"Kinetiscope Simulation Reaction:\n"
                f"  Kinetiscope Name: {self.kinetiscope_name}\n"
                f"  Kinetiscope Index: {self.kinetiscope_index}\n"
                f"  Selection frequency: {self.selection_freq}\n"
                f"  Reactants: {self.reactants}\n"
                f"  Products: {self.products}")
