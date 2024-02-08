# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 13:40:39 2024

@author: jacob
"""

import unittest
from monty.serialization import loadfn, dumpfn
from monty.json import MSONable
import copy
import pickle
import time
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash
import sqlite3
from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from mc_analysis import ReportGenerator
import glob
import pickle
from monty.serialization import loadfn, dumpfn
import os
import copy
import networkx as nx
import time
import csv
import sqlite3
import copy
import operator

class Reaction(MSONable):
    def __init__(self, mpculecule_id_name, reaction_dict, kinetiscope_name, rate_coefficient, reaction_order):
        self.mpculecule_id_name = mpculecule_id_name
        self.reaction_dict = reaction_dict
        self.kinetiscope_name = kinetiscope_name
        self.rate_coefficient = rate_coefficient
        self.reaction_order = reaction_order

    def as_dict(self):
        # Convert reaction attributes to dictionary
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "mpculecule_id_name": self.mpculecule_id_name,
            "reaction_dict": self.reaction_dict,
            "kinetiscope_name": self.kinetiscope_name,
            "rate_coefficient": self.rate_coefficient,
            "reaction_order": self.reaction_order
        }

    @classmethod
    def from_dict(cls, d):
        # Create a Reaction object from a dictionary
        return cls(**d)
    
class ChemicalReaction(Reaction):
    def __init__(self, molecule_id_name, reaction_dict, kinetiscope_name, rate_coefficient, reaction_order, transition_state, classification):
        super().__init__(molecule_id_name, reaction_dict, kinetiscope_name, rate_coefficient, reaction_order)
        self.transition_state = transition_state
        self.classification = classification

    def as_dict(self):
        # Convert chemical reaction attributes to dictionary
        chemical_dict = super().as_dict()
        chemical_dict.update({
            "transition_state": self.transition_state,
            "classification": self.classification
        })
        return chemical_dict

class LightReaction(Reaction):
    pass  # No additional attributes for light reactions

class ElectronReaction(Reaction):
    pass  # No additional attributes for electron reactions

#charge neutralization: reactants of opposite charges, products neutral
#radical reactions require a total spin of 2
#charge transfer: charge does not change
 

def find_molecules_from_mpculeids(participants, mpculeid_to_mol_entry):
    """
    Converts a dictionary associating reaction "side" descriptors ("reactants" and/or "products") to
    reactant and product mpculeids into a dictionary associating those descriptors with
    lists of pymatgen Molecule objects.

    participants: a dictionary with keys limited to "reactants" and/or "products," 
                  where each value is a list of mpculeids representing reactants or products.
    mpculeid_to_mol_entry: a dictionary associating each mpculeid with its corresponding pymatgen Molecule object.

    Returns: new_dict, a new dictionary with the same keys ("reactants" and/or "products") 
             but updated values, where each value is a list of pymatgen Molecule objects associated with the mpculeids.
             Raises a KeyError if an mpculeid is not found in the dictionary.
    """
    new_dict = {}
    for key, side in participants.items():
        molecule_list = [mpculeid_to_mol_entry[mpculeid] for mpculeid in side]
        new_dict[key] = molecule_list
    return new_dict

# Function has been tested and passed the provided test cases
# Test cases covered scenarios such as basic conversion, missing mpculeids, and empty input
    
def name_molecule(mol, func_group_dict):
    """
    Takes a pymatgen Molecule object and generates a "name" for that object based
    on the functional groups present. Any atoms that are not associated with a 
    functional group will be be added to the name; this also means that species
    with no functional groups present will simply return their composition.
    
    Parameters
    ----------
    mol : pymatgen Molecule object
        The species we wish to name
    func_group_dict : dictionary
        A dictionary whos keys are the names of functional groups and values are
        the Molecule objects associated with those groups

    Returns
    -------
    name : string
        The generated "name" of the species

    """
    mol_graph = nx.Graph(mol.graph)
    mol_copy = copy.deepcopy(mol_graph)
    name = ""
    for n, group in func_group_dict.items(): 
        group_undirected_graph = nx.Graph(MoleculeGraph.with_local_env_strategy(group, OpenBabelNN(order = False)).graph) #Creates a Networkx (undirected) Graph object from func Molecule object
        while functional_group_present(mol_copy, group_undirected_graph)[0]: #want to repeat until func group is no longer present
            if n not in name:
                name = name + n + '_' #if present, add func group name to name
            else: #if the functional group is already in the name, we just increase the number present
                if name[-2].isnumeric():
                    new_number = int(name[-2]) + 1
                    name = name[0:len(name)-2]+ str(new_number) + '_'
                else:
                    name = name[0:len(name)-1] + '2_'
            mappings = list(functional_group_present(mol_copy, group_undirected_graph)[1]) #points to location of functional group in molecule
            for atom_index in mappings[0].keys():
                mol_copy.remove_node(atom_index)  #and remove the group
    if name and mol_copy: #i.e. if we have added some functional groups to the name but there remain atoms not in functional group
        comp = ""
        for node in mol_copy.nodes(data = True): #each node is a 2-tuple containing the index of that node and a dictionary with data related to that node
            species = node[1]['specie']
            if species in comp:
                species_index = comp.index(species)
                new_number = int(comp[species_index+1]) + 1
                comp = comp[0:len(comp)-1]+str(new_number)
            else:
                comp = comp + species + "1"
        name = name + comp + '_'
    elif name and not mol_copy: #i.e. our functional groups encompass all atoms in the molecule
        pass
    else: #i.e. the molecule contains no functional groups
        name = str(mol.molecule.composition).replace(" ", "")
        name = name + '_'
    if mol.molecule.charge == 1: #add charges to the end of the names
        name = name + "+" + str(mol.molecule.charge)
    else:
        name = name + str(mol.molecule.charge)
    return name

def functional_group_present(mol_graph, func):
    """
    Tests whether or not a given functional group is present in a molecule via
    testing if the Networkx graphical representation of a molecule contains a 
    subgraph that is isomorphic to the functional group.
    
    Parameters
    ----------
    mol_graph : Networkx Undirected Graph
        Graph of the molecule you want to name
    func : Networkx Undirected Graph
        Graph of the functional group

    Returns
    -------
    Boolean
        True if a subgraph of the molecule is isomorphic to the functional group
        graph, False otherwise.
    Generator over isomorphisms between a subgraph of G1 and G2.
        This is a generator of mappings between a subgraph of the molecule and
        the functional group--this helps us deterimine where the functional
        group is in the original molecule

    """
    nm = nx.isomorphism.categorical_node_match("specie", None) #ensures isomorphic graphs must have the same atoms
    isomorphism = nx.isomorphism.GraphMatcher(mol_graph, func, node_match = nm)
    return isomorphism.subgraph_is_isomorphic(), isomorphism.subgraph_isomorphisms_iter()

def stereoisomer_test(name_dict, mol_entry_id, mol_name):
    """
    Structural isomers--species with the same formula but different graphs--will return
    the same name when the name_molecule function is called. This function
    numbers each stereoisomer such that the new names are unique.

    Parameters
    ----------
    name_dict : dictionary
        dictionary whose keys are names and values are entry_ids of the species
        associated with that name
    mol_entry_id : string
        the unique entry id associated with a species
    mol_name : string
        the name generated for a species from the name_molecule function, which
        may not yet be unique

    Returns
    -------
    stereos : dictionary
        a dictionary of length 1 or 2, whose keys are now unique names of stereoisomers
        and values are entryids

    """
    stereos = {}
    if mol_name in name_dict and mol_entry_id != name_dict[mol_name]: #i.e. if name already in dictionary but graphs of molecules are different
        old_id = name_dict[mol_name]
        name_1 = mol_name + '_#1'
        name_2 = mol_name + '_#2'
        stereos.update({name_1: old_id, name_2: mol_entry_id})
    elif mol_name + '_#1' in name_dict and mol_entry_id != name_dict[mol_name + '_#1']: #i.e. if multiple stereoisomers are already in the dictionary
        current_max = 3
        while mol_name + '_#' + str(current_max) in name_dict and mol_entry_id != name_dict[mol_name + '_#' + str(current_max)]:
            current_max += 1
        new_name = mol_name + '_#' + str(current_max)
        stereos[new_name] = mol_entry_id
    return stereos

# def write_reaction(reaction_dict, mpculid_dict): #convert to strings of the appropriate format for kinetiscope
#     """
#     Kinetiscope can read in reactions of the format name1 + name2 => name3
#     + name4. This function uses the names we've generated to write strings for
#     reactions of this format.

#     Parameters
#     ----------
#     reaction_dict : dictionary
#         dictionary whose keys are "reactants" and "products" and whose values
#         are a list of mpculeids associated with a given reaction
#     mol_entry_id : dictionary
#         dictionary associating each mpculeid with its name for fast lookup

#     Returns
#     -------
#     rxn_eqn: string
#         the desired reaction written in the form: name1 + name2 => name3
#         + name4
#     """
#     rxn_eqn = ""
#     if len(reaction_dict["reactants"]) == 1:
#         name = mpculid_dict[reaction_dict["reactants"][0]]
#         rxn_eqn = name + " "
#     else:
#         first_reactant = reaction_dict["reactants"][0]
#         first_name = mpculid_dict[first_reactant]
#         rxn_eqn = first_name + " + "
#         second_reactant = reaction_dict["reactants"][1]
#         second_name = mpculid_dict[second_reactant]
#         if first_name == second_name:
#             # print(first_name)
#             rxn_eqn = "2 " + first_name + " " #if reaction is A + A -> B + C we have to write 2A rather than A + A
#             # print(rxn_eqn)
#         else:
#             rxn_eqn = rxn_eqn + second_name + " "
#     rxn_eqn = rxn_eqn + "=> "
#     if len(reaction_dict["products"]) == 1:
#         name = mpculid_dict[reaction_dict["products"][0]]
#         rxn_eqn = rxn_eqn + name
#     else:
#          first_product = reaction_dict["products"][0]
#          first_name = mpculid_dict[first_product]
#          # rxn_eqn = rxn_eqn + first_name + " + "
#          second_product = reaction_dict["products"][1]
#          second_name = mpculid_dict[second_product]
#          if first_name == second_name:
#              # print(first_name)
#              rxn_eqn += "2 " + first_name #if reaction is A + A -> B + C we have to write 2A rather than A + A
#              # print(rxn_eqn)
#          else:
#              rxn_eqn = rxn_eqn + first_name + " + " + second_name
#     return rxn_eqn

#ionization reactions have electron as a product, attachment have electrons as a reactant
#think I can find starting species with their indicies
#write light and electron reactions for starting species
#write LEE reactions for starting species
#output everything to a kinetiscope-readable file