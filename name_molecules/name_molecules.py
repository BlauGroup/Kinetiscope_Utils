# -*- coding: utf-8 -*-
"""
Created on Wed May 10 13:26:19 2023

@author: JRMilton
"""

import copy
import networkx as nx
from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
import glob
from monty.serialization import dumpfn
import os
import pickle
import sys

"""

This module takes the mol_entries.pickle file and generates a unique name 
for each species in our reaction network such that each name is somewhat human 
interpretable and fits requirements for names in Kinetiscope. Its output is a 
json of the form name:mpculeid.

"""

def functional_group_present(mol_graph, func): 
    """
    Tests whether or not a given functional group is present in a molecule via
    testing if the Networkx graphical representation of a molecule contains a 
    subgraph that is isomorphic to the functional group.
    
    Parameters
    ----------
    mol_graph : Networkx undirected graph
        Graph of the molecule you want to name
    func : Networkx undirected graph
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

def update_species_name(species_name, func_group_name, func_name_already_added):
    """
    Adds the name of the functional group present within the molecule to its
    species name, accounting for the name of the functional group already 
    having been added to the name. 

    Parameters
    ----------
    species_name : string
        current name generated for the species
    func_group_name : string
        name of the functional group being added to the species name
    func_name_already_added : Boolean
        True if the name of this functional group is already in the species name,
        False otherwise.

    Returns
    -------
    species_name : string
        the updated name of this species

    """
    
    if func_name_already_added:
        
        if species_name[-2].isnumeric():
            
            current_number_present = int(species_name[-2])
            new_number = current_number_present + 1
            species_name = species_name[:-2] + str(new_number)
            
        else:
            
            species_name = species_name[:-1] + str(2)  
            
    else:
        
        species_name += func_group_name 
        
    species_name = species_name + "_"
    
    return species_name

def update_remaining_species_graph(remaining_species_graph, func_group_mappings):
    """
    Networkx returns a mapping between a subgraph of our species and a functional
    group graph. This function deletes the atoms from a copy of the species graph
    to "remove" that instance of the functional group.

    Parameters
    ----------
    remaining_species_graph : networkx Graph object
        The graph of a given species after the nodes associated with a
        functional group have been removed.
    func_group_locations : generator
        generator over isomorphisms between a subgraph of two graphs

    Returns
    -------
    remaining_species_graph : networkx Graph object
        the graph of the species with the functional group removed.

    """
    
    to_remove = func_group_mappings[0].keys() #we only care about one possible mapping b/c we remove iteratively
    
    for atom_index in to_remove:
        
        remaining_species_graph.remove_node(atom_index)
    
    return remaining_species_graph

def update_name_and_graph(species_name, func_group_name, remaining_species_graph, func_group_mappings):
    """
    Converts the subgraph isomorphism mapping to a list, tests the name of the
    functional group is already present in species_name, and calls functions
    to update the species name and graph before finally returning them.

    Parameters
    ----------
    species_name : string
        current name of the species
    func_group_name : string
        name of the functional group to be added to species_name
    remaining_species_graph : networkx Graph object
        current graph of the species
    func_group_mappings : generator
        generator over isomorphisms between a subgraph of one graph and
        another graph

    Returns
    -------
    species_name : string
        the species name with the functional group added
    remaining_species_graph : networkx Graph object
        graph of the species after the functional group is removed

    """
    
    mappings_list = list(func_group_mappings)
    func_name_already_added = func_group_name in species_name
    species_name = update_species_name(species_name, func_group_name, func_name_already_added)
    remaining_species_graph = update_remaining_species_graph(remaining_species_graph, mappings_list)
    
    return species_name, remaining_species_graph

def add_composition(species_name, remaining_species_graph):
    """
    After all functional group names have been added, this function adds to the
    name the composition of everything that remains.

    Parameters
    ----------
    species_name : str
        the current name
    remaining_species_graph : networkx undirected graph
        the graph of the species we're naming with any functional groups
        removed
    species_charge: int
        charge of the species
    Returns
    -------
    species_name : str
        the name updated with composition and charge

    """   
    
    atom_composition = {}
    
    for node in remaining_species_graph.nodes(data=True):
        
        element = node[1]['specie']
        
        if element in atom_composition:
            
            atom_composition[element] += 1
            
        else:
            
            atom_composition[element] = 1

    for element, count in atom_composition.items():
        species_name += f"{element}{count}"
        
    species_name = species_name + "_"
    
    return species_name

def add_charge(species_name, charge):
    """
    Simple function for adding the charge of a species to its name.

    Parameters
    ----------
    species_name : string
        current name of the species
    charge : int
        charge of the species

    Returns
    -------
    species_name : string
        species name with the charge added to the end

    """
    
    charge_str = str(charge)
    charge_suffix = "+" + charge_str if charge >= 1 else charge_str
    species_name += charge_suffix
    
    return species_name

def generate_species_name(species_graph, species_pymatgen_mol, func_group_dict):
    """
    Generates a name for a species based on the functional groups it contains.
    
    Parameters
    ----------
    species_dict : dict
        a dictionary containing information related to the species
    func_group_dict : dict
        a dictionary whose keys are names (strings) and whose values are 
        networkx undirected graphs associated with a given functional group
    
    Returns
    -------
    species_name : string
        the name generated for the species
    
    """
    
    remaining_species_graph = copy.deepcopy(species_graph)
    species_name = ""

    for func_group_name, func_group_graph in func_group_dict.items():
        
        func_group_is_present, func_group_mappings = \
            functional_group_present(remaining_species_graph, func_group_graph)
        
        while func_group_is_present:
            
            species_name, remaining_species_graph = \
                update_name_and_graph(species_name, func_group_name, remaining_species_graph, func_group_mappings)
            
            func_group_is_present, func_group_mappings = \
                functional_group_present(remaining_species_graph, func_group_graph)
    
    if remaining_species_graph: #i.e. not all atoms are accounted for in the name
    
        species_name = add_composition(species_name, remaining_species_graph)
        
    species_charge = species_pymatgen_mol.charge
    species_name = add_charge(species_name, species_charge)
    
    return species_name

def name_stereoisomers(stereoisomer_list, name):
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
    
    new_stereos_list = []
    
    if len(stereoisomer_list) == 1: #i.e. only one other species with this name
    
        name_1 = name + '_#1'
        name_2 = name + '_#2'
        new_stereos_list.append(name_1)
        new_stereos_list.append(name_2)
        
    else:
        
        current_max_num = len(stereoisomer_list)
        new_max_num_str = str(current_max_num + 1)
        new_isomer = name + "_#" + new_max_num_str
        new_stereos_list = new_stereos_list + stereoisomer_list
        new_stereos_list.append(new_isomer)
        
    return new_stereos_list

def update_names(test_name, stereo_dict, name_mpcule_dict, mpculeid):
    """
    When we generate a name for a species that is a stereoisomer of a species
    we have already named, this function updates the dictionary containing names
    to include the new stereoisomer.

    Parameters
    ----------
    test_name : string
        the name generated for this new stereoisomer
    stereo_dict : dict
        dictionary associating test_name with all of its current stereoisomers
    name_mpcule_dict : dict
        dictionary associating the name of a species with its mpculeid
    mpculeid : string
        an identifier assigned to each species

    Returns
    -------
    None, updates stereo_dict and name_mpcule_dict

    """
    
    current_stereos = stereo_dict.get(test_name)
    new_stereos = name_stereoisomers(current_stereos, test_name)
    stereo_dict[test_name] = new_stereos 
    
    old_mpculeid = name_mpcule_dict.get(test_name, False)
    
    if old_mpculeid: 
        
        name_mpcule_dict[new_stereos[-2]] =  old_mpculeid #we assign this one just so the following line is valid regardless of whether or not this test fired
        name_mpcule_dict.pop(test_name)
        
    name_mpcule_dict[new_stereos[-1]] = mpculeid

def reorder_dict_keys(input_dict, key_order):
    """
    Some of the functional groups contain each other--for example, a carboxylic
    acid contains a carbonyl group, so we need to add more specific functional
    groups to names before less specific ones. This function reorders the 
    functional group dictionary based on the list give with the desired order.
    If order doesn't matter for a given functional group, it's just added as
    a key in alphabetical order after the groups for which order does matter.

    Parameters
    ----------
    input_dict : dict
        the dictionary we're reordering
    key_order : list
        list of keys in input_dict, given in the order we want those keys to be
        in in the dict we return

    Returns
    -------
    reordered_dict : dict
        the dictionary with our new order

    """
    ordered_keys = [key for key in key_order if key in input_dict]
    remaining_keys = sorted([key for key in input_dict if key not in key_order]) 
    new_order = ordered_keys + remaining_keys
    reordered_dict = {key: input_dict[key] for key in new_order}
    
    return reordered_dict
    
#Change directory to the functional groups folder
func_groups_dir = r"func_groups"
os.chdir(func_groups_dir)

# Process functional groups' XYZ files
print('Associating functional groups with their Molecule objects...')
func_group_dict = {}

for filename in glob.glob('*.xyz'):
    func_group = Molecule.from_file(filename)  # Load functional group as a pymatgen Molecule
    func_group_mol_graph = MoleculeGraph.with_local_env_strategy(func_group, OpenBabelNN(order=False)) #build pymatgen MoleculeGraph
    func_group_undirected_graph = nx.Graph(func_group_mol_graph.graph) #convert graph to networkx undirected graph
    
    name = filename.replace('.xyz', '')
    func_group_dict[name] = func_group_undirected_graph

order_list = ["PtBMAb", "PHSb", "TPS", "phol", "ester", "COOH", "COO"]
func_group_dict = reorder_dict_keys(func_group_dict, order_list)

pickle_directory = \
    "C:/Users/JRMilton/Kinetiscope_Utils/name_molecules"
    
os.chdir(pickle_directory)

#Note: the pickle file needs to be in the same folder as HiPRGen for the 
#mol_entry objects to be loaded

print('Loading mol_entries.pickle...')

with open('mol_entries.pickle', 'rb') as f: #loads HiPRGen mol_entry objs
        
    mol_entries = pickle.load(f)
    
print('Done!')

number_to_name = len(mol_entries)

print("Generating names...")
name_mpcule_dict = {}
stereo_dict = {} #assocites a given base name with all of its stereoisomers

# Process each molecule and generate a species name
for mol_entry in mol_entries:
    
    graph = mol_entry.graph
    molecule = mol_entry.molecule
    
    species_name = generate_species_name(graph, molecule, func_group_dict)
    
    mpcule_id = mol_entry.entry_id
    
    if species_name in stereo_dict:
        
        update_names(species_name, stereo_dict, name_mpcule_dict, mpcule_id)
        
    else:
        
        name_mpcule_dict[species_name] = mpcule_id
        stereo_dict[species_name] = [species_name]
        
    if len(species_name) >= 33: #names in kinetiscope are limited to 32 chars
        
        print(f"Generated species name is too long! {species_name}")
        print(f"Length of name is: {len(species_name)}")
        print("consider changing the names of functional groups to be shorter")
        print("Aborting")
        sys.exit()
        

print("Done!")

number_named = len(name_mpcule_dict)
total_num_names = len(name_mpcule_dict)

if number_named != number_to_name:
    
    print("Failed to name the correct number of species")
    print(f"Expected number of names: {number_to_name}")
    print(f"Actual number of names generated: {number_named}" )
    print("Aborting")
    sys.exit()

print("saving names to files...")

dumpfn(name_mpcule_dict, "name_full_mpculeid_081924.json")

print("Done!")