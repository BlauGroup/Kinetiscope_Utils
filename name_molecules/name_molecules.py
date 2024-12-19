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
from monty.serialization import loadfn, dumpfn
import os
import pickle
import sys
sys.path.append('../common')
from utilities import correct_path_change_dir

#TODO incorporate spin into names

"""

This module takes the mol_entries.pickle file and generates a unique name 
for each species in our reaction network such that each name is somewhat human 
interpretable and fits requirements for names in Kinetiscope. Its output is a 
json of the form name:mpculeid.

"""


class NameLengthError(Exception):
    """
    Custom exception for names exceeding the length limit in Kinetiscope.
    """
    def __init__(self, name, length, message=None):
        self.name = name
        self.length = length
        if message is None:
            message = (
                f"Generated species name is too long! {name}\n"
                f"Length of name is: {length}\n"
                "Consider changing the names of functional groups to be shorter."
            )
        super().__init__(message)
        
        
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
    
    # ismorphism.subgraph_is_isomorphic returns True if a subgraph of G1 is 
    # isomorphic to G2.
    
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

def generate_species_name(species_graph, charge, func_group_dict):
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
        
    # species_charge = species_pymatgen_mol.charge
    species_name = add_charge(species_name, charge)
    
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
    
def NameSpecies(graph, charge, func_group_dict, mpculeid,  stereo_dict, name_mpcule_dict):

    test_name = generate_species_name(graph, charge, func_group_dict)

    if test_name in stereo_dict:
        
        update_names(test_name, stereo_dict, name_mpcule_dict, mpculeid)
        
    else:
        
        name_mpcule_dict[test_name] = mpculeid
        stereo_dict[test_name] = [test_name]
        
    if len(test_name) >= 33:  # Names in Kinetiscope are limited to 32 chars
        raise NameLengthError(test_name, len(test_name))

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

def build_recomb_func_group_dict(func_group_dict):
    """
    Builds a recombination functional group dictionary where the keys in
    func_group_dict are replaced with shorthand versions from shorthand_dict
    if available. If no shorthand exists, the original key is retained. The 
    names for the recombinants are too long for kinetiscope, while those we
    wrote previously (for other species in our network) are not, so this
    just generates special, shorter names for recombinants.

    Parameters
    ----------
    func_group_dict : dict
        The original dictionary of functional groups.

    Returns
    -------
    dict
        A new dictionary with updated keys using shorthand where applicable.
    """
    shorthand_dict = {
        "COO": "CX",
        "ester": "e",
        "PHSb": "HS",
        "phol": "po",
        "phyl": "py",
        "tbut": "t",
        "tBMb": "TBM"
    }
    
    return {shorthand_dict.get(key, key): value for key, value in func_group_dict.items()}


kinetiscope_files_dir = (
    r"G:\My Drive\Kinetiscope\production_simulations_092124"
)

correct_path_change_dir(kinetiscope_files_dir)

recomb_mol_entries = loadfn("recomb_mol_entries_121624.json")
# sys.exit()
# Change directory to the functional groups folder
func_groups_dir = r"C:\Users\jacob\Kinetiscope_Utils\name_molecules\func_groups"
correct_path_change_dir(func_groups_dir)

# Process functional groups' XYZ files
print('Associating functional groups with their Molecule objects...')
func_group_dict = {}

for filename in glob.glob('*.xyz'):
    func_group = Molecule.from_file(filename)  # Load functional group as a pymatgen Molecule
    func_group_mol_graph = MoleculeGraph.with_local_env_strategy(func_group, OpenBabelNN(order=False)) #build pymatgen MoleculeGraph
    func_group_undirected_graph = nx.Graph(func_group_mol_graph.graph) #convert graph to networkx undirected graph
    
    name = filename.replace('.xyz', '')
    func_group_dict[name] = func_group_undirected_graph
    
print("Done!")
print("Naming recombinants...")

order_list = ["tBMb", "PHSb", "TPS", "phol", "ester", "COOH", "COO"]
func_group_dict = reorder_dict_keys(func_group_dict, order_list)
recomb_func_group_dict = build_recomb_func_group_dict(func_group_dict)
name_mpcule_dict = {}
stereo_dict = {} #assocites a given base name with all of its stereoisomers
 
for recomb_mol_entry in recomb_mol_entries:
    graph = nx.Graph(recomb_mol_entry.graph)
    charge = recomb_mol_entry.charge
    mpculeid = recomb_mol_entry.entry_id
    NameSpecies(graph, charge, recomb_func_group_dict, mpculeid,  stereo_dict, name_mpcule_dict)
    
for name in name_mpcule_dict.keys():
    if len(name) > 32:
        raise NameLengthError(name, len(name))

print("Done!")
# sys.exit()
   
pickle_directory = \
    "C:/Users/jacob/Kinetiscope_Utils/name_molecules"
    
os.chdir(pickle_directory)

#Note: the pickle file needs to be in the same folder as HiPRGen for the 
#mol_entry objects to be loaded. Also need to make sure to use the pickle file
#for phase 2, because the file from phase 1 will cause an error, presumably
#due to the presence of the "electron species"

print('Loading mol_entries.pickle...')

with open('mol_entries.pickle', 'rb') as f: #loads HiPRGen mol_entry objs
        
    try:
        
        mol_entries = pickle.load(f)
    
    except AttributeError:
        
        print("Need to use the pickle file from phase 2, not phase 1")
        print("The 'electron species' in phase 1 causes issues here")
        sys.exit()
    
print('Done!')

# number_to_name = len(mol_entries)

print("Generating other species names...")

# Process each molecule and generate a species name
for mol_entry in mol_entries:
    
    graph = mol_entry.graph
    molecule = mol_entry.molecule
    charge = molecule.charge
    mpculeid = mol_entry.entry_id
    
    NameSpecies(graph, charge, func_group_dict, mpculeid,  stereo_dict, name_mpcule_dict)
    
    # species_name = generate_species_name(graph, charge, func_group_dict)
    
    # # mpcule_id = mol_entry.entry_id
    
    # if species_name in stereo_dict:
        
    #     update_names(species_name, stereo_dict, name_mpcule_dict, mpcule_id)
        
    # else:
        
    #     name_mpcule_dict[species_name] = mpcule_id
    #     stereo_dict[species_name] = [species_name]
        
    # if len(species_name) >= 33: #names in kinetiscope are limited to 32 chars
        
    #     print(f"Generated species name is too long! {species_name}")
    #     print(f"Length of name is: {len(species_name)}")
    #     print("consider changing the names of functional groups to be shorter")
    #     print("Aborting")
    #     sys.exit()
for name in name_mpcule_dict.keys():
    if len(name) > 32:
        raise NameLengthError(name, len(name))

print("Done!")

number_named = len(name_mpcule_dict)
number_to_name = len(mol_entries) + len(recomb_mol_entries)

if number_named != number_to_name:
    
    print("Failed to name the correct number of species")
    print(f"Expected number of names: {number_to_name}")
    print(f"Actual number of names generated: {number_named}" )
    print("Aborting")
    sys.exit()

print("saving names to files...")

dumpfn(name_mpcule_dict, "name_mpculeid_withrecombs_121824.json")

print("Done!")