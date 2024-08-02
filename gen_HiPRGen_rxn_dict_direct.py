# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:54:29 2024

@author: JRMilton
"""
import os
from monty.serialization import loadfn, dumpfn
from Rxn_classes import HiPRGen_Reaction
from classify_ionization_reactions import (
    determine_broad_ionization_tag,
    narrow_down_ionization_type, 
    reaction_is_ionization
)
from classify_chemical_reactions import (
    determine_chemical_reaction_tag
)
from reaction_classification_utilities import (
    add_reaction_if_new,
    find_mpculeid_charge
    )

__version__ = '1.1.0'

def process_P2_reaction(rxn_dict, rxns_for_simulation, rxns_already_added):
    """
    Takes a reaction dictionary from phase 2, creates a HiPRGen reaction object
    from it, checks to see if it's new, and if it is, adds it to the list of 
    reactions for simulation and adds its name to the set of added names

    Parameters
    ----------
    rxn_dict : dict
        dictionary of the form: {"reactants":[mpculeids], "proudcts":[mpculeids]}
    rxns_for_simulation : list
        list of HiPRGen reaction objects we'll be simulating
    rxns_already_added : set
        names for reactions already added

    Returns
    -------
    tuple
        tuple containing (potentially) modified rxns_for_simulation and 
        rxns_already added with the new reaction added

    """
    new_rxn = HiPRGen_Reaction(rxn_dict,phase=2)
    tag = determine_chemical_reaction_tag(new_rxn)
    return add_reaction_if_new(new_rxn, tag, rxns_for_simulation, rxns_already_added)

def sort_pathways_list(pathways_list):
    """
    Takes the list of pathways, which are already sorted by weight, and then
    sorts by frequency, meaning if two paths have the same weight their order
    depends on their frequencies, with higher frequencies first

    Parameters
    ----------
    pathways_list : list
        unsorted list of pathways

    Returns
    -------
    sorted_list : list
        the list now sorted by weight and then frequency

    """
    sorted_list = sorted(pathways_list, key=lambda x: (x['weight'], x['frequency']))
    return sorted_list

def add_high_frequency_P2_reactions(directory, rxns_for_simulation, rxns_already_added, frequency_threshold=500):
    
    try:
        os.chdir(directory)  

        data = loadfn("reaction_tally.json")
        index_frequency_dict = data.get("pathways", {})
        index_rxn_dict = data.get("reactions", {})

        for rxn_index, rxn_dict in index_rxn_dict.items():
            rxn_frequency = index_frequency_dict.get(rxn_index, 0)
            
            if rxn_frequency >= frequency_threshold:
                rxns_for_simulation, rxns_already_added = process_P2_reaction(
                    rxn_dict, rxns_for_simulation, rxns_already_added
                )

    except FileNotFoundError:
        print(f"Error: The directory {directory} or the file 'reaction_tally.json' does not exist.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

def add_reaction_to_nested_dict(rxn_dict, keys, reaction):
    """Helper function to add a reaction to a nested dictionary."""
    d = rxn_dict
    for key in keys[:-1]:
        if key not in d:
            d[key] = {}
        d = d[key]
    if keys[-1] not in d:
        d[keys[-1]] = []
    d[keys[-1]].append(reaction)

def add_reaction_to_dict(rxn, rxn_dict, category, reaction_type=None, subclass=None, subtype=None, subcategory=None):
    keys = [category]
    possible_keys = [reaction_type, subclass, subtype, subcategory]
    for possible_key in possible_keys:
        if possible_key is None:
            break
        keys.append(possible_key)
    keys.append(reaction.tag)
    add_reaction_to_nested_dict(rxn_dict, keys, reaction)
    return rxn_dict

def add_unimolecular_reaction(rxn, rxn_dict):
    subclass = "fragmentation" if "fragmentation" in reaction.tag else "isomerization"
    return add_reaction_to_dict(rxn, rxn_dict, "chemical", "unimolecular", subclass)

def add_bimolecular_reaction(rxn, rxn_dict):
    subtype = determine_bimolecular_subtype(rxn)
    if "combination" in reaction.tag:
        subclass = "combination"
        return add_reaction_to_dict(rxn, rxn_dict, "chemical", "bimolecular", subclass, subtype)
    
    subclass = "biproduct"
    subcategories = [
        "proton_transfer",
        "hydride_abstraction",
        "H_atom_abstraction",
        "proton_coupled_electron_transfer",
        "electron_transfer",
        "reaction"
    ]
    
    for subcategory in subcategories:
        if subcategory in reaction.tag:
            return add_reaction_to_dict(rxn, rxn_dict, "chemical", "bimolecular", subclass, subtype, subcategory)
    
    return rxn_dict

def add_chemical_reactions(rxn, rxn_dict):
    if "fragmentation" in reaction.tag or "isomerization" in reaction.tag:
        return add_unimolecular_reaction(rxn, rxn_dict)
    return add_bimolecular_reaction(rxn, rxn_dict)

def add_ionization_reactions(rxn, rxn_dict):
    return add_reaction_to_dict(rxn, rxn_dict, "ionization")

def add_to_rxn_dict(rxn, rxn_dict, ionization_classifications):
    if reaction.tag in ionization_classifications:
        return add_ionization_reactions(rxn, rxn_dict)
    return add_chemical_reactions(rxn, rxn_dict)

def find_reactant_charges(rxn):
    reactant_1_charge = find_mpculeid_charge(rxn.reactants[0])
    reactant_2_charge = find_mpculeid_charge(rxn.reactants[1])
    return reactant_1_charge, reactant_2_charge

def reaction_is_neutral(rxn):
    reactant_1_charge, reactant_2_charge = find_reactant_charges(rxn)
    return reactant_1_charge == reactant_2_charge == 0

def reaction_is_ion_molecule(rxn):
    reactant_1_charge, reactant_2_charge = find_reactant_charges(rxn)
    charges_are_different = reactant_1_charge != reactant_2_charge
    one_charge_zero = reactant_1_charge == 0 or reactant_2_charge == 0
    return charges_are_different and one_charge_zero

def determine_bimolecular_subtype(rxn):
    if reaction_is_neutral(rxn):
        return "neutral"
    elif reaction_is_ion_molecule(rxn):
        return "ion-molecule"
    else:
        return "ion-ion"
    
# add all P1 reactions
P1_directory = "G:/My Drive/CRNs/071924_test_p1"
# P1_directory = "G:/My Drive/CRNs/071924_p1"
os.chdir(P1_directory)
P1_pathways_and_reactions = loadfn("reaction_tally.json") 
P1_rxn_dicts = P1_pathways_and_reactions["reactions"].values()
rxns_already_added = set()
rxns_for_simulation = []

for rxn_dict in P1_rxn_dicts:
    
    new_rxn = HiPRGen_Reaction(rxn_dict,phase=1)
    
    if reaction_is_ionization(new_rxn):
        
        tag = determine_broad_ionization_tag(new_rxn) 
    
    else:
        
        tag = determine_chemical_reaction_tag(new_rxn)
        
    rxns_for_simulation, rxns_already_added = \
        add_reaction_if_new(new_rxn, tag, rxns_for_simulation, rxns_already_added)
        
for reaction in rxns_for_simulation: #can narrow down here when we have all
                                      #reactions dealt with
    
    if reaction.tag == "attachment_or_recombination":
        
        reaction.tag = narrow_down_ionization_type(reaction, rxns_for_simulation)
 
#add all P2 reactions that fired >= 500 times

P2_directory = "G:/My Drive/CRNs/071924_test_p2"
# P2_directory = "G:/My Drive/CRNs/071924_p2"

add_high_frequency_P2_reactions(P2_directory, rxns_for_simulation, rxns_already_added, frequency_threshold=500)

network_products = loadfn("sink_report.json")

for product_dict in network_products.values():
    
    species_index = product_dict["species_index"]
    reactions_and_pathways = \
        loadfn(str(species_index) + "_pathway.json")
    all_pathways_list = list(reactions_and_pathways["pathways"])
    all_reactions = reactions_and_pathways["reactions"]
    sorted_pathways_list = sort_pathways_list(all_pathways_list)
    to_save_pathways = []
    
    for pathway_dict in sorted_pathways_list:
        pathway = pathway_dict["pathway"]
        to_save_pathways.append(pathway)
        number_pathways_saved = len(to_save_pathways)
        
        if number_pathways_saved >= 10:
            break
    
    for pathway in to_save_pathways:
        for reaction in pathway:
            rxn_dict = all_reactions.get(str(reaction), None)
            rxns_for_simulation, rxns_already_added = \
                        process_P2_reaction(rxn_dict, rxns_for_simulation, rxns_already_added)
                        
tagged_rxn_dict = {"ionization":{}, "chemical":{"unimolecular":{"fragmentation":{},"isomerization":{}}, "bimolecular":{"combination":{}, "biproduct":{}}}}
ionization_classifications = set(["positive_ionization", "electron_attachment",
                                  "electron_cation_recombination"])

for reaction in rxns_for_simulation:
    
    tagged_rxn_dict = \
        add_to_rxn_dict(
            reaction, tagged_rxn_dict, ionization_classifications
            )
    
dumpfn(tagged_rxn_dict,"HiPRGen_rxns_to_name.json")

# def print_reaction_dict_summary(dictionary):
#     """
#     Prints a summary of reactions in the nested dictionary.
    
#     Parameters:
#     - dictionary: A nested dictionary with a structure including categories,
#                   reaction types, subclasses, subtypes, and tags, each containing lists of reactions.
#     """
#     total_number_rxns = 0
    
#     def count_reactions(d):
#         nonlocal total_number_rxns
#         if isinstance(d, list):
#             total_number_rxns += len(d)
#         elif isinstance(d, dict):
#             for key, value in d.items():
#                 if isinstance(value, (list, dict)):
#                     count_reactions(value)
    
#     for category, types_dict in dictionary.items():
#         print(f"Category '{category}':")
#         count_reactions(types_dict)
    
#     print(f"Total number of reactions in the dictionary: {total_number_rxns}")
def print_reaction_dict_summary(dictionary):
    """
    Prints a detailed summary of reactions in the nested dictionary,
    including the values of each list and all associated dictionary keys.
    
    Parameters:
    - dictionary: A nested dictionary with categories, reaction types, subclasses,
                  subtypes, and tags, each containing lists of reactions.
    """
    def print_dict(d, level=0):
        """
        Recursively prints the dictionary structure, lists, and their values.
        
        Parameters:
        - d: The dictionary or list to print.
        - level: The current level of recursion, used for formatting.
        """
        indent = "  " * level
        if isinstance(d, list):
            # Print the length of the list and its items
            print(f"{indent}List with {len(d)} reactions")
            # for item in d:
            #     print(f"{indent}  {item}")  # Print each reaction in the list
        elif isinstance(d, dict):
            for key, value in d.items():
                print(f"{indent}Key '{key}':")
                print_dict(value, level + 1)  # Recursively print nested dictionary or list
    
    print("Reaction Dictionary Summary:")
    print_dict(dictionary)
    print(f"Total number of reactions in the dictionary: {count_reactions(dictionary)}")

def count_reactions(d):
    """ 
    Counts the total number of reactions in the nested dictionary.

    Parameters:
    - d: The dictionary or list to count reactions from.
    
    Returns:
    - The total number of reactions in the nested structure.
    """
    if isinstance(d, list):
        return len(d)
    elif isinstance(d, dict):
        return sum(count_reactions(value) for value in d.values())
    return 0

# def print_dict_lengths(dictionary):
#     total_number_rxns = 0
#     for superclass, tags_dict in dictionary.items():
#         print(f"Superclass '{superclass}':")
#         for tag, reactions in tags_dict.items():
#             value_length = len(reactions)
#             total_number_rxns += value_length
#             print(f"  Tag '{tag}' has {value_length} reactions")
#     print(f"total number of reactions in test set: {total_number_rxns}")

print_reaction_dict_summary(tagged_rxn_dict)