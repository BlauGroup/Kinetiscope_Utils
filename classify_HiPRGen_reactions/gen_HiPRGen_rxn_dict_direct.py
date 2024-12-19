# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:54:29 2024

@author: JRMilton
"""

import os
from monty.serialization import loadfn, dumpfn
from classify_ionization_reactions import (
    reaction_is_ionization,
    narrow_down_ionization_type,
    determine_broad_ionization_tag
)
from classify_chemical_reactions import determine_chemical_reaction_tag
from reaction_classification_utilities import (
    add_reaction_if_new,
    write_reaction_classification
)
import sys
sys.path.append('../common')
from Rxn_classes import HiPRGen_Reaction

__version__ = '1.1.0'


def process_P2_reaction(rxn_dict, rxns_for_simulation, rxns_already_added):
    """
    Takes a reaction dictionary from phase 2, creates a HiPRGen reaction object
    from it, checks to see if it's new, and if it is, adds it to the list of
    reactions for simulation and adds its name to the set of added names

    Parameters
    ----------
    rxn_dict : dict
        dictionary of the form: {"reactants":[mpculeids],
                                 "proudcts":[mpculeids]}
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

    new_rxn = HiPRGen_Reaction(rxn_dict, phase=2)
    tag = determine_chemical_reaction_tag(new_rxn)
    reaction = add_reaction_if_new(
        new_rxn, tag, rxns_for_simulation, rxns_already_added
    )
    return reaction


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
    sorted_list = sorted(
        pathways_list, key=lambda x: (x['weight'], x['frequency'])
    )

    return sorted_list


def add_high_frequency_P2_reactions(
        directory, rxns_for_simulation, rxns_already_added,
        frequency_threshold):

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

        print(
            f"Error: The directory {directory} or the file does not exist."
        )

    except Exception as e:

        print(f"An unexpected error occurred: {e}")


def add_value_to_nested_dict(value, keys, nested_dict):
    """
    Adds a value to a nested dictionary, creating dictionaries and lists as
    needed.

    Parameters:
    value (Any): The value to add.
    keys (list): A list of keys representing the path in the nested dictionary.
    nested_dict (dict): The nested dictionary to update.
    """
    current_level = nested_dict
    # Traverse the dictionary up to the second to last key
    for key in keys[:-1]:
        if key not in current_level:
            current_level[key] = {}
        current_level = current_level[key]

    # Handle the last key
    last_key = keys[-1]
    if last_key not in current_level:
        current_level[last_key] = []

    current_level[last_key].append(value)


def add_reaction_to_dictionary(reaction, reaction_dict):
    classification_list = reaction.classification_list
    add_value_to_nested_dict(reaction, classification_list, reaction_dict)
    return reaction_dict


# # add all P1 reactions
# # P1_directory = "G:/My Drive/CRNs/071924_test_p1"
# P1_directory = "G:/My Drive/CRNs/euvl_092024_phase1"
# os.chdir(P1_directory)
# P1_pathways_and_reactions = loadfn("reaction_tally.json")
# P1_rxn_dicts = P1_pathways_and_reactions["reactions"].values()
rxns_already_added = set()
rxns_for_simulation = []

# for rxn_dict in P1_rxn_dicts:

#     new_rxn = HiPRGen_Reaction(rxn_dict, phase=1)

#     if reaction_is_ionization(new_rxn):

#         tag = determine_broad_ionization_tag(new_rxn)

#     else:

#         tag = determine_chemical_reaction_tag(new_rxn)

#     rxns_for_simulation, rxns_already_added = add_reaction_if_new(
#             new_rxn, tag, rxns_for_simulation, rxns_already_added
#         )

# for reaction in rxns_for_simulation:  # can narrow down now

#     if reaction.tag == "attachment_or_recombination":

#         reaction.tag = narrow_down_ionization_type(
#             reaction, rxns_for_simulation
#         )


# # add all P2 reactions that fired >= 500 times

# # P2_directory = "G:/My Drive/CRNs/071924_test_p2"
# P2_directory = "G:/My Drive/CRNs/euvl_092024_phase2"

# frequency_threshold = 100

# add_high_frequency_P2_reactions(
#     P2_directory, rxns_for_simulation, rxns_already_added, frequency_threshold
#     )

# network_products = loadfn("sink_report.json")

# for product_dict in network_products.values():

#     species_index = product_dict["species_index"]

#     reactions_and_pathways = loadfn(str(species_index) + "_pathway.json")
#     all_pathways_list = list(reactions_and_pathways["pathways"])
#     all_reactions = reactions_and_pathways["reactions"]
#     sorted_pathways_list = sort_pathways_list(all_pathways_list)
#     to_save_pathways = []

#     for pathway_dict in sorted_pathways_list:

#         pathway = pathway_dict["pathway"]
#         to_save_pathways.append(pathway)
#         number_pathways_saved = len(to_save_pathways)

#         if number_pathways_saved >= 10:
#             break

#     for pathway in to_save_pathways:
#         for reaction in pathway:
#             rxn_dict = all_reactions.get(str(reaction), None)
#             rxns_for_simulation, rxns_already_added = process_P2_reaction(
#                             rxn_dict, rxns_for_simulation, rxns_already_added
#                         )

# add recombination recations

recombination_directory = "G:/My Drive/Kinetiscope/production_simulations_092124/"
os.chdir(recombination_directory)
recombination_reactions = loadfn("recombination_reactions_121724.json")
for recomb_reaction in recombination_reactions:
    rxns_for_simulation, rxns_already_added = process_P2_reaction(recomb_reaction, rxns_for_simulation, rxns_already_added)

tagged_rxn_dict = {
    "ionization": {},
    "chemical": {
        "unimolecular": {
            "fragmentation": {},
            "isomerization": {}},
        "bimolecular": {
            "combination": {},
            "biproduct": {}
            }
        }
    }

ionization_classifications = set(["positive_ionization", "electron_attachment",
                                  "electron_cation_recombination"])
tagged_rxn_dict = {}

for reaction in rxns_for_simulation:

    reaction.classification_list = write_reaction_classification(reaction)
    add_reaction_to_dictionary(reaction, tagged_rxn_dict)

json_name = "Recombination_rxns_121824.json"
dumpfn(tagged_rxn_dict, json_name)
