# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 12:55:07 2024

@author: jacob
"""
import pickle
from monty.serialization import loadfn
import os
import sys
import networkx as nx
sys.path.append('../analyze_kinetiscope_data')
from associate_indicies_with_frequencies import build_index_freq_dict
from create_kinetiscope_simulation_reactions import build_index_reaction_dict
sys.path.append('../common')
from utilities import correct_path_change_dir
from Rxn_classes import HiPRGen_Reaction
sys.path.append('../Classify_HiPRGen_reactions')
from reaction_classification_utilities import (
    find_mpculeid_charge,
    find_mpculeid_spin
)

def find_relevant_reactions():
    pass


def find_kscope_rxns(select_freq_file, reaction_name_file, start_index, end_index):
    index_freq_dict = (
        build_index_freq_dict(select_freq_file, start_index, end_index)
    )
    
    index_reaction_dict = (
        build_index_reaction_dict(reaction_name_file, index_freq_dict)
    )

    return index_reaction_dict


def extract_radical_mpculeids():
    pass


def collect_lists_from_nested_dict(d):
    """
    Collect all list items from a nested dictionary.

    This function traverses a nested dictionary and collects all items
    from lists contained within the dictionary. It recurses through nested
    dictionaries and aggregates items from all lists found.

    Parameters:
    -----------
    d : dict
        The nested dictionary from which lists are to be collected.

    Returns:
    --------
    list
        A list containing all items collected from lists within the nested
        dictionary.
    """
    def recurse_through_dict(d, collected_items):
        """
        Recursively traverse a nested dictionary to collect items from lists.

        This helper function is used by `collect_lists_from_nested_dict` to
        perform a recursive traversal of the dictionary. It extends the list
        of collected items with elements from any lists found in the dictionary.

        Parameters:
        -----------
        d : dict
            The nested dictionary to traverse.
        collected_items : list
            The list to which items from nested lists are added.
        """

        for value in d.values():

            if isinstance(value, dict):

                recurse_through_dict(value, collected_items)

            elif isinstance(value, list):

                collected_items.extend(value)

    collected_items = []
    recurse_through_dict(d, collected_items)

    return collected_items

def build_radicals_sets(mol_entries, relevant_radical_set):
    def find_starting_species_graphs(mol_entries):
        starting_species_mpculeids = {
            "TPS": "bfed458e642b8daa1eab6bc02d5e5682-C18H15S1-1-1",
            "Nf": "9a8a88b8b92c714d7f65b8526ffabc7a-C4F9O3S1-m1-1",
            "4-CNBZ": "17f31f89123edbaa0e3b9c7eb49d26f3-C8H4N1O2-m1-1",
            "PHS": "00a7dcc352b0d613f58e850935bf5609-C10H14O1-0-1",
            "PtBMA": "4864aee73a83d357c31fadd50b81e3cd-C10H20O2-0-1"
        }
        starting_species_graphs = []
        
        for mol_entry in mol_entries:
            if mol_entry.entry_id in starting_species_mpculeids.values():
                starting_species_graphs.append(mol_entry.graph)
        return starting_species_graphs
    
    def is_neutral_radical(mol_entry):
        return mol_entry.spin_multiplicity == 2 and mol_entry.charge == 0
    
    def radical_is_relevant(mol_entry, relevant_radical_set):
        return mol_entry.entry_id in relevant_radical_set
    
    def radical_isnt_fragment(mol_entry, starting_species_graphs):
        nm = nx.isomorphism.categorical_node_match("specie", None) # ensures isomorphic graphs must have the same atoms
        for graph in starting_species_graphs:
            isomorphism = nx.isomorphism.GraphMatcher(mol_entry.graph, graph, node_match = nm)
            if isomorphism.subgraph_is_isomorphic():
                return False
        return True

    all_radicals = []
    all_pi_radicals = []
    potentially_problematic = []

    starting_species_graphs = find_starting_species_graphs(mol_entries)

    for mol_entry in mol_entries:

        if is_neutral_radical(mol_entry) and radical_is_relevant(mol_entry, relevant_radical_set):

            all_radicals.append(mol_entry)
            if radical_isnt_fragment(mol_entry, starting_species_graphs):
                potentially_problematic.append(mol_entry)
    total_number_radicals = len(all_radicals)
    number_problematic_radicals = len(potentially_problematic)
    print(f"total number of radicals: {total_number_radicals}")
    print(f"total number of problematic radicals: {number_problematic_radicals}")
    print(f"estimated number of reactions: {total_number_radicals * number_problematic_radicals}")
    sys.exit()

            

    # for mol_entry in mol_entries:
    #     if 

# mongodb_info = {
#     "database": "sb_qchem",
#     "collection": "new_tasks",
#     "admin_user": "smblau_lbl.gov_readWrite",
#     "admin_password": "baffler-underranger-sanguinely-distent-flukeworm",
#     "host": "mongodb03.nersc.gov",
#     "port": 27017,
#     "aliases": {},
#     "authSource": "sb_qchem"
# }

# client = pymongo.MongoClient( #connect to mongo db
#     host=mongodb_info["host"], username = mongodb_info["admin_user"],
#     password = mongodb_info["admin_password"], authSource = mongodb_info["authSource"])

# database = client.sb_qchem
# collection = database.euvl_mar_summary

# load current names and mpculeids 

os.chdir("C:/Users/jacob/Kinetiscope_Utils/classify_HiPRGen_reactions")
current_mpculeids_and_names = loadfn("name_full_mpculeid_092124.json")

print("Loading mol entries...")

with open('mol_entries.pickle', 'rb') as f: #loads HiPRGen mol_entry objs
        
    try:
        
        mol_entries = pickle.load(f)
    
    except AttributeError:
        
        print("Need to use the pickle file from phase 2, not phase 1")
        print("The 'electron species' in phase 1 causes issues here")
        sys.exit()
    
print('Done!')

# full_rxns = "HiPRGen_rxns_to_name_full_092124.json"
# HiPRGen_reaction_list = loadfn(full_rxns)
# HiPRGen_ionization_reactions = list(HiPRGen_reaction_list["ionization"].values())[0]
# chemical_reaction_list = (
#     collect_lists_from_nested_dict(HiPRGen_reaction_list["chemical"])
# )

# all_HiPRGen_rxns = HiPRGen_ionization_reactions + chemical_reaction_list

# relevant_radical_set = set()

# for HiPRGen_rxn in all_HiPRGen_rxns:
#     reactants = HiPRGen_rxn.reactants
#     products = HiPRGen_rxn.products
    
#     species = reactants + products
    
#     for specie in species:
#         charge = find_mpculeid_charge(specie)
#         spin = find_mpculeid_spin(specie)
        
#         if charge == 0 and spin == 2:
#             relevant_radical_set.add(specie)

# make a set of only radicals that actually are involved in HiPRGen reactions 

# built radical and recomb_radical sets

kinetiscope_files_dir = (
    r"G:\My Drive\Kinetiscope\production_simulations_092124"
)

correct_path_change_dir(kinetiscope_files_dir)

select_freq_file = "excitation_selection_freq_092524.txt"
reaction_name_file = "excitation_reactions_092524.txt"

kscope_rxns = find_kscope_rxns(
    select_freq_file,
    reaction_name_file,
    start_index=7,
    end_index=5384
)

relevant_radical_list = []

for kscope_rxn in kscope_rxns.values():
    if kscope_rxn.selection_freq > 0:
        species_list = kscope_rxn.kinetiscope_name.split()
        chemical_dict = loadfn("name_full_mpculeid_092124.json")
        for species in species_list:
            mpculeid = chemical_dict.get(species, None)
            if mpculeid:
                charge = find_mpculeid_charge(mpculeid)
                spin = find_mpculeid_spin(mpculeid)
                if charge == 0 and spin == 2:
                    radical_dict = {mpculeid:kscope_rxn}
                    if mpculeid not in relevant_radical_list:
                        relevant_radical_list.append(radical_dict)
                    elif 

# radical_reactions = list()

# for kscope_rxn in kscope_rxns.values():
#     if kscope_rxn.selection_freq > 0:
#         species_list = kscope_rxn.kinetiscope_name.split()
#         chemical_dict = loadfn("name_full_mpculeid_092124.json")
#         for species in species_list:
#             mpculeid = chemical_dict.get(species, None)
#             if mpculeid:
#                 charge = find_mpculeid_charge(mpculeid)
#                 spin = find_mpculeid_spin(mpculeid)
#                 if charge == 0 and spin == 2:
#                     if kscope_rxn not in radical_reactions:
#                         radical_reactions.append(kscope_rxn)

# sorted_reactions = sorted(radical_reactions, key=lambda reaction: reaction.selection_freq, reverse=True)

# sum = 0

# for reaction in sorted_reactions:
#     sum += reaction.selection_freq
# print(sum)
# print(0.01*sum)

# all_radicals, radicals_to_recombine = build_radicals_sets(
#     mol_entries, relevant_radical_set
# )

# for each recomb radical,
# generate a recomb rxn,
#  save the recomb mpculid and
# name
# if nbo_spins[atom_ind] > 0.09 and str(mol_graph.molecule[int(atom_ind)].specie) != "H": #how did we choose 0.09 as the cutoff?