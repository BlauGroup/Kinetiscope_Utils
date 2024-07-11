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
from reaction_classification_utilities import add_reaction_if_new
import time
from monty.serialization import loadfn, dumpfn
import copy
import operator
import pickle
import sys

def process_P2_reaction(rxn_dict, rxns_for_simulation, rxns_already_added):
    new_rxn = HiPRGen_Reaction(rxn_dict,phase=2)
    tag = determine_chemical_reaction_tag(new_rxn)
    return add_reaction_if_new(new_rxn, tag, rxns_for_simulation, rxns_already_added)

def sort_pathways_list(pathways_list):
    sorted_list = sorted(pathways_list, key=lambda x: (x['weight'], x['frequency']))
    return sorted_list
  
# add all P1 reactions
start=time.time()

os.chdir("G:/My Drive/CRNs/041323_p1")
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
        
    rxns_already_added, rxns_already_added = \
        add_reaction_if_new(new_rxn, tag, rxns_for_simulation, rxns_already_added)
        
for reaction in rxns_for_simulation: #can narrow down here when we have all
                                      #reactions dealt with
    
    if reaction.tag == "attachment_or_recombination":
        
        reaction.tag = narrow_down_ionization_type(reaction, rxns_for_simulation)
 
#add all P2 reactions that fired >= 500 times

os.chdir("G:/My Drive/CRNs/041423_p2")  

P2_pathways_and_reactions = loadfn("reaction_tally.json")
P2_index_frequency_dict =  P2_pathways_and_reactions["pathways"]
P2_index_rxn_dict = P2_pathways_and_reactions["reactions"]

for rxn_index, rxn_dict in P2_index_rxn_dict.items():
    
    rxn_frequency = P2_index_frequency_dict.get(rxn_index, None)
    
    if rxn_frequency >= 500:
        
        rxns_already_added, rxns_already_added = \
            process_P2_reaction(rxn_dict, rxns_for_simulation, rxns_already_added)

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
            rxns_already_added, rxns_already_added = \
                        process_P2_reaction(rxn_dict, rxns_for_simulation, rxns_already_added)
                        
# tagged_rxn_dict = {}

# for reaction in rxns_for_simulation:
                
#     if reaction.tag not in tagged_rxn_dict:
        
#         tagged_rxn_dict[reaction.tag] = []
    
#     tagged_rxn_dict[reaction.tag].append(reaction)
    
# dumpfn(tagged_rxn_dict,"HiPRGen_rxns_to_name.json")

# def print_dict_lengths(dictionary):
#     for key, value in dictionary.items():
#         value_length = len(value)
#         print(f"Key '{key}' has a value with length {value_length}")

# print_dict_lengths(tagged_rxn_dict)