# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:54:29 2024

@author: JRMilton
"""
import os
from monty.serialization import loadfn, dumpfn
from Rxn_classes import HiPRGen_Reaction
from classify_ionization_reactions import (
    process_ionization_reactions,
    narrow_down_ionization_type
)
from classify_chemical_reactions import tag_chemical_reaction

os.chdir("G:/My Drive/CRNs/041323_p1")
P1_rxn_tally = loadfn("reaction_tally.json")
all_P1_reactions = P1_rxn_tally["reactions"].values()
P1_rxn_hashes = set()
ionization_hashes = set()
ionization_rxn_list = []

for rxn_dict in all_P1_reactions:
    
    new_rxn = HiPRGen_Reaction(rxn_dict,phase=1)
    
    #we didn't ask for TSs for ionzation rxns but need them in our simulations
    #so we save them here
    
    ionization_rxn_list, ionization_hashes = \
    process_ionization_reactions(new_rxn, ionization_rxn_list, ionization_hashes)
    
    P1_rxn_hashes.add(new_rxn.reaction_hash) #lets us test if other reactions 
                                             #are in P1 or P2 later

for reaction in ionization_rxn_list:
    
    if reaction.tag == "attachment_or_recombination":
        
        reaction.tag = narrow_down_ionization_type(reaction, ionization_rxn_list)
            
#read in HiPRGen rxns we want to write names for

os.chdir("G:/My Drive/CRNs/Kinetiscope_rxn_naming052824")
all_requested_reactions = loadfn("euvl_TSreactions_041823.json") #list of all 
                                                                 #reactions we
                                                                 #asked for TSs

rxns_for_simulation = ionization_rxn_list
rxns_already_added = ionization_hashes

for rxn_dict in all_requested_reactions:
    
    new_rxn = HiPRGen_Reaction(rxn_dict,tag="chemical")
    new_rxn.tag = tag_chemical_reaction(new_rxn)
    
    new_rxn.phase = 1 if new_rxn.reaction_hash in P1_rxn_hashes else 2
    
    if new_rxn.reaction_hash not in rxns_already_added:
        
        rxns_for_simulation.append(new_rxn)
        rxns_already_added.add(new_rxn.reaction_hash)
        
tagged_rxn_dict = {}

for reaction in rxns_for_simulation:
                
    if reaction.tag not in tagged_rxn_dict:
        
        tagged_rxn_dict[reaction.tag] = []
    
    tagged_rxn_dict[reaction.tag].append(reaction)
    
dumpfn(tagged_rxn_dict,"HiPRGen_rxns_to_name.json")

def print_dict_lengths(dictionary):
    for key, value in dictionary.items():
        value_length = len(value)
        print(f"Key '{key}' has a value with length {value_length}")


print_dict_lengths(tagged_rxn_dict)