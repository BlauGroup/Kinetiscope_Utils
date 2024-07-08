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
from monty.serialization import loadfn, dumpfn
import copy
import operator
import pickle

# reaction_list = [] #list of reactions we want to save of the form: [{'reactants': [mpculid, mpculeid], 'products': etc.}]
# rxns_added = set() #keeps track of reactions we have already added to prevent redundancy
# added_hashes = {}
first_name = 'reaction_tally_p1'
first_entries = loadfn(first_name + ".json") #loads json as a dictionary whose keys are mol ids and values are
                                                #dictionaries whose keys are labels of values

# print('Loading mol_entries.pickle...')
# with open('mol_entries.pickle', 'rb') as f: #loads mol_entries from pickle file
#     mol_entries = pickle.load(f)
# print('Done!')

# # Create a mapping from molecule id to its properties for quick lookup
# mol_entry_dict = {entry.entry_id: entry for entry in mol_entries}
# reaction_indicies = first_entries["pathways"].keys()
# reaction_index_dict = first_entries["reactions"]
    
# print('Adding reactions from phase 1...')   
# for reaction_index in reaction_indicies:
#     reaction = reaction_index_dict.get(reaction_index, False)
#     assert bool(reaction) == True
#     reaction.pop(None, None) #removes "electron species" which has no mpculeid
#     # is_ionization_reaction = False
#     # for species in reaction.values():
#     #     if species == None:
            
#             # is_ionization_reaction = True
#     # if not is_ionization_reaction:
#     new_rxn = HiPRGen_Reaction(reaction,phase=1)
#     if new_rxn.reaction_hash not in rxns_added:
#         # reaction_dict = create_reaction_dict(reaction, mol_entry_dict)
#         # if not resonant_reaction(reaction_dict, added_hashes):
#         rxns_added.add(new_rxn.reaction_hash)
#         reaction_list.append(reaction)
#                 # added_hashes.update(reaction_dict.items())
            
#     # for rxn in first_entries["reactions"].keys():
#     #     if str(reaction) == rxn:
#              #removes reactions occuring between a reactant/product and an "electron"
#             # for entry in first_entries["reactions"][rxn].values():
#             #     for species in entry:
#             #         if species == None:
#             #             electron_test = True
#             # if not electron_test:
#                 # if reaction not in added: #don't add reactions twice
#                 #     reactants = first_entries["reactions"][rxn]["reactants"]
#                 #     products = first_entries["reactions"][rxn]["products"]
#                 #     participants = [reactants, products] 
                    
                    
                        
                        
print('Done! ', len(mpcule_ids), ' added')


print('Adding reactions for phase 2 network products...')                                             
second_name = 'sink_report'
second_entries = loadfn(second_name + ".json") 

for dictionary in second_entries.values(): #iterates through each network product
    species_index = dictionary["species_index"]
    reaction_json = loadfn(str(species_index) + "_pathway.json") #loads each network_product.json file as dictionary containing  
    pathways = reaction_json["pathways"]                         #two keys: pathways and reactions.
    pathways.sort(key = operator.itemgetter('frequency'), reverse = True) #sorts by frequency
    top_pathways = []
    n = 1
    ten_saved = False
    
    while not ten_saved: #makes sure we save top ten reactions even for network products that don't have 10 1-step reactions forming them
        for dictionary in pathways:
            if int(dictionary['weight']) == n:
                top_pathways.append(dictionary['pathway'])
                if len(top_pathways) == 10:
                    ten_saved = True
                    break
        n += 1
    
    for rxn in top_pathways:
        for reaction in rxn:
            for num in reaction_json["reactions"].keys():
                if num == str(reaction):
                    reactants = reaction_json["reactions"][num]["reactants"]
                    products = reaction_json["reactions"][num]["products"]
                    participants = [reactants, products] 
                    reaction_dict = create_reaction_dict(participants, mol_entry_dict)
                    if not resonant_reaction(reaction_dict, added_hashes):
                        if not reverse_reaction(reaction_dict, added_hashes):
                            if not charge_transfer_reaction(reaction_dict):
                                mpcule_ids.append(reaction_json["reactions"][num])
                                added.append(num)
                                added_hashes.update(reaction_dict.items())
    n = 1
                
print('Done! ', len(mpcule_ids), ' reactions total')

print('Adding reactions for phase 2 tally...')  
third_name = 'reaction_tally'
third_entries = loadfn(third_name + ".json")

for reaction in third_entries["pathways"].keys():
    if third_entries["pathways"][reaction] > 500: #only add network products found >500 times
        for rxn in third_entries["reactions"].keys():
            if str(reaction) == rxn:
                reactants = third_entries["reactions"][reaction]["reactants"]
                products = third_entries["reactions"][reaction]["products"]
                participants = [reactants, products]
                reaction_dict = create_reaction_dict(participants, mol_entry_dict)
                if not resonant_reaction(reaction_dict, added_hashes):
                    if not reverse_reaction(reaction_dict, added_hashes):
                        if not charge_transfer_reaction(reaction_dict):
                            mpcule_ids.append(third_entries["reactions"][reaction])
                            added.append(reaction)
                            added_hashes.update(reaction_dict.items())

print('Done! ', len(mpcule_ids), ' reactions total')
dumpfn(mpcule_ids, 'euvl_TSreactions_041823.json')

os.chdir("G:/My Drive/CRNs/041323_p1")
P1_rxn_tally = loadfn("reaction_tally.json")
all_P1_reactions = P1_rxn_tally["reactions"].values()
already_added = set()
rxns_for_simulation = []

for rxn_dict in all_P1_reactions:
    
    new_rxn = HiPRGen_Reaction(rxn_dict,phase=1)
    
    #we didn't ask for TSs for ionzation rxns but need them in our simulations
    #so we save them here
    
    ionization_rxn_list, ionization_hashes = \
    process_ionization_reactions(new_rxn, ionization_rxn_list, ionization_hashes)
    
    rxns_for_simulation.append(new_rxn)
    rxns_already_added.add(new_rxn.reaction_hash)
    
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