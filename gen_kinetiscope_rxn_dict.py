# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:54:29 2024

@author: JRMilton
"""
__version__ = '1.1.0'

import os
import sys
from monty.serialization import loadfn, dumpfn
from Rxn_classes import HiPRGen_Reaction, Kinetiscope_Reaction
import sys
from write_kinetiscope_ionization_reactions import handle_ionization_reactions
from write_kinetiscope_chemical_reactions import handle_chemical_rxns
# from reaction_classification_utilities import find_mpculeid_charge

# def build_rxn_object(HiPRGen_rxn, kinetiscope_name, rate_constant, order, marker_species):
#     return Kinetiscope_Reaction(HiPRGen_rxn, kinetiscope_name, rate_constant, order, marker_species)

def collect_lists_from_nested_dict(d):
    collected_items = []
    
    def recurse_through_dict(d):
        for value in d.values():
            if isinstance(value, dict):
                # If the value is a dictionary, recurse into it
                recurse_through_dict(value)
            elif isinstance(value, list):
                # If the value is a list, extend the collected_items list with the items
                collected_items.extend(value)

    recurse_through_dict(d)
    return collected_items
              
#builds expected number of ionization reactions

kinetiscope_reaction_list = []
os.chdir("G:/My Drive/Kinetiscope/new_kinetiscope_naming_080224")
HiPRGen_reaction_list = loadfn("HiPRGen_rxns_to_name.json")
name_mpculeid_file = "name_test_mpculeid_080624.json"
name_mpculeid_dict = loadfn(name_mpculeid_file)
mpculeid_name_dict = {mpculeid: name for name, mpculeid in name_mpculeid_dict.items()}

rate_constant_dict = { #values calculated in excel and taken as scientific notation values
    "ionization":
        {"absorption":{
            "0aade5ee5263fd1ad77490688fb37d0e-C10H20O2-0-1":1.4e-01,
            "1d64d1a96b5db50b9fdf583dc18f90cf-C10H14O1-0-1":1.0e-01,
            "94be97269b03e361ba1b344f48f19e44-C18H15S1-1-1":1.8e-01,
            "f357352f5dd5488b54d50242d228ed6d-C4F9O3S1-m1-1":5.8e-01},
            # "4-CNBZ":1.5E-01},
         "electron_ionization":{
             "eV_80":1.4e+15,
             "eV_55":1.2e+15,
             "eV_30":8.5e+14,
             "LEE":1.1E+14},
         "recombination":{
             "all": 4.1e+13},
         "attachment":{
             "positive":4.1e+13,
             "neutral":3.1e+12}
         },
    "chemical":
        {"2nd_order":6.2e+12,
         "1st_order":1.3e+12}
        }
    
marker_species_shorthand = {
    "radical_cation":"rc",
    "radical_anion":"ra",
    "neutral_radical":"nr",
    "cation":"c",
    "anion":"a",
    # "radical":"r",
    "neutral":"n",
    "proton_coupled_electron_transfer":"PCET"}


excitation_set = set()
# LEE_H_rxn = None
# LEE_collision_name = "LEE => TE"
# LEE_rate_constant = 5.1E+14
# LEE_order = 1
# LEE_marker_species = None

# LEE_rxn = \
#     build_rxn_object(LEE_H_rxn, LEE_collision_name, LEE_rate_constant, LEE_order, LEE_marker_species)

# for rxn_list in HiPRGen_reaction_list["ionization"].values():
#     for HiPRGen_rxn in rxn_list:
#         ionization_reaction_list = \
#             handle_ionization_reactions(HiPRGen_rxn, mpculeid_name_dict, rate_constant_dict)
#         kinetiscope_reaction_list.extend(ionization_reaction_list)

chemical_reaction_list = collect_lists_from_nested_dict(HiPRGen_reaction_list["chemical"])
for H_rxn in chemical_reaction_list:
    if H_rxn.phase == 1:
        chemical_reaction_list, excitation_set = handle_chemical_rxns(H_rxn, mpculeid_name_dict, rate_constant_dict, marker_species_shorthand, excitation_set)
        kinetiscope_reaction_list.extend(chemical_reaction_list)
        
            
# supercategory_order = [
#     "absorption", "electron_ionization", "recombination", "attachment", 
#      "fragmentation", "isomerization", "ion-ion", "ion-molecule", "neutral"]

# supercategoeries_with_subcategories = set(["ion-ion", "ion-molecule", "neutral"])

# subcategory_order = [
#     "proton_transfer", "H_atom_abstraction", "hydride_abstraction", 
#     "proton_coupled_electron_transfer", "electron_transfer", "reaction"
#     ] 

# dumpfn(kinetiscope_reaction_list, "kinetiscope_reaction_list.json")
        
