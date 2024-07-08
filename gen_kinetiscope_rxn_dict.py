# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:54:29 2024

@author: JRMilton
"""
import os
import sys
from monty.serialization import loadfn
from Rxn_classes import HiPRGen_Reaction, Kinetiscope_Reaction
from build_kinetiscope_rxns import write_light_ionization

os.chdir("G:/My Drive/CRNs/Kinetiscope_rxn_naming052824")
HiPRGen_reaction_dict = loadfn("HiPRGen_rxns_to_name.json")
name_mpculeid_dict = loadfn("name_mpculeid_association_060424.json")
mpculeid_name_dict = {value:key for key, value in name_mpculeid_dict.items()}

for reaction in HiPRGen_reaction_dict["positive_ionization"]:
    print(write_light_ionization(reaction, mpculeid_name_dict))

# kinetiscope_dict = {
#     "ionization":{
#         "photon:":[],
#         "electron":[],
#         "attachment":[]
#         },
#     "recombination":[],
#     "chemical":{
#         "H-rxn":{
#             "proton":[],
#             "h_atom":[],
#             "hydride":[]
#         },
#         "ion-ion":[],
#         "aromatic_sub":[],
#         "ion-molecule":[],
#         "neutral-radical":[]
#     }
# }