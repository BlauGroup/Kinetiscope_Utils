# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:54:29 2024

@author: JRMilton
"""
import os
import sys
from monty.serialization import loadfn
from Rxn_classes import HiPRGen_Reaction, Kinetiscope_Reaction

def write_kinetiscope_name():
    pass

def attach_rate_constant():
    pass

def build_rxn_objects():
    pass

os.chdir("G:/My Drive/CRNs/Kinetiscope_rxn_naming052824")
HiPRGen_reaction_list = loadfn("HiPRGen_rxns_to_name.json")
print(HiPRGen_reaction_list)
#build Kinetiscope dict from HiPRGen dict

kinetiscope_dict = {
    "ionization":{},
    "recombination":{},
    "chemical":{
        "ion-ion":[],
        "aromatic_sub":[],
        "H-rxn":{
            "proton":[],
            "h_atom":[],
            "hydride":[]
        },
        "ion-molecule":[]
    }
}