# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:54:29 2024

@author: JRMilton
"""
import os
import sys
from monty.serialization import loadfn
from Rxn_classes import HiPRGen_Reaction, Kinetiscope_Reaction

#TODO account for 1st vs. second-order rate constants when using kT/h
def write_light_ionization():
    pass

def write_electron_ionization():
    pass

def write_positive_ionization_reactions():
    pass

def write_electron_attachment():
    pass

def write_recombination():
    pass

def handle_ionization_reactions():
    pass

def determine_order():
    pass

def handle_chemical_reactions():
    pass

def write_kinetiscope_name():
    pass

def attach_rate_constant():
    pass

def build_rxn_objects():
    pass

os.chdir("G:/My Drive/CRNs/Kinetiscope_rxn_naming052824")
HiPRGen_reaction_list = loadfn("HiPRGen_rxns_to_name.json")

kinetiscope_dict = {
    "ionization":{
        "photon:":[],
        "electron":[],
        "attachment":[]
        },
    "recombination":[],
    "chemical":{
        "H-rxn":{
            "proton":[],
            "h_atom":[],
            "hydride":[]
        },
        "ion-ion":[],
        "aromatic_sub":[],
        "ion-molecule":[],
        "neutral-radical":[]
    }
}