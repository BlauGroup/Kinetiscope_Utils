# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:54:29 2024

@author: JRMilton
"""
import os
import sys
from monty.serialization import loadfn
from Rxn_classes import HiPRGen_Reaction, Kinetiscope_Reaction
# from openpyxl import load_workbook

def add_marker_species(name_without_marker_species, marker_species):
    name_with_marker_species = name_without_marker_species
    
    for marker_species in marker_species:
        name_with_marker_species = \
            name_with_marker_species + "+ " + marker_species
    
    return name_with_marker_species
    
def write_kinetiscope_name(HiPRGen_name, mpculeid_dict, marker_species):
    reaction_as_list = HiPRGen_name.split()
    
    for index, obj in reaction_as_list[:]:
        kinetiscope_name = mpculeid_dict.get(obj, None)
        if kinetiscope_name:
            reaction_as_list[index] = kinetiscope_name
    
    name_without_marker_species = reaction_as_list.join()
    name_with_marker_species = add_marker_species(name_without_marker_species)
    return name_with_marker_species        
    
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

def build_rxn_object():
    pass

# def open_book_change_sheet(book, sheet):
#     workbook = load_workbook(book)
#     rate_constant_sheet = workbook[sheet]
#     return rate_constant_sheet

# def populate_rate_constant_dict(rate_constant_dict, rate_constant_book, rate_constant_sheet):
#     rate_constant_sheet = \
#         open_book_change_sheet(rate_constant_book, rate_constant_sheet)
        

kinetiscope_reaction_list = []
os.chdir("G:/My Drive/Kinetiscope/new_kinetiscope_naming_080224")
HiPRGen_reaction_list = loadfn("HiPRGen_rxns_to_name.json")

rate_constant_dict = { #values calculated in excel and taken as scientific notation values
    "ionization":
        {"absorption":{
            "PtBMA":1.4E-01,
            "PHS":1.0E-01,
            "TPS":1.8E-01,
            "Nf":5.8E-01,
            "4-CNBZ":1.5E-01},
         "electron":{
             "80_eV":1.4E+15,
             "55_eV":1.2E+15,
             "30_eV":8.5E+14},
         "recombination":4.1E+13,
         "attachment":3.1E+12
         },
    "chemical":
        {"2nd_order":6.2E+12,
         "1st_order":1.3E+12}
        }

supercategory_order = [
    "absorption", "electron_ionization", "recombination", "attachment", 
     "fragmentation", "isomerization", "ion-ion", "ion-molecule", "neutral"]

supercategoeries_with_subcategories = set(["ion-ion", "ion-molecule", "neutral"])

subcategory_order = [
    "proton_transfer", "H_atom_abstraction", "hydride_abstraction", 
    "proton_coupled_electron_transfer", "electron_transfer", "reaction"
    ]    

for rxn_list in HiPRGen_reaction_list["ionization"]:
    for HiPRGen_rxn in rxn_list:
        if HiPRGen_rxn.tag == "positive_ionization":
            print(HiPRGen_rxn)


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