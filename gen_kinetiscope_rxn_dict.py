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
    
    name_without_marker_species = "".join(reaction_as_list)
    name_with_marker_species = add_marker_species(name_without_marker_species)
    
    return name_with_marker_species  

def build_ionization_reactions(H_reaction, mpculeid_dict, rate_constant_dict, ionization_type, order, specific_list=None):
    reaction_list = []
    marker_species_list = H_reaction.classification_list[:-1]  # does not include tag

    if specific_list is None:
        specific_list = [H_reaction.reactants[0]]

    for item in specific_list:
        K_marker_species_list = marker_species_list + [ionization_type]
        kinetiscope_name = write_kinetiscope_name(H_reaction.name, mpculeid_dict, K_marker_species_list)
        rate_constant = rate_constant_dict["ionization"][ionization_type].get(item, None)
        reaction_list.append(build_rxn_object(H_reaction, kinetiscope_name, rate_constant, order, K_marker_species_list))
    
    return reaction_list if len(reaction_list) > 1 else reaction_list[0]

def build_light_ionization(H_reaction, mpculeid_dict, rate_constant_dict):
    return build_ionization_reactions(H_reaction, mpculeid_dict, rate_constant_dict, "absorption", 1)

def build_electron_ionization_reactions(H_reaction, mpculeid_dict, rate_constant_dict):
    electron_list = ["eV_80", "eV_55", "eV_30"]
    return build_ionization_reactions(H_reaction, mpculeid_dict, rate_constant_dict, "electron_ionization", 2, electron_list)

def build_positive_ionization_reactions(H_reaction, mpculeid_dict, rate_constant_dict):
    positive_ionization_reaction_list = []
    potential_absorber = H_reaction.reactants[0]
    absorbers = rate_constant_dict["ionization"]["absorption"].keys()
    if potential_absorber in absorbers:
        light_ionization_reaction = \
            build_light_ionization(H_reaction, mpculeid_dict, rate_constant_dict)
        positive_ionization_reaction_list.append(light_ionization_reaction)
    electron_ionization_reactions = \
        build_electron_ionization_reactions(H_reaction, mpculeid_dict, rate_constant_dict)
    return positive_ionization_reaction_list + electron_ionization_reactions
    
# def build_light_ionization(H_reaction, mpculeid_dict, rate_constant_dict):
#     absorber = H_reaction.reactants[0]
#     H_marker_species_list = H_reaction.classification_list[:-1] #does not include tag
#     K_marker_species_list =  H_marker_species_list.append("absorption")
#     kinetiscope_name = \
#         write_kinetiscope_name(H_reaction.name, mpculeid_dict, K_marker_species_list)
#     rate_constant = rate_constant_dict["ionization"]["absorption"].get(absorber, None)
#     order = 1
#     return build_rxn_object(H_reaction, kinetiscope_name, rate_constant, order, K_marker_species_list)
    
# def build_electron_ionization_reactions(H_reaction, mpculeid_dict, rate_constant_dict):
#     electron_list = ["eV_80", "eV_55", "eV_30"]
#     reaction_list = []
#     for electron in electron_list:
#         # reactant = H_reaction.reactants[0]
#         H_marker_species_list = H_reaction.classification_list[:-1] #does not include tag
#         K_marker_species_list =  H_marker_species_list.append("electron_ionization")
#         kinetiscope_name = \
#             write_kinetiscope_name(H_reaction.name, mpculeid_dict, K_marker_species_list)
#         rate_constant = \
#             rate_constant_dict["ionization"]["electron"].get(electron, None)
#         order = 2
#         reaction_list.append(build_rxn_object(H_reaction, kinetiscope_name, rate_constant, order, K_marker_species_list))
#     return reaction_list
    


def build_attachment_or_recombination():
    ionization_type = "recombination" if "recombination" in H_rxn.tag else "attachment"
    
    H_reaction, mpculeid_dict, rate_constant_dict, ionization_type, order, specific_list=None

def build_excitation_reactions():
    pass

# def write_recombination():
#     pass

def handle_ionization_reactions(H_reaction, mpculeid_dict, rate_constant_dict):
    if HiPRGen_rxn.tag == "positive_ionization":
        return build_positive_ionization_reactions(H_reaction, mpculeid_dict, rate_constant_dict)
    return build_attachment_or_recombination(H_reaction, mpculeid_dict, rate_constant_dict)

def determine_order():
    pass

def handle_chemical_reactions():
    pass

def attach_rate_constant(rate_constant_dict):
    pass

def build_rxn_object(HiPRGen_rxn, kinetiscope_name, rate_coefficient, order, marker_species):
    return Kinetiscope_Reaction(HiPRGen_rxn, kinetiscope_name, rate_coefficient, order, marker_species)


kinetiscope_reaction_list = []
os.chdir("G:/My Drive/Kinetiscope/new_kinetiscope_naming_080224")
HiPRGen_reaction_list = loadfn("HiPRGen_rxns_to_name.json")

mpculeid_easy_name_dict = {
    "0aade5ee5263fd1ad77490688fb37d0e-C10H20O2-0-1": "PtBMA",
    "1d64d1a96b5db50b9fdf583dc18f90cf-C10H14O1-0-1": "PHS",
    "94be97269b03e361ba1b344f48f19e44-C18H15S1-1-1": "TPS",
    "f357352f5dd5488b54d50242d228ed6d-C4F9O3S1-m1-1": "Nf",
    # "": "4-CNBZ" not in test set
}

rate_constant_dict = { #values calculated in excel and taken as scientific notation values
    "ionization":
        {"absorption":{
            "0aade5ee5263fd1ad77490688fb37d0e-C10H20O2-0-1":1.4E-01,
            "1d64d1a96b5db50b9fdf583dc18f90cf-C10H14O1-0-1":1.0E-01,
            "94be97269b03e361ba1b344f48f19e44-C18H15S1-1-1":1.8E-01,
            "f357352f5dd5488b54d50242d228ed6d-C4F9O3S1-m1-1":5.8E-01},
            # "4-CNBZ":1.5E-01},
         "electron":{
             "eV_80":1.4E+15,
             "eV_55":1.2E+15,
             "eV_30":8.5E+14},
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
        ionization_reaction_list = handle_ionization_reactions()
        kinetiscope_reaction_list = kinetiscope_reaction_list + ionization_reaction_list
