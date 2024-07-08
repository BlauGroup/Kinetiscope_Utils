# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:54:29 2024

@author: JRMilton
"""
import os
import sys
from monty.serialization import loadfn
from Rxn_classes import HiPRGen_Reaction, Kinetiscope_Reaction

def write_light_ionization(HiPRGen_rxn, mpculeid_name_dict):
    reactant = HiPRGen_rxn.reactants[0]
    print(reactant)
    reactant_name = mpculeid_name_dict.get(reactant, None)
    print(reactant_name)
    product = HiPRGen_rxn.products[0]
    product_name = mpculeid_name_dict.get(product, None)
    return reactant_name + " => " + product_name + " + eV_80"

# def calc_molecular_cross_section(reactant_formula, cross_section_dict):
#     pass

# def calc_light_rate_constant(reactant_mpculeid, concentration_dict):
#     reactant_concentration = concentration_dict.find(reactant_mpculeid, 0)
#     if reactant_concentration:
#         reactant_formula = reactant_mpculeid.split("-")[1]
#         molecular_cross_section = calc_molecular_cross_section(reactant_formula)
#         cross_sectional_area = molecular_cross_section * reactant_concentration
#         total_volume = 3.00e-14
        
    

# def build_light_reaction(HiPRGen_rxn, mpculeid_name_dict, light_constant_dict):
#     name = write_light_ionization(HiPRGen_rxn, mpculeid_name_dict)
#     absorber = HiPRGen_rxn.reactants[0]
#     abs_rate_constant = light_constant_dict.find(absorber, None)
#     order = 1
#     tag = "absorption"
#     kinetiscope_rxn = Kinetiscope_Reaction(name, abs_rate_constant, order, tag)
#     return kinetiscope_rxn
    
# def write_electron_ionization():
#     pass

# def build_electron_ionization():
#     electrons = ["ev_80", "eV_55", "eV_30", "LEE", "TE"]
    
# def write_electron_attachment():
#     pass

# def write_recombination():
#     pass

# def write_ionization_reactions():
#     pass

# def write_kinetiscope_name():
#     pass

# def attach_rate_constant():
#     pass

# def build_rxn_objects():
#     pass
