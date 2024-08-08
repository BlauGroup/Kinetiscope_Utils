# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 14:24:24 2024

@author: jacob
"""

from Rxn_classes import Kinetiscope_Reaction

def add_reactants(base_name, added_reactants):
    return " + ".join(added_reactants + [base_name])

def add_products(base_name, added_products):
    return " + ".join([base_name] + added_products)
    
def add_reacting_species(initial_name, added_reactants, added_products):
    name_after_reactants_alias = initial_name

    if added_reactants:
        name_after_reactants_alias = add_reactants(initial_name, added_reactants)
    
    return add_products(name_after_reactants_alias, added_products)

def replace_mpculeids_with_names(H_name, mpculeid_dict):
    reaction_as_list = H_name.split()
    
    for index, obj in enumerate(reaction_as_list[:]):
        kinetiscope_name = mpculeid_dict.get(obj, None)
        if kinetiscope_name:
            reaction_as_list[index] = kinetiscope_name
    
    return " ".join(reaction_as_list)
    
def write_kinetiscope_name(H_name, mpculeid_dict, added_reactants, added_products):
    reaction_with_names = replace_mpculeids_with_names(H_name, mpculeid_dict)
    kinetiscope_name = \
        add_reacting_species(reaction_with_names, added_reactants, added_products)
    
    return kinetiscope_name 

def build_rxn_object(HiPRGen_rxn, kinetiscope_name, rate_constant, order, marker_species):
    return Kinetiscope_Reaction(HiPRGen_rxn, kinetiscope_name, rate_constant, order, marker_species)

def replace_tag_with_shorthand(marker_species, shorthand_dict):
    tag = marker_species[-1]
    for longhand, shorthand in shorthand_dict.items():
        if longhand in tag:
            tag = tag.replace(longhand, shorthand)
    marker_species[-1] = tag
    return marker_species