# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 14:24:02 2024

@author: jacob
"""
from reaction_classification_utilities import find_mpculeid_charge
from kinetiscope_reaction_writing_utilities import (
write_kinetiscope_name,
build_rxn_object)

def reaction_is_positive(H_reaction):
    reactant_mpculeid = H_reaction.reactants[0]
    reactant_charge = find_mpculeid_charge(reactant_mpculeid)
    reactant_is_positive = reactant_charge > 0
    return  reactant_is_positive  

def build_ionization_reactions(H_reaction, mpculeid_dict, rate_constant_key, rate_constant_dict, ionization_type, order, added_reactants, added_products):
    added_products.append(ionization_type)
    kinetiscope_name = write_kinetiscope_name(H_reaction.name, mpculeid_dict, added_reactants, added_products)
    rate_constant =  rate_constant_dict["ionization"][ionization_type].get(rate_constant_key, None)
    
    return build_rxn_object(H_reaction, kinetiscope_name, rate_constant, order, added_products)

def build_light_ionization(H_reaction, mpculeid_dict, rate_constant_dict):
    added_products =  H_reaction.classification_list[:-1] # does not include tag
    absorber = H_reaction.reactants[0]
    return build_ionization_reactions(H_reaction, mpculeid_dict, absorber, rate_constant_dict, "absorption", 1, None, added_products)

def build_electron_ionization_reactions(H_reaction, mpculeid_dict, rate_constant_dict):
    EI_reaction_list = []
    reactant_product_dict = {
        "eV_80": ["eV_55", "LEE"],
        "eV_55": ["eV_30", "LEE"],
        "eV_30": ["2LEE"]
                            }
    for reactant, product_list in reactant_product_dict.items():
        product_list = H_reaction.classification_list[:-1] + product_list
        EI_reaction = build_ionization_reactions(H_reaction, mpculeid_dict, reactant, rate_constant_dict, "electron_ionization", 2, [reactant], product_list)
        EI_reaction_list.append(EI_reaction)
        
    return EI_reaction_list

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

def set_ionization_type_and_rate_constant_key(H_reaction):
    if "recombination" in H_reaction.tag:
        
        ionization_type = "recombination"
        rate_constant_key = "all"
        
    else:
        
        ionization_type = "attachment"
        
        if reaction_is_positive(H_reaction):
            
            rate_constant_key = "positive"
            
        else:
            
            rate_constant_key = "neutral"
    
    return ionization_type, rate_constant_key

def build_attachment_or_recombination(H_reaction, mpculeid_name_dict, rate_constant_dict):
    ionization_type, rate_constant_key = set_ionization_type_and_rate_constant_key(H_reaction) 
    product_list = H_reaction.classification_list[:-1]
    reactant_list = ["TE"]   
    reaction = build_ionization_reactions(H_reaction, mpculeid_name_dict, rate_constant_key, rate_constant_dict, ionization_type, 2, reactant_list, product_list)
    return [reaction]

def handle_ionization_reactions(H_reaction, mpculeid_dict, rate_constant_dict):
    if H_reaction.tag == "positive_ionization":
        return build_positive_ionization_reactions(H_reaction, mpculeid_dict, rate_constant_dict)
    return build_attachment_or_recombination(H_reaction, mpculeid_dict, rate_constant_dict)