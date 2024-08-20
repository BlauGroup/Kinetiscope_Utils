# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 14:23:37 2024

@author: jacob
"""

from kinetiscope_reaction_writing_utilities import (
build_rxn_object,
write_kinetiscope_name,
replace_tag_with_shorthand
)

def write_name_with_excited_reactant(kinetiscope_name, reactant):
    reactant = reactant.strip()
    new_name_list = [name + "*" if name == reactant else name for name in kinetiscope_name.split()]
    return " ".join(new_name_list)

def determine_ordinal_number_order(H_rxn):
    if len(H_rxn.reactants) == 1:
        return "1st_order"
    return "2nd_order"

def add_reaction_with_excited_reactant(name_without_excitation, H_rxn, reaction_list, reactant, rate_constant_dict, shorthand_dict):
    ordinal_number_order = determine_ordinal_number_order(H_rxn)
    order = 2 if ordinal_number_order == "2nd_order" else 1 
    rate_constant = rate_constant_dict["chemical"].get(ordinal_number_order, None)
    name_with_excited_reactant = write_name_with_excited_reactant(name_without_excitation, reactant)
    marker_species = replace_tag_with_shorthand(H_rxn.classification_list, shorthand_dict)
    excited_reaction = build_rxn_object(H_rxn, name_with_excited_reactant, rate_constant, order, marker_species)
    reaction_list.append(excited_reaction)
    return reaction_list
    
def write_dexcitation_reaction_name(excitation_reaction_name, reactant):
    reactant = reactant.strip(" ")
    dexcitation_list = [name for name in excitation_reaction_name.split() if reactant in name]
    dexcitation_list.reverse()
    return " => ".join(dexcitation_list)

def write_excitation_reaction_name(reactant):
    excited_reactant = reactant + "* "
    return reactant + " + LEE => " + excited_reactant + "+ TE + excitation"

def build_dexcitation_reaction(H_rxn, dexcitation_reaction_name, rate_constant_dict, reaction_list):
    de_order = 1
    de_rate_constant = 2.0E+06
    de_species = ["dexcitation"]
    dexcitation_reaction_name = dexcitation_reaction_name + " + dexcitation"
    dexcitation_reaction = build_rxn_object(H_rxn, dexcitation_reaction_name, de_rate_constant, de_order, de_species)
    return dexcitation_reaction

def add_dexcitation_reaction(excitation_name, reactant, H_rxn, rate_constant_dict, reaction_list):
    dexcitation_name = write_dexcitation_reaction_name(excitation_name, reactant)
    dexcitation_reaction = build_dexcitation_reaction(H_rxn, dexcitation_name, rate_constant_dict, reaction_list)
    reaction_list.append(dexcitation_reaction)
    return reaction_list

def add_excitation_reaction(H_rxn, rate_constant_dict, excitation_set, excitation_name, reaction_list):
    excitation_set.add(excitation_name)
    ex_order = 2
    ex_rate_constant = rate_constant_dict["chemical"].get("2nd_order", None)
    ex_species = ["excitation"]
    excitation_reaction = build_rxn_object(H_rxn, excitation_name, ex_rate_constant, ex_order, ex_species)
    reaction_list.append(excitation_reaction)
    return reaction_list, excitation_set

def add_excitation_and_dexcitation(reactant, H_rxn, excitation_set, rate_constant_dict, reaction_list, excitation_name):
    reaction_list, excitation_set = add_excitation_reaction(H_rxn, rate_constant_dict, excitation_set, excitation_name, reaction_list)
    reaction_list = add_dexcitation_reaction(excitation_name, reactant, H_rxn, rate_constant_dict, reaction_list)
    
    return reaction_list, excitation_set
    
def add_excitation_dexcitation_if_new(reactant, H_rxn, excitation_set, rate_constant_dict, reaction_list):
    excitation_name = write_excitation_reaction_name(reactant)
    excitation_not_added = excitation_name not in excitation_set
    if excitation_not_added:
       reaction_list, excitation_set = \
           add_excitation_and_dexcitation(reactant, H_rxn, excitation_set, rate_constant_dict, reaction_list, excitation_name)
    return reaction_list, excitation_set
 
def find_reactants_to_excite(name_without_excitation):
    reactants = name_without_excitation.split(" => ")[0]
    reactant_list = reactants.split()
    if "+" in reactant_list:
        reactant_list.remove("+")
    return reactant_list

def write_name_without_excitation(H_rxn, mpculeid_name_dict, marker_species_shorthand):
    marker_species = H_rxn.classification_list[:]
    marker_species = replace_tag_with_shorthand(marker_species, marker_species_shorthand)
    added_reactants = None
    H_name = H_rxn.name
    kinetiscope_name = \
        write_kinetiscope_name(H_name, mpculeid_name_dict, added_reactants, marker_species)
    return kinetiscope_name

def handle_phase1_chemical_reactions(H_rxn, mpculeid_name_dict, rate_constant_dict, marker_species_shorthand, excitation_set):
    reaction_list = []
    name_without_excitation = write_name_without_excitation(H_rxn, mpculeid_name_dict, marker_species_shorthand)
    
    reactants_to_excite = find_reactants_to_excite(name_without_excitation)
    
    for reactant in reactants_to_excite:
        reaction_list, excitation_set = \
            add_excitation_dexcitation_if_new(reactant, H_rxn, excitation_set, rate_constant_dict, reaction_list)
        reaction_list = \
            add_reaction_with_excited_reactant(name_without_excitation, H_rxn, reaction_list, reactant, rate_constant_dict, marker_species_shorthand)
        
    return reaction_list, excitation_set

def handle_phase2_chemical_reactions(H_rxn, mpculeid_name_dict, rate_constant_dict, marker_species_shorthand):
    added_reactants = None
    added_products = replace_tag_with_shorthand(H_rxn.classification_list, marker_species_shorthand)
    H_name = H_rxn.name
    kinetiscope_name = write_kinetiscope_name(H_name, mpculeid_name_dict, added_reactants, added_products)
    ordinal_number_order = determine_ordinal_number_order(H_rxn)
    order = 2 if ordinal_number_order == "2nd_order" else 1 
    rate_constant = rate_constant_dict["chemical"].get(ordinal_number_order, None)
    # print(build_rxn_object(H_rxn, kinetiscope_name, rate_constant, order, added_products))
    return build_rxn_object(H_rxn, kinetiscope_name, rate_constant, order, added_products)

# def handle_phase_1_chemical_rxns(H_rxn, mpculeid_name_dict, rate_constant_dict, marker_species_shorthand, excitation_set):
    
# def handle_chemical_rxns(H_rxn, mpculeid_name_dict, rate_constant_dict, marker_species_shorthand, excitation_set):
#     if H_rxn.phase == 1:
#         return handle_phase_1_chemical_reactions(H_rxn, mpculeid_name_dict, rate_constant_dict, marker_species_shorthand, excitation_set)
#     return handle_phase_2_chemical_reactions(H_rxn, mpculeid_name_dict, rate_constant_dict, marker_species_shorthand)