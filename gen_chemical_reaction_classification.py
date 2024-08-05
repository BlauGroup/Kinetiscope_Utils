# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 11:01:50 2024

@author: jacob
"""
from reaction_classification_utilities import find_mpculeid_charge

def find_reactant_charges(reaction):
    reactant_1_charge = find_mpculeid_charge(reaction.reactants[0])
    reactant_2_charge = find_mpculeid_charge(reaction.reactants[1])
    return reactant_1_charge, reactant_2_charge

def reaction_is_neutral(reaction):
    reactant_1_charge, reactant_2_charge = find_reactant_charges(reaction)
    return reactant_1_charge == reactant_2_charge == 0

def reaction_is_ion_molecule(reaction):
    reactant_1_charge, reactant_2_charge = find_reactant_charges(reaction)
    charges_are_different = reactant_1_charge != reactant_2_charge
    one_charge_zero = reactant_1_charge == 0 or reactant_2_charge == 0
    return charges_are_different and one_charge_zero

def determine_bimolecular_subtype(reaction):
    if reaction_is_neutral(reaction):
        return "neutral"
    elif reaction_is_ion_molecule(reaction):
        return "ion-molecule"
    else:
        return "ion-ion"

def write_bimolecular_classification(reaction):
    subtype = determine_bimolecular_subtype(reaction)
    if "combination" in reaction.tag:
        subclass = "combination"
        return ["chemical", "bimolecular", subclass, subtype, reaction.tag]
    
    subclass = "biproduct"
    subcategories = [
        "proton_transfer",
        "hydride_abstraction",
        "H_atom_abstraction",
        "proton_coupled_electron_transfer",
        "electron_transfer",
        "reaction"
    ]
    
    for subcategory in subcategories:
        if subcategory in reaction.tag:
            return ["chemical", "bimolecular", subclass, subtype, subcategory, reaction.tag]

def write_unimolecular_classification(reaction):
    subclass = "fragmentation" if "fragmentation" in reaction.tag else "isomerization"
    return ["chemical", "unimolecular", subclass, reaction.tag]

def write_chemical_classification(reaction):
    number_of_reactants = len(reaction.reactants)
    if number_of_reactants == 1:
        return write_unimolecular_classification(reaction)
    return write_bimolecular_classification(reaction)

def write_ionization_classification(reaction):
    return ["ionization", reaction.tag]

# def add_to_rxn_dict(rxn, rxn_dict, ionization_classifications):
def write_reaction_classification(reaction):
    ionization_classifications = set(["positive_ionization", "electron_attachment",
                                      "electron_cation_recombination"])
    if reaction.tag in ionization_classifications:
        return write_ionization_classification(reaction)
    return write_chemical_classification(reaction)



# def add_reaction_to_nested_dict(rxn_dict, keys, reaction):
#     """Helper function to add a reaction to a nested dictionary."""
#     d = rxn_dict
#     for key in keys[:-1]:
#         if key not in d:
#             d[key] = {}
#         d = d[key]
#     if keys[-1] not in d:
#         d[keys[-1]] = []
#     d[keys[-1]].append(reaction)

# def add_reaction_to_dict(rxn, rxn_dict, category, reaction_type=None, subclass=None, subtype=None, subcategory=None):
#     keys = [category]
#     possible_keys = [reaction_type, subclass, subtype, subcategory]
#     for possible_key in possible_keys:
#         if possible_key is None:
#             break
#         keys.append(possible_key)
#     keys.append(reaction.tag)
#     add_reaction_to_nested_dict(rxn_dict, keys, reaction)
#     return reaction_dict