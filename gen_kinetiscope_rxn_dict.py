# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:54:29 2024

@author: JRMilton
"""
__version__ = '1.1.0'

import os
import sys
from monty.serialization import loadfn, dumpfn
# from Rxn_classes import HiPRGen_Reaction, Kinetiscope_Reaction
# import sys
from write_kinetiscope_ionization_reactions import handle_ionization_reactions
from write_kinetiscope_chemical_reactions import (
handle_phase1_chemical_reactions,
handle_phase2_chemical_reactions
)
import csv

def collect_lists_from_nested_dict(d):
    collected_items = []
    
    def recurse_through_dict(d):
        for value in d.values():
            if isinstance(value, dict):
                # If the value is a dictionary, recurse into it
                recurse_through_dict(value)
            elif isinstance(value, list):
                # If the value is a list, extend the collected_items list with the items
                collected_items.extend(value)

    recurse_through_dict(d)
    return collected_items

def find_kinetiscope_reactants_and_products(k_rxn):
    k_rxn_name = k_rxn.kinetiscope_name
    reactant_str, product_str = k_rxn_name.split("=>")
    reactant_list = reactant_str.split()
    product_list = product_str.split()[:3]
    if "+" in reactant_list:
        reactant_list.remove("+")
    if "+" in product_list:
        product_list.remove("+")
    return reactant_list, product_list

def are_species_duplicates(species_list):
    return species_list[0] == species_list[1]

def reaction_has_duplicate_species(k_rxn):
    reactant_list, product_list = find_kinetiscope_reactants_and_products(k_rxn)
    reactants_are_duplicates = len(reactant_list) == 2 and are_species_duplicates(reactant_list)
    products_are_duplicates = len(product_list) == 2 and are_species_duplicates(product_list)
    # if not products_are_duplicates and len(product_list) == 2:
    #     print(product_list)
    return reactants_are_duplicates or products_are_duplicates

def modify_kinetiscope_name(k_rxn, species_list, side_identifier):
    index = 0 if side_identifier == "reactants" else 1
    k_rxn_name = k_rxn.kinetiscope_name
    old_rxn_name = k_rxn_name.split("=>")
    old_rxn_name[index] = " + ".join(species_list)
    new_rxn_name = " =>".join(old_rxn_name)
    return new_rxn_name

def star_test(species_list):
    for species in species_list:
        if "*" in species:
            return True
    return False

def replace_second_star_reactant(reactant_list):
    reactant_list[1] = reactant_list[1].replace("*", "")
    return reactant_list

def replace_species_with_2_species(species_list):
    species = species_list[0]
    two_species = "2" + species
    return two_species

def correct_kinetiscope_name(k_rxn):
    reactant_list, product_list = find_kinetiscope_reactants_and_products(k_rxn)
    reactants_are_duplicates = len(reactant_list) == 2 and are_species_duplicates(reactant_list)
    if reactants_are_duplicates:
        if star_test(reactant_list):
            reactant_list = replace_second_star_reactant(reactant_list)
        else:
            reactant_list = [replace_species_with_2_species(reactant_list)]
        k_rxn.kinetiscope_name = modify_kinetiscope_name(k_rxn, reactant_list, "reactants")
    products_are_duplicates = len(product_list) == 2 and are_species_duplicates(product_list)
    if products_are_duplicates: #should be true at least once for full
        reactant_list = [replace_species_with_2_species(product_list)]
        k_rxn.kinetiscope_name = modify_kinetiscope_name(k_rxn, product_list, "products")
    return k_rxn

def order_kinetiscope_reactions(kinetiscope_reactions, supercategory_order, supercategories_with_subcategories, subcategory_order):
    def get_supercategory_index(reaction):
        for i, supercategory in enumerate(supercategory_order):
            if supercategory in reaction.marker_species:
                return i
        # Return a large number if no supercategory is found
        return len(supercategory_order)

    def get_subcategory_index(reaction):
        # Check if the reaction's supercategory has subcategories
        for supercategory in supercategory_order:
            if supercategory in reaction.marker_species and supercategory in supercategories_with_subcategories:
                for i, subcategory in enumerate(subcategory_order):
                    if subcategory in reaction.marker_species:
                        return i
        # Return a large number if no subcategory is found
        return len(subcategory_order)

    # Sort first by supercategory, then by subcategory
    ordered_reactions = sorted(
        kinetiscope_reactions,
        key=lambda reaction: (get_supercategory_index(reaction), get_subcategory_index(reaction))
    )
    
    return ordered_reactions

def shorten_PCET(reaction):
    reaction.kinetiscope_name = reaction.kinetiscope_name.replace("proton_coupled_electron_transfer", "PCET")
    return reaction

kinetiscope_reaction_list = []
os.chdir("G:/My Drive/Kinetiscope/new_kinetiscope_naming_080224")
# test_rxns = "HiPRGen_rxns_to_name.json"
# HiPRGen_reaction_list = loadfn(test_rxns)
full_rxns = "HiPRGen_rxns_to_name_full.json"
HiPRGen_reaction_list = loadfn(full_rxns)
# name_mpculeid_file = "name_test_mpculeid_080624.json"
name_mpculeid_file = "name_full_mpculeid_080624.json"
name_mpculeid_dict = loadfn(name_mpculeid_file)
mpculeid_name_dict = {mpculeid: name for name, mpculeid in name_mpculeid_dict.items()}

rate_constant_dict = { #values calculated in excel and taken as scientific notation values
    "ionization":
        {"absorption":{
            "4864aee73a83d357c31fadd50b81e3cd-C10H20O2-0-1":1.4e-01,
            "00a7dcc352b0d613f58e850935bf5609-C10H14O1-0-1":1.0e-01,
            "bfed458e642b8daa1eab6bc02d5e5682-C18H15S1-1-1":1.8e-01,
            "9a8a88b8b92c714d7f65b8526ffabc7a-C4F9O3S1-m1-1":5.8e-01,
            "17f31f89123edbaa0e3b9c7eb49d26f3-C8H4N1O2-m1-1":1.5E-01},
         "electron_ionization":{
             "eV_80":1.4e+15,
             "eV_55":1.2e+15,
             "eV_30":8.5e+14,
             "LEE":1.1E+14},
         "recombination":{
             "all": 4.1e+13},
         "attachment":{
             "positive":4.1e+13,
             "neutral":3.1e+12}
         },
    "chemical":
        {"1st_order":6.2e+12,
         "2nd_order":1.3e+12}
        }
    
marker_species_shorthand = {
    "radical_cation":"rc",
    "radical_anion":"ra",
    "neutral_radical":"nr",
    "cation":"c",
    "anion":"a",
    # "radical":"r",
    "neutral":"n",
    "proton_coupled_electron_transfer":"PCET"}



excitation_set = set()
# LEE_H_rxn = None
# LEE_collision_name = "LEE => TE"
# LEE_rate_constant = 5.1E+14
# LEE_order = 1
# LEE_marker_species = None

# LEE_rxn = \
#     build_rxn_object(LEE_H_rxn, LEE_collision_name, LEE_rate_constant, LEE_order, LEE_marker_species)

for rxn_list in HiPRGen_reaction_list["ionization"].values():
    for HiPRGen_rxn in rxn_list:
        ionization_reaction_list = \
            handle_ionization_reactions(HiPRGen_rxn, mpculeid_name_dict, rate_constant_dict)
        kinetiscope_reaction_list.extend(ionization_reaction_list)

chemical_reaction_list = collect_lists_from_nested_dict(HiPRGen_reaction_list["chemical"])

for H_rxn in chemical_reaction_list:
    if H_rxn.phase == 1:
        chemical_reaction_list, excitation_set = \
            handle_phase1_chemical_reactions(H_rxn, mpculeid_name_dict, rate_constant_dict, marker_species_shorthand, excitation_set)
        kinetiscope_reaction_list.extend(chemical_reaction_list)
    else:
        reaction = handle_phase2_chemical_reactions(H_rxn, mpculeid_name_dict, rate_constant_dict, marker_species_shorthand)
        kinetiscope_reaction_list.append(reaction)


for index, reaction in enumerate(kinetiscope_reaction_list):
    if reaction_has_duplicate_species(reaction):
        kinetiscope_reaction_list[index] = correct_kinetiscope_name(reaction)
        
    if "proton_coupled_electron_transfer" in reaction.marker_species:
        kinetiscope_reaction_list[index] = shorten_PCET(reaction)

supercategory_order = [
    "absorption", "electron_ionization", "recombination", "attachment", "excitation",
    "dexcitation", "fragmentation", "isomerization", "ion-ion", "ion-molecule", "neutral"]

supercategories_with_subcategories = set(["ion-ion", "ion-molecule", "neutral"])

subcategory_order = [
    "proton_transfer", "H_atom_abstraction", "hydride_abstraction", 
    "proton_coupled_electron_transfer", "electron_transfer", "reaction"
    ]
    
ordered_reactions = order_kinetiscope_reactions(
    kinetiscope_reaction_list, 
    supercategory_order, 
    supercategories_with_subcategories, 
    subcategory_order
)

# dumpfn(ordered_reactions, "ordered_kinetiscope_reactions.json")

print('Writing reactions to csv file...')
# with open('Kinetiscope_rxn_template.csv', newline = "") as csvfile:
    # reader = csv.reader(csvfile)
    # fields = list(next(reader))

dict_list = []

for reaction in ordered_reactions:
    rate_coefficient_format = 3 if "absorption" in reaction.marker_species else 0
    csv_dict = {}
    csv_dict['# equation'] = reaction.kinetiscope_name
    csv_dict['fwd_A'] = 1
    csv_dict['fwd_temp_coeff'] = 0
    csv_dict['fwd_Ea'] = 0
    csv_dict['fwd_k'] = reaction.rate_coefficient if "absorption" not in reaction.marker_species else 1
    csv_dict['rev_A'] = 1
    csv_dict['rev_temp_coeff'] = 0
    csv_dict['rev_Ea'] = 0
    csv_dict['rev_k'] = 1
    csv_dict['fwd_k0'] = 1
    csv_dict['rev_k0'] = 1
    csv_dict['alpha_alv'] = 0.5
    csv_dict['equil_potential'] = 0.5
    csv_dict['num_electrons'] = 0
    csv_dict['fwd_prog_k'] = reaction.rate_coefficient if "absorption" in reaction.marker_species else 1
    csv_dict['rev_prog_k'] = 1
    csv_dict['non_stoichiometric'] = 0
    csv_dict['rate_constant_format'] = rate_coefficient_format
    dict_list.append(csv_dict)

# os.chdir(new_dir)
    
with open("euvl_full_reactions.csv", 'w', newline = "") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames = dict_list[0].keys())
    writer.writeheader()
    for reaction in dict_list:
        writer.writerow(reaction)
  
print('Done!') 
    
#         A dictionary with kinetiscope_names as keys and their counts as values.
#     """
#     # Dictionary to store the count of each kinetiscope_name
#     name_counts = {}

#     # Iterate over the list of reactions
#     for reaction in kinetiscope_reaction_list:
#         # Retrieve kinetiscope_name
#         kinetiscope_name = reaction.kinetiscope_name
        
#         # Update the count in the dictionary
#         if kinetiscope_name in name_counts:
#             name_counts[kinetiscope_name] += 1
#         else:
#             name_counts[kinetiscope_name] = 1

#     for name, count in name_counts.items():
#         if count > 1:
#             print(f"{name}: {count}")

# counts = count_kinetiscope_names(kinetiscope_reaction_list)
# dumpfn(kinetiscope_reaction_list, "kinetiscope_full_reaction_list.json")
        
