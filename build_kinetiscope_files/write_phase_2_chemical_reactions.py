# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 12:29:38 2024

@author: jacob
"""

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