# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 12:29:38 2024

@author: jacob
"""

from kinetiscope_reaction_writing_utilities import (
    determine_ordinal_number_order,
    replace_tag_with_shorthand,
    build_rxn_object,
    write_kinetiscope_name)

def write_phase_2_chemical_reactions(HiPRGen_rxn, reaction_writing_data):
    added_reactants = None #no reactant tags added to these rxns
    marker_species_shorthand = reaction_writing_data.marker_species_dict
    mpculeid_name_dict = reaction_writing_data.mpculeid_dict
    rate_constant_dict = reaction_writing_data.rate_constant_dict
    
    added_products = (
        replace_tag_with_shorthand(HiPRGen_rxn.classification_list, marker_species_shorthand)
    )
    
    HiPRGen_name = HiPRGen_rxn.name
    
    kinetiscope_name = (
        write_kinetiscope_name(HiPRGen_name, mpculeid_name_dict, added_reactants, added_products)
    )
    
    ordinal_number_order = determine_ordinal_number_order(HiPRGen_rxn)
    order = 2 if ordinal_number_order == "2nd_order" else 1 
    rate_constant = (
        rate_constant_dict["chemical"].get(ordinal_number_order, None)
    )
    
    #we return the reaction as the single element in a list just because the 
    #function call expects a list
    
    rxn_list = (
        [build_rxn_object(HiPRGen_rxn, kinetiscope_name, rate_constant, order, added_products)]
    )
    
    #we also return reaction_writing_data because some other functions alter it

    return rxn_list, reaction_writing_data