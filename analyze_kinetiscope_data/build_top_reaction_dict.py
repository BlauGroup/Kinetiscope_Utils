# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 13:29:27 2024

@author: jacob
"""

from associate_indicies_with_frequencies import build_index_freq_dict
from create_kinetiscope_simulation_reactions import build_index_reaction_dict

def should_update_reaction(top_formation_reaction, product, reaction):
    """
    Determines if a reaction should be updated in the top reactions dictionary.

    Parameters
    ----------
    top_reactions : dict
        A dictionary of top reactions where keys are products and values are reaction objects.
    product : str
        The product whose reaction is being considered for update.
    reaction : object
        The reaction object to be checked for update.

    Returns
    -------
    bool
        True if the reaction should be updated, False otherwise.
    """
    
    product_is_new = product not in top_formation_reaction
    
    if product_is_new:
        
        return True

    current_reaction = top_formation_reaction[product]
    
    return reaction.selection_freq > current_reaction.selection_freq

def build_top_reaction_dict(index_reaction_dict):
    """
    Builds a dictionary of top reactions for each product based on selection 
    frequency.

    Parameters
    ----------
    index_reaction_dict : dict
        A dictionary where keys are reaction indices and values are reaction 
        objects containing product information and selection frequencies.

    Returns
    -------
    dict
        A dictionary where keys are products and values are the reaction object 
        with the highest selection frequency associated with that product.
    """
    
    top_formation_reactions = {}
    
    for reaction in index_reaction_dict.values():
        products = reaction.products
        for product in products:
            if should_update_reaction(top_formation_reactions, product, reaction):
                top_formation_reactions[product] = reaction

    return top_formation_reactions

def find_top_formation_reactions(select_freq_file, reaction_name_file, start_index, end_index):
    """
    Finds the top formation reactions for each product based on selection
    frequencies.
    
    Parameters
    ----------
    select_freq_file : str
        The path to the file containing selection frequencies for reactions. 
        This file is used to build a dictionary mapping reaction indices to 
        their selection frequencies.
    reaction_name_file : str
        The path to the file containing reaction names and details. This file 
        is used to build a dictionary mapping reaction indices to reaction 
        objects.
    start_index : int
        The starting index for the range of lines to be considered in the 
        selection frequencies file.
    end_index : int
        The ending index for the range of lines to be considered in the 
        selection frequencies file.
    
    Returns
    -------
    dict
        A dictionary where keys are products and values are the reaction object 
        with the highest selection frequency associated with each product. This 
        dictionary represents the top formation reactions for each product.
    """
    
    index_freq_dict = (
        build_index_freq_dict(select_freq_file, start_index, end_index)
    )
    
    index_reaction_dict = (
        build_index_reaction_dict(reaction_name_file, index_freq_dict)
    )
    
    top_formation_reactions = build_top_reaction_dict(index_reaction_dict)
    
    return top_formation_reactions