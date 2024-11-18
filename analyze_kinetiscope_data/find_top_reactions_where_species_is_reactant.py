# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 13:29:27 2024

@author: jacob
"""

from associate_indicies_with_frequencies import build_index_freq_dict
from create_kinetiscope_simulation_reactions import build_index_reaction_dict
import sys
sys.path.append('../common')
from utilities import correct_path_change_dir

"""
This script defines a function, find_top_formation_reactions, that creates a 
dictionary linking each species' chemical name in the simulations to the 
reaction where it is a product with the highest selection frequency.
"""

# def should_update_reaction(top_formation_reaction, product, reaction):
#     """
#     Determines if a reaction should be updated in the top reactions dictionary.

#     Parameters
#     ----------
#     top_reactions : dict
#         A dictionary of top reactions where keys are products and values are reaction objects.
#     product : str
#         The product whose reaction is being considered for update.
#     reaction : object
#         The reaction object to be checked for update.

#     Returns
#     -------
#     bool
#         True if the reaction should be updated, False otherwise.
#     """
    
#     product_is_new = product not in top_formation_reaction
    
#     if product_is_new:
        
#         return True

#     current_reaction = top_formation_reaction[product]
    
#     return reaction.selection_freq > current_reaction.selection_freq

def find_all_reactions_forming_reactant(reactant, index_reaction_dict, index_freq_dict):
    def sort_reactions_by_selection_freq(reactions_list):
        """
        Sorts a list of dictionaries of the form {reaction: selection_frequency}
        in descending order by selection frequency.
    
        Parameters
        ----------
        reactions_list : list
            A list of dictionaries where each dictionary has one key-value pair
            with the reaction as the key and its selection frequency as the value.
    
        Returns
        -------
        list
            The sorted list of dictionaries by selection frequency in descending order.
        """
        return sorted(reactions_list, key=lambda x: list(x.values())[0], reverse=True)

    all_species_is_reactant = []

    for index, reaction in index_reaction_dict.items():
        reactants = reaction.reactants
        if reactant in reactants:
            reactant_removed_punctuation =reactant.replace("-","").replace("_","")
            reactant_is_real = not reactant_removed_punctuation.isalpha()
            if reactant_is_real:
                selection_freq = index_freq_dict.get(index, None)
                if selection_freq > 0:
                    all_species_is_reactant.append({reaction: selection_freq})
    
    # Use the helper function to sort by selection frequency
    all_reactions_species_is_reactant = sort_reactions_by_selection_freq(all_species_is_reactant)
    
    previous_selection_freq = 0

    for reaction_dict in all_reactions_species_is_reactant:
        for selection_freq in reaction_dict.values():
            if previous_selection_freq != 0:
                try:
                    assert previous_selection_freq >= selection_freq
                except AssertionError:
                    print(previous_selection_freq, selection_freq)
                    sys.exit()
            previous_selection_freq = selection_freq

    return all_reactions_species_is_reactant
       
# def build_top_reaction_dict(index_reaction_dict, index_freq_dict):
#     """
#     Builds a dictionary of top reactions for each product based on selection 
#     frequency.

#     Parameters
#     ----------
#     index_reaction_dict : dict
#         A dictionary where keys are reaction indices and values are reaction 
#         objects containing product information and selection frequencies.

#     Returns
#     -------
#     dict
#         A dictionary where keys are products and values are the reaction object 
#         with the highest selection frequency associated with that product.
#     """
    
#     top_formation_reactions = {}
    
#     for reaction in index_reaction_dict.values():
#         products = reaction.products
#         for product in products:
#             if product not in top_formation_reactions:
#                 top_formation_reactions[product] = find_all_reactions_forming_product(product, index_reaction_dict, index_freq_dict)

#     return top_formation_reactions

# def find_top_formation_reactions(select_freq_file, reaction_name_file, start_index, end_index):
#     """
#     Finds the top formation reactions for each product based on selection
#     frequencies.
    
#     Parameters
#     ----------
#     select_freq_file : str
#         The path to the file containing selection frequencies for reactions. 
#         This file is used to build a dictionary mapping reaction indices to 
#         their selection frequencies.
#     reaction_name_file : str
#         The path to the file containing reaction names and details. This file 
#         is used to build a dictionary mapping reaction indices to reaction 
#         objects.
#     start_index : int
#         The starting index for the range of lines to be considered in the 
#         selection frequencies file.
#     end_index : int
#         The ending index for the range of lines to be considered in the 
#         selection frequencies file.
    
#     Returns
#     -------
#     dict
#         A dictionary where keys are products and values are the reaction object 
#         with the highest selection frequency associated with each product. This 
#         dictionary represents the top formation reactions for each product.
#     """
    
#     index_freq_dict = (
#         build_index_freq_dict(select_freq_file, start_index, end_index)
#     )
    
#     index_reaction_dict = (
#         build_index_reaction_dict(reaction_name_file, index_freq_dict)
#     )
    
#     top_formation_reactions = build_top_reaction_dict(index_reaction_dict, index_freq_dict)
    
#     return top_formation_reactions

def build_top_reaction_dict_for_reactants(index_reaction_dict, index_freq_dict):
    """
    Builds a dictionary of top reactions for each reactant based on selection 
    frequency.

    Parameters
    ----------
    index_reaction_dict : dict
        A dictionary where keys are reaction indices and values are reaction 
        objects containing reactant information and selection frequencies.

    Returns
    -------

    dict
        A dictionary where keys are reactants and values are the reaction object 
        with the highest selection frequency associated with that reactant.
    """
    
    top_reactant_reactions = {}
    
    for reaction in index_reaction_dict.values():
        reactants = reaction.reactants
        for reactant in reactants:
            if reactant not in top_reactant_reactions:
                top_reactant_reactions[reactant] = find_all_reactions_forming_reactant(reactant, index_reaction_dict, index_freq_dict)

    return top_reactant_reactions

def find_top_reactant_reactions(select_freq_file, reaction_name_file, start_index, end_index):
    """
    Finds the top reactions for each reactant based on selection frequencies.
    
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
        A dictionary where keys are reactants and values are the reaction object 
        with the highest selection frequency associated with each reactant. This 
        dictionary represents the top reactions for each reactant.
    """
    
    index_freq_dict = ( 
        build_index_freq_dict(select_freq_file, start_index, end_index)
    )
    
    index_reaction_dict = (
        build_index_reaction_dict(reaction_name_file, index_freq_dict)
    )
    
    top_reactant_reactions = build_top_reaction_dict_for_reactants(index_reaction_dict, index_freq_dict)
    
    return top_reactant_reactions

if __name__ == "__main__":
    kinetiscope_files_dir = r"G:\My Drive\Kinetiscope\production_simulations_092124"
    correct_path_change_dir(kinetiscope_files_dir)
    select_freq_file = "excitation_selection_freq_092524.txt"
    reaction_name_file= "excitation_reactions_092524.txt"
    top_reaction_dict = find_top_reactant_reactions(
        select_freq_file,
        reaction_name_file,
        start_index=7,
        end_index=5384
    )
    for reaction_dict in top_reaction_dict["PHSb_phol_H1_+1_#2"]:
        reaction = list(reaction_dict.keys())[0]
        selection_frequency = list(reaction_dict.values())[0]
        print(f"reaction: {reaction}")
        print(f"selection frequency: {selection_frequency}")
        print()
