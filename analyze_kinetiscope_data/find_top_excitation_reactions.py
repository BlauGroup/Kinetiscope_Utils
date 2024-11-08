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

def find_all_reactions_forming_product(product, index_reaction_dict, index_freq_dict):
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

    all_reactions_forming_product = []

    for index, reaction in index_reaction_dict.items():
        products = reaction.products
        if product in products:
            product_removed_punctuation = product.replace("-","").replace("_","")
            product_is_real = not product_removed_punctuation.isalpha()
            if product_is_real:
                selection_freq = index_freq_dict.get(index, None)
                if selection_freq > 0:
                    all_reactions_forming_product.append({reaction: selection_freq})
    
    # Use the helper function to sort by selection frequency
    all_reactions_forming_product = sort_reactions_by_selection_freq(all_reactions_forming_product)
    
    previous_selection_freq = 0

    for reaction_dict in all_reactions_forming_product:
        for selection_freq in reaction_dict.values():
            if previous_selection_freq != 0:
                try:
                    assert previous_selection_freq >= selection_freq
                except:
                    print(previous_selection_freq, selection_freq)
                    sys.exit()
            previous_selection_freq = selection_freq

    return all_reactions_forming_product
       
def build_top_reaction_dict(index_reaction_dict, index_freq_dict):
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
            if product not in top_formation_reactions:
                top_formation_reactions[product] = find_all_reactions_forming_product(product, index_reaction_dict, index_freq_dict)

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
    
    top_formation_reactions = build_top_reaction_dict(index_reaction_dict, index_freq_dict)
    
    return top_formation_reactions


if __name__ == "__main__":
    kinetiscope_files_dir = (
        r"G:\My Drive\Kinetiscope\production_simulations_092124"
    )

    correct_path_change_dir(kinetiscope_files_dir)

    select_freq_file = "excitation_selection_freq_092524.txt"
    reaction_name_file = "excitation_reactions_092524.txt"

    top_reaction_dict = find_top_formation_reactions(
        select_freq_file,
        reaction_name_file,
        start_index=7,
        end_index=5384
    )
    max_length = 0  # Initialize to zero
    
    starting_species = [
        "PtBMAb_COO_tbut_0",
        "Nf_-1",
        "TPS_H1_+1_#1",
        "COO_CN_C6H4_-1",
        "PHSb_phol_0"]

    # for product, reaction_list in top_reaction_dict.items():
    #     if product not in starting_species:
    #         max_length = max(max_length, len(reaction_list))
    #         if len(reaction_list) == 111:
    #             print(product)

    # print("The maximum length of lists in top_reaction_dict is:", max_length)
    # for reaction_list in top_reaction_dict.values():
    #     for dictionary in reaction_list:
    #         reaction = list(dictionary.keys())[0]
    #         if "*" in reaction.kinetiscope_name:
    #             if list(dictionary.values())[0] > 100:
    #                 print(reaction.kinetiscope_name)
    #                 print(list(dictionary.values())[0])
    # Initialize a list to store (frequency, reaction, dictionary) tuples
    filtered_reactions = []
    unique_reactions = set()
    
    # Iterate through the top_reaction_dict
    for reaction_list in top_reaction_dict.values():
        for dictionary in reaction_list:
            reaction = list(dictionary.keys())[0]
            frequency = list(dictionary.values())[0]
            
            # Check for "*" in the reaction name and frequency > 100
            if "*" in reaction.kinetiscope_name and frequency > 100:
                if reaction.kinetiscope_name not in unique_reactions:
                    # Append tuple of (frequency, reaction, dictionary) to the list
                    unique_reactions.add(reaction.kinetiscope_name)
                    filtered_reactions.append((frequency, reaction, dictionary))
    
    # Sort the list by frequency in descending order
    filtered_reactions.sort(reverse=True, key=lambda x: x[0])
    
    # Print each reaction and dictionary in sorted order
    for frequency, reaction, dictionary in filtered_reactions:
        print("Reaction:", reaction.kinetiscope_name)
        print("Frequency:", frequency)
        print("Dictionary:", dictionary)
        print()  # For better readability between entries
