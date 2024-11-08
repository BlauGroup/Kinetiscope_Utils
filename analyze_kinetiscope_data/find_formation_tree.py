# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 09:51:24 2024

@author: jacob
"""

from build_top_reaction_dict import find_top_formation_reactions
import os
import sys
sys.path.append('../common')
from utilities import correct_path_change_dir

"""
This script uses the Pathfinding class's method from HiPRGen to find reaction 
pathways leading to the formation of a specific species in our simulations. It 
follows the highest-frequency reactions that form the species and recursively 
searches for the reactants of each reaction. It turns out, if we created a
graph representing our kinetiscope network, with species as nodes connected
by reaction edges, the graph contains cycles; thus, we track whether we have
visited a product before, and if we have, we handle the cycle by choosing the
next most frequent reaction.
"""


def find_reactant_pathways(
     highest_frequency_reaction, highest_select_dict, starting_species, visited,
     reactant_pathways_cache=None
):
    """
    Finds and returns the pathways for each reactant of a given reaction by 
    recursively exploring the pathways of the reactants.

    This function uses recursion to find pathways for each reactant. For each 
    reactant, it calls `find_most_selected_pathway` to get the pathway from the
    reactant to the initial starting species. The recursion continues until it 
    reaches a starting species or a pathway with no further reactions.

    Parameters
    ----------
    highest_frequency_reaction : object
        The reaction object containing the reactants and products.
    highest_select_dict : dict
        A dictionary of reactions keyed by product names.
    starting_species : set
        A set containing names of the starting species.
    reactant_pathways_cache : dict, optional
        A dictionary for caching reactant pathways to avoid recalculation.
    
    Returns
    -------
    list
        A list of reactant pathways for the reactants.
    """
    if reactant_pathways_cache is None:
        reactant_pathways_cache = {}

    reactant_pathways = []

    # Use reaction reactants and products as a cache key (converting to tuples)
    # reaction_key = (tuple(highest_frequency_reaction.reactants),
    #                 tuple(highest_frequency_reaction.products))

    # if reaction_key in reactant_pathways_cache:
    #     return reactant_pathways_cache[reaction_key]

    for reactant in highest_frequency_reaction.reactants:
        # Skip non-real species (e.g., electrons) and starting species
        is_real_species = "eV" not in reactant and not reactant.isalpha()

        if reactant not in starting_species and is_real_species:
            pathway = find_most_selected_pathway(
                reactant, highest_select_dict, starting_species, 
                reactant_pathways_cache, visited
            )
            reactant_pathways.append(pathway)     

    # Cache the computed pathways before returning them
    # reactant_pathways_cache[reaction_key] = reactant_pathways
    return reactant_pathways


def find_most_selected_pathway(
    product, highest_select_dict, starting_species,
    most_selected_pathways_cache=None, visited=None
):
    """
    Finds the most selected pathway for a given product based on the highest
    selection frequencies.

    This function uses recursion to build the pathway from the given product to
    the starting species. It checks if the final reaction in the pathway
    matches the reactants of the highest frequency reaction. If so, it returns
    the combined pathway. If not, it recursively combines pathways from the
    reactants with the current pathway.

    Parameters
    ----------
    product : str
        The name of the product to find the pathway for.
    highest_select_dict : dict
        A dictionary where keys are products and values are the highest
        frequency reactions.
    starting_species : set
        A set containing names of the starting species.
    most_selected_pathways_cache : dict, optional
        A dictionary for caching pathways to avoid redundant calculations.

    Returns
    -------
    list
        The most selected pathway for the given product.
    """

    def final_reaction_forms_reactants(highest_frequency_reaction, pathway):
        """
        Checks if the final reaction in a pathway forms the reactants of the
        highest frequency reaction.

        This function examines the last reaction in the given pathway and
        compares its products with the reactants of the highest frequency
        reaction. It helps to verify if the pathway correctly leads to the
        highest frequency reaction.

        Parameters
        ----------
        highest_frequency_reaction : object
            The reaction object containing the reactants to be checked.
        pathway : list
            A list of reactions representing the pathway.

        Returns
        -------
        bool
            True if the final reaction in the pathway forms the reactants of
            the highest frequency reaction, False otherwise.
        """
        if not pathway:
            return False

        final_reaction = pathway[-1]
        return highest_frequency_reaction.reactants == final_reaction.products

    def combine_reactant_pathways_with_current_pathway(
        reactant_pathways, current_pathway
    ):
        """
        Combines all reactant pathways with the current pathway.

        This function appends the current pathway to each reactant pathway. The
        result is a combined pathway that includes the pathways of the
        reactants leading up to the current pathway.

        Parameters
        ----------
        reactant_pathways : list
            A list of pathways for the reactants.
        current_pathway : list
            The current pathway to combine with reactant pathways.

        Returns
        -------
        list
            The combined pathway with all reactant pathways.
        """

        for pathway in reactant_pathways:
            current_pathway = pathway + current_pathway

        return current_pathway
    
    # problematic_species = {"PHSb_CO_C5H5_0_#2", "PtBMAb_PHSb_phol_COO_0_#1", "PtBMAb_PHSb_COO2_C5H4_0_#2"}

    if visited is None:
        visited = set()

    # if product in problematic_species:
    #     print(visited)
    
    if product in visited:
        return []

    visited.add(product)
    
    index = 0

    if most_selected_pathways_cache is None:
        most_selected_pathways_cache = {}

    if product in most_selected_pathways_cache:
        return most_selected_pathways_cache[product]

    if product not in highest_select_dict:
        # TODO make this not just quit--ask if the name is wrong
        raise KeyError("Product does not have a reaction with a nonzero "
                       "selection frequency.")

    current_pathway = []
    
    all_reactions_forming_product = list(highest_select_dict[product]) #list of dicts
    highest_frequency_reaction = list(all_reactions_forming_product[index].keys())[0]
    
    if product == "PtBMAb_PHSb_phol_COO_0_#1":
        highest_frequency_reaction = list(all_reactions_forming_product[index+1].keys())[0]
    # print(highest_frequency_reaction)
    # # if product in problematic_species:
    # #     print(highest_frequency_reaction)

    # while True:
    #     altered = False  # Track if we change highest_frequency_reaction
    #     for reactant in highest_frequency_reaction.reactants:
    #         if reactant in visited:
    #             print(all_reactions_forming_product)
    #             index += 1
    #             print(highest_frequency_reaction)
    #             highest_frequency_reaction = list(all_reactions_forming_product[index].keys())[0]
    #             print(highest_frequency_reaction)
    #             altered = True  # Set to True since we altered highest_frequency_reaction
    #             break  # Exit the for loop to restart the check
    
    #     if not altered:  # If no change was made, exit the while loop
    #         break
    current_pathway.append(highest_frequency_reaction)

    # returns pathways leading to the formation of the reactants of the
    # highest selection frequency reaction by calling this function recursively

    reactant_pathways = find_reactant_pathways(
        highest_frequency_reaction, highest_select_dict, starting_species,visited,
        most_selected_pathways_cache
    )

    for pathway in reactant_pathways:
        if final_reaction_forms_reactants(highest_frequency_reaction, pathway):

            full_pathway = pathway + current_pathway
            most_selected_pathways_cache[product] = full_pathway

            return full_pathway

    final_pathway = combine_reactant_pathways_with_current_pathway(
        reactant_pathways, current_pathway
    )

    most_selected_pathways_cache[product] = final_pathway

    return final_pathway


def save_pathway_to_file(filename, reactions):
    """
    Saves the pathway information to a text file.
    
    This function creates a text file named after the given chemical name and writes 
    the reaction information to it. Each line in the file contains the reaction name 
    (from `kinetiscope_name`) and its selection frequency.
    
    Parameters
    ----------
    chemical_name : str
        The name of the chemical, used to create the filename for the output file.
    reactions : list
        A list of reaction objects to be saved. Each reaction object must have attributes 
        `kinetiscope_name` (the name of the reaction) and `selection_freq` (the selection 
        frequency of the reaction).
    
    Returns
    -------
    None
        This function does not return a value. It writes the reaction information to a file 
        and prints a confirmation message indicating the file has been saved.
    """

    with open(filename, 'w') as file:
        
        for rxn in reactions:
        
            file.write(f"{rxn.kinetiscope_name}: {rxn.selection_freq}\n")
        
    print(f"Reaction information saved to {filename}")

# name of problematic species is PtBMAb_PHSb_phol_COO_0_#1


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

starting_species = [
    "PtBMAb_COO_tbut_0",
    "Nf_-1",
    "TPS_H1_+1_#1",
    "COO_CN_C6H4_-1",
    "PHSb_phol_0"]

starting_species = set(starting_species)

while True:

    chemical_name = input("Enter a chemical name (or 'end' to stop): ")

    if chemical_name.lower() == 'end':

        print("Exiting the script.")
        break

    filename = f"{chemical_name}.txt"

    if os.path.exists(filename):
        raise OSError(f"The file '{filename}' already exists.")

    pathway = (
        find_most_selected_pathway(chemical_name,
                                   top_reaction_dict,
                                   starting_species)
    )

    save_pathway_to_file(filename, pathway)
