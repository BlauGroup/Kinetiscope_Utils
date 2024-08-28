# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 09:51:24 2024

@author: jacob
"""

from utilities import correct_path_change_dir
from build_top_reaction_dict import find_top_formation_reactions
import os

"""
This script uses the Pathfinding class's method from HiPRGen to find reaction 
pathways leading to the formation of a specific species in our simulations. It 
follows the highest-frequency reactions that form the species and recursively 
searches for the reactants of each reaction.
"""

def find_reactant_pathways(highest_frequency_reaction, highest_select_dict, starting_species):
    """
    Finds and returns the pathways for each reactant of a given reaction by 
    recursively exploring the pathways of the reactants.

    This function uses recursion to find pathways for each reactant. For each 
    reactant, it calls `find_most_selected_pathway` to get the pathway from the
    reactant to the initial starting species. The recursion continues until it 
    reaches a starting speciesor a pathway with no further reactions.

    Parameters
    ----------
    highest_frequency_reaction : object
        The reaction object containing the reactants and products.
    highest_select_dict : dict
        A dictionary of reactions keyed by product names.
    starting_species : set
        A set containing names of the starting species.
    
    Returns
    -------
    list
        A list of reactant pathways for the reactants.
    """
    
    reactant_pathways = []
    
    for reactant in highest_frequency_reaction.reactants:
        if reactant not in starting_species:
            
            pathway = find_most_selected_pathway(reactant, highest_select_dict, starting_species)
            reactant_pathways.append(pathway)     
            
    return reactant_pathways

def combine_reactant_pathways_with_current_pathway(reactant_pathways, current_pathway):
    """
    Combines all reactant pathways with the current pathway.

    This function appends the current pathway to each reactant pathway. The 
    result is a combined pathway that includes the pathways of the reactants 
    leading up to the current pathway.

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

def final_reaction_forms_reactants(highest_frequency_reaction, pathway):
    """
    Checks if the final reaction in a pathway forms the reactants of the highest
    frequency reaction.

    This function examines the last reaction in the given pathway and compares
    its products with the reactants of the highest frequency reaction. It helps
    to verify if the pathway correctly leads to the highest frequency reaction.

    Parameters
    ----------
    highest_frequency_reaction : object
        The reaction object containing the reactants to be checked.
    pathway : list
        A list of reactions representing the pathway.
    
    Returns
    -------
    bool
        True if the final reaction in the pathway forms the reactants of the 
        highest frequency reaction, False otherwise.
    """
    if not pathway:
        return False

    final_reaction = pathway[-1]
    return highest_frequency_reaction.reactants == final_reaction.products

def find_most_selected_pathway(product, highest_select_dict, starting_species):
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
        A set containing names of the starting species

    Returns
    -------
    list
        The most selected pathway for the given product.
    """
    
    if product not in highest_select_dict:
        raise KeyError("Product does not have a reaction with a nonzero selection frequency.")
    
    current_pathway = []
    highest_frequency_reaction = highest_select_dict[product]
    current_pathway.append(highest_frequency_reaction)
    
    reactant_pathways = (
        find_reactant_pathways(highest_frequency_reaction, highest_select_dict, starting_species)
    )
    
    for pathway in reactant_pathways:
        
        if final_reaction_forms_reactants(highest_frequency_reaction, pathway):
            
            return pathway + current_pathway
    
    final_pathway = (
        combine_reactant_pathways_with_current_pathway(reactant_pathways, current_pathway)
        )
    
    return final_pathway

def save_pathway_to_file(chemical_name, reactions):
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
    
    filename = f"{chemical_name}.txt"
    
    if os.path.exists(filename):
       raise OSError(f"The file '{filename}' already exists.")
    
    with open(filename, 'w') as file:
        
        for rxn in reactions:
            
            file.write(f"{rxn.kinetiscope_name}: {rxn.selection_freq}\n")
            
    print(f"Reaction information saved to {filename}")
             
kinetiscope_files_dir = r"G:\My Drive\Kinetiscope\full_simulations_081224"

correct_path_change_dir(kinetiscope_files_dir)

select_freq_file = "10^5_selection_freqs.txt"
reaction_name_file="10^5_reaction_steps.txt"

top_reaction_dict        = find_top_formation_reactions(
    select_freq_file, 
    reaction_name_file, 
    start_index=7, 
    end_index=1315
)

starting_species = [
    "COO_pMMA_t-butyl_0",
    "SO3_C4F9_-1",
    "phenyl3_S1_+1_#1",
    "COO_CN_C6H4_-1",
    "PHS_CO_C5H5_0_#1"]

starting_species = set(starting_species)

while True:
    
    chemical_name = input("Enter a chemical name (or 'end' to stop): ")
    
    if chemical_name.lower() == 'end':
        
        print("Exiting the script.")
        break
       
    pathway = find_most_selected_pathway(chemical_name, top_reaction_dict, starting_species)
    save_pathway_to_file(chemical_name, pathway)