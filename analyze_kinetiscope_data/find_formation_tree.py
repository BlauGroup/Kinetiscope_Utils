# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 09:51:24 2024

@author: jacob
"""

import os
import sys
from Rxn_classes import Kinetiscope_simulation_reaction

def remove_whitespace(dirty_name):
    name_without_left_space = dirty_name.lstrip()
    return name_without_left_space.rstrip()

def find_reaction_name(line):
    split_line = line.split(":")
    dirty_name = split_line[1]
    return remove_whitespace(dirty_name)

def build_k_sim_reaction(line, index, selection_freq):
    name = find_reaction_name(line)
    return Kinetiscope_simulation_reaction(name, index, selection_freq)

def extract_reactions(reaction_name_file, reaction_dict, index_freq_dict):
    index = 1
    
    with open(reaction_name_file, 'r') as file:
        for line in file:
            if "Equation" in line:
                selection_freq = index_freq_dict.get(int(index), False)
                k_sim_reaction = build_k_sim_reaction(line, index, selection_freq)
                reaction_dict[index] = k_sim_reaction
                index += 1
    return reaction_dict

def find_selection_freqs(all_lines, last_point_index):
    last_point_line_content = all_lines[last_point_index].strip().split()
    selection_freqs = last_point_line_content[2::] #first two columns are point # and time
    selection_freqs = [int(selection_freq) for selection_freq in selection_freqs]
    return selection_freqs

def find_step_numbers(all_lines, step_index):
    step_numbers = []
    
    step_line_content = all_lines[step_index].split()
    
    for title in step_line_content:
        if title.isdigit():
            step_number = int(title)
            step_numbers.append(step_number)
            
    return step_numbers
    
def convert_to_indicies(step_line, last_point_line):
    return step_line-1, last_point_line-1  

def find_numbers_and_frequencies(all_lines, step_line, last_point_line):
    step_index, last_point_index = convert_to_indicies(step_line, last_point_line)
    step_numbers = find_step_numbers(all_lines, step_index)
    selection_freqs = find_selection_freqs(all_lines, last_point_index)
    return step_numbers, selection_freqs

def build_index_freq_dict(select_freq_file, step_line, last_point_line):

    with open(select_freq_file, 'r') as file:
        all_lines = file.readlines()
        
    if step_line > len(all_lines) or last_point_line > len(all_lines):
        raise IndexError("One or both of the line numbers are out of range.")
        
    step_numbers, selection_freqs = (
        find_numbers_and_frequencies(all_lines, step_line, last_point_line)
    )
    
    if len(step_numbers) != len(selection_freqs):
       raise ValueError("Mismatch in the number of step numbers and selection frequencies.")
    
    index_freq_dict = dict(zip(step_numbers, selection_freqs))

    return index_freq_dict

def build_highest_select_dict(reaction_dict):
    highest_select_dict = {}
    
    for reaction in reaction_dict.values():
        products = reaction.products
        for product in products:
            if product not in highest_select_dict and reaction.selection_freq != 0:
                highest_select_dict[product] = reaction
            elif product not in highest_select_dict and reaction.selection_freq == 0:
                pass
            else:
                current_most_frequent_reaction = highest_select_dict.get(product, False)
                if reaction.selection_freq > current_most_frequent_reaction.selection_freq:
                    highest_select_dict[product] = reaction
                    
    return highest_select_dict

def find_most_selected_pathway(product, highest_select_dict, starting_species):
    
    if product not in highest_select_dict:
        raise KeyError("product does not have a reaction with a nonzero selection frequency")
    
    pathway = []
    highest_frequency_reaction = highest_select_dict[product]
    pathway.append(highest_frequency_reaction)
    
    prefixes = []
    for reactant in highest_frequency_reaction.reactants:
        if reactant not in starting_species:
            prefix = find_most_selected_pathway(reactant, highest_select_dict, starting_species)
            prefixes.append(prefix)     
            prefix_final_reaction = prefix[-1]
            
            if highest_frequency_reaction.reactants == prefix_final_reaction.products:
                return prefix + pathway
    
    for prefix in prefixes:
        pathway = prefix + pathway
    return pathway
              
kinetiscope_files_dir = r"G:\My Drive\Kinetiscope\full_simulations_081224"
corrected_path = os.path.normpath(kinetiscope_files_dir)
os.chdir(corrected_path)

select_freq_file = "10^5_selection_freqs.txt"

index_freq_dict = build_index_freq_dict(select_freq_file, 7, 1315)

reaction_dict = {}
reaction_name_file = "10^5_reaction_steps.txt"  # Replace with your actual file name

reaction_dict = extract_reactions(reaction_name_file, reaction_dict, index_freq_dict)
highest_select_dict = build_highest_select_dict(reaction_dict)
starting_species = [
    "COO_pMMA_t-butyl_0",
    "SO3_C4F9_-1",
    "phenyl3_S1_+1_#1",
    "COO_CN_C6H4_-1",
    "PHS_CO_C5H5_0_#1"]

starting_species = set(starting_species)

while True:
    # Query user for a chemical name
    chemical_name = input("Enter a chemical name (or 'end' to stop): ")
    
    # Check if the user wants to end the loop
    if chemical_name.lower() == 'end':
        print("Exiting the script.")
        break
    
    test = find_most_selected_pathway(chemical_name, highest_select_dict, starting_species)
    for rxn in test:
        print(rxn.kinetiscope_name)
        print(rxn.selection_freq)
    # Check if the chemical name is in the dictionary
    # if chemical_name in chemical_dict:
    #     print(f"The mpculeid for {chemical_name} is {chemical_dict[chemical_name]}.")
    # else:
    #     print(f"Error: {chemical_name} not found in the dictionary.")
        

