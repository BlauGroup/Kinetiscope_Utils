# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 15:05:58 2024

@author: jacob
"""

import sys
from plot_kinetiscope_graphs import extract_time_concentrations
import scipy.stats as stats
import os
from monty.serialization import loadfn
from collections import defaultdict
# import time
sys.path.append('../common')
from utilities import correct_path_change_dir
# import matplotlib.pyplot as plt
# import numpy as np
# from matplotlib.ticker import ScalarFormatter


def retain_last_entry(d):
    if d:
        last_key = list(d.keys())[-1]
        last_value = d[last_key]
        return {last_key: last_value}
    return {}


def convert_mols_to_particles(mol_dict):
    particles_per_mol = 2.765e-17

    particle_dict = {
        name: int(1/particles_per_mol * moles) for name, moles in mol_dict.items()
    }

    return particle_dict


def find_top_10_products(product_dict, chemical_names):
    """
    Find the 10 keys with the greatest values in a dictionary, excluding
    specific products.

    Args:
        product_dict (dict): A dictionary where the keys are product names
                             and the values are concentrations or counts.

    Returns:
        list: A list of the 10 keys with the greatest values, excluding
              specified products.
    """
    # Define the products to exclude
    exclude_products = {"PtBMAb_COO_tbut_0", "PHSb_phol_0"}

    # Filter the dictionary to exclude the specified products
    filtered_dict = {
        k: v for k, v in product_dict.items() if k not in exclude_products and k in chemical_names
    }

    top_10_products = sorted(
        filtered_dict, key=filtered_dict.get, reverse=True
    )[:10]

    return top_10_products


# def plot_histogram_comparison(species_list, sim1_dict, sim2_dict):
#     """
#     Plot a histogram comparing the number of particles formed for the top species
#     in two simulations.

#     Args:
#         species_list (list): A list of species to compare (e.g., top 10 concentrated products).
#         sim1_dict (dict): Dictionary of species counts in simulation 1.
#         sim2_dict (dict): Dictionary of species counts in simulation 2.
#     """
#     # Extract the data for species in the list from the two dictionaries
#     sim1_values = [sim1_dict.get(species, 0) for species in species_list]
#     sim2_values = [sim2_dict.get(species, 0) for species in species_list]

#     # Set up the plot
#     bar_width = 0.35  # Width of the bars
#     index = np.arange(len(species_list))  # Index for the x-axis

#     # Create a figure and axis
#     fig, ax = plt.subplots()

#     # Plotting the bars for both simulations
#     bars1 = ax.bar(index, sim1_values, bar_width, label='excitation')
#     bars2 = ax.bar(
#         index + bar_width, sim2_values, bar_width, label='no_excitation'
#     )

#     # Adding labels and title
#     ax.set_xlabel('Species')
#     ax.set_ylabel('Number of Particles')
#     ax.set_title('Comparison of Particle Counts for Top 10 Species')
#     ax.set_xticks(index + bar_width / 2)
#     ax.set_xticklabels(species_list, rotation=45, ha="right")

#     # Add legend
#     ax.legend()

#     # Display the plot
#     plt.tight_layout()  # Adjust layout to fit x-labels
#     plt.show()
def process_all_txt_files(base_dir, subfolder):
    # Change to the base directory
    os.chdir(base_dir)
    
    # Change to the specified subfolder
    os.chdir(subfolder)

    # Get all .txt files in the subfolder
    txt_files = [f for f in os.listdir() if f.endswith('.txt')]

    # Initialize a list to store the dictionaries from each run
    dicts_from_runs = []

    # Loop through each txt file and process it
    for filepath in txt_files:
        # Extract data from the file
        excitation_dict = extract_time_concentrations(filepath)
        
        # Retain the last entry in the dictionary
        final_dict = retain_last_entry(excitation_dict)

        # Add the processed dictionary to the list
        dicts_from_runs.append(final_dict)

    return dicts_from_runs


def check_inner_keys_consistency(list_of_dicts):
    # Get the inner keys from the first dictionary (the inner dictionary of the first float key)
    first_inner_keys = None
    for inner_dict in list_of_dicts[0].values():
        first_inner_keys = set(inner_dict.keys())
        break  # We only need the first set of inner keys
    
    # Compare the inner keys of all other dictionaries to the first one
    for idx, single_dict in enumerate(list_of_dicts[1:], start=1):
        for inner_dict in single_dict.values():
            current_inner_keys = set(inner_dict.keys())
            if current_inner_keys != first_inner_keys:
                # Print the differences in the inner keys
                missing_in_current = first_inner_keys - current_inner_keys
                extra_in_current = current_inner_keys - first_inner_keys
                
                print(f"Mismatch found in dictionary {idx}:")
                if missing_in_current:
                    print(f"  Missing keys in this dictionary: {missing_in_current}")
                if extra_in_current:
                    print(f"  Extra keys in this dictionary: {extra_in_current}")
                return False
    
    print("All dictionaries have consistent inner keys.")
    return True


def calculate_mean_from_dicts(list_of_dicts):
    # Initialize dictionaries to store the sum of values and count of occurrences for each inner key
    summed_values = defaultdict(float)

    # Iterate over the list of dictionaries
    for single_dict in list_of_dicts:
        for inner_dict in single_dict.values():
            # Iterate over the inner dictionary {string: float}
            for string_key, float_value in inner_dict.items():
                # Sum the values for each string key
                summed_values[string_key] += float_value

    # Calculate the mean for each key (we assume there are 3 dictionaries)
    mean_dict = {k: summed_values[k] / len(list_of_dicts) for k in summed_values}

    return mean_dict


kinetiscope_files_dir = (
    r"G:\My Drive\Kinetiscope\production_simulations_092124"
)


files_to_compare = []
excitation_folder = "excitation_amount_files"
excitation_dicts = process_all_txt_files(kinetiscope_files_dir, excitation_folder)
files_to_compare.append(excitation_dicts)
no_excitation_folder = "no_excitation_amount_files"
no_excitation_dicts = process_all_txt_files(kinetiscope_files_dir, no_excitation_folder)
files_to_compare.append(no_excitation_dicts)

check_inner_keys_consistency(excitation_dicts)
check_inner_keys_consistency(no_excitation_dicts)
os.chdir(kinetiscope_files_dir)
chemical_names = loadfn("name_full_mpculeid_092124.json")


dicts_to_compares = []
for list_of_dicts in files_to_compare:
    dicts_to_compares.append(calculate_mean_from_dicts(list_of_dicts))

# correct_path_change_dir(kinetiscope_files_dir)
# excitation_folder = "excitation_amount_files"
# no_excitation_folder = "no_excitation_amount_files"

# excitation_dict = extract_time_concentrations(excitation_filepath)
# no_excitation_dict = extract_time_concentrations(no_excitation_filepath)

# excitation_dict_last = retain_last_entry(excitation_dict)
# no_excitation_dict_last = retain_last_entry(no_excitation_dict)

# excitation_mol_dict = list(excitation_dict_last.values())[0]
# no_excitation_mol_dict = list(no_excitation_dict_last.values())[0]

particle_dicts = []

for dictionary in dicts_to_compares:
    particle_dicts.append(convert_mols_to_particles(dictionary))

excitation_concs = particle_dicts[0]
no_excitation_concs = particle_dicts[1]
# excitation_concs = convert_mols_to_particles(excitation_mol_dict)
# no_excitation_concs = convert_mols_to_particles(no_excitation_mol_dict)

species_set = set(find_top_10_products(particle_dicts[0], chemical_names))

mass_spec_set = set(
    ["C4H8_0_#2", "phyl_H1_0", "t_but_H1_0", "COO_0", "C1H4_0"]
)

species_to_check = species_set.union(mass_spec_set)
# plot_histogram_comparison(
#     species_list, excitation_concs, no_excitation_concs
# )
different_list = []
same_list = []

# total_num_particles = int(1e7)

for chemical, num_particles in excitation_concs.items():
    if chemical in species_to_check:
        no_excitation_conc = no_excitation_concs.get(chemical, None)
        if no_excitation_conc:
            res = (
                stats.poisson_means_test(num_particles, 3, no_excitation_conc, 3)
            )
            if res.pvalue < 0.01:
                print(num_particles)
                print(no_excitation_conc)
                different_list.append(chemical)
            else:
                same_list.append(chemical)

print(f"Total number of different reactions: {len(different_list)}")
print(different_list)
print(f"Total number of same reactions: {len(same_list)}")
print(same_list)
