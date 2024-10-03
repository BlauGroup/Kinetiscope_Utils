# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 15:05:58 2024

@author: jacob
"""

import sys
from plot_kinetiscope_graphs import extract_time_concentrations
import scipy.stats as stats
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
    # Define Avogadro's number
    particles_per_mol = 2.765e-17

    particle_dict = {
        name: int(1/particles_per_mol * moles) for name, moles in mol_dict.items()
    }

    return particle_dict


def find_top_10_products(product_dict):
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
    # exclude_products = {"PtBMAb_COO_tbut_0", "PHSb_phol_0"}

    # Filter the dictionary to exclude the specified products
    filtered_dict = {
        k: v for k, v in product_dict.items()
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


kinetiscope_files_dir = (
    r"G:\My Drive\Kinetiscope\production_simulations_092124"
)

correct_path_change_dir(kinetiscope_files_dir)
excitation_filepath = "euvl_withexcitation_mol_conc_time.txt"
no_excitation_filepath = "euvl_noexcitation_mol_conc_time.txt"

excitation_dict = extract_time_concentrations(excitation_filepath)
no_excitation_dict = extract_time_concentrations(no_excitation_filepath)

excitation_dict_last = retain_last_entry(excitation_dict)
no_excitation_dict_last = retain_last_entry(no_excitation_dict)

excitation_mol_dict = list(excitation_dict_last.values())[0]
no_excitation_mol_dict = list(no_excitation_dict_last.values())[0]

excitation_concs = convert_mols_to_particles(excitation_mol_dict)
no_excitation_concs = convert_mols_to_particles(no_excitation_mol_dict)

species_set = set(find_top_10_products(excitation_concs))

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
            res = stats.poisson_means_test(num_particles, 1, no_excitation_conc, 1)
            if res.pvalue < 0.01:
                different_list.append(chemical)
            else:
                same_list.append(chemical)

print(f"Total number of different reactions: {len(different_list)}")
print(f"Total number of same reactions: {len(same_list)}")
