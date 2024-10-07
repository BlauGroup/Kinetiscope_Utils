# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 10:51:48 2024

@author: jacob
"""
import sys
import matplotlib.pyplot as plt
sys.path.append('../common')
from utilities import correct_path_change_dir
# from compare_simulations import convert_mols_to_particles


class LineNotFoundError(Exception):
    """Custom exception raised when a line with the specified prefix is not
    found."""
    pass


def find_lines_starting_with(file_path, prefix, return_all=False):
    """
    Opens the provided file and finds the first line starting with the given
    prefix. If return_all is True, returns all lines starting from that line
    to the end of the file.

    Args:
    file_path (str): Path to the text file.
    prefix (str): The prefix to search for at the start of a line.
    return_all (bool): If True, returns all lines from the line found
                       to the end of the file.

    Returns:
    list: A list of lines found, either a single line or all lines
          starting from the specified line.

    Raises:
    LineNotFoundError: If no line starts with the specified prefix.
    """

    with open(file_path, "r") as file:
        for line in file:
            if line.strip().startswith(prefix):

                this_line = [line.strip()]

                if return_all:

                    all_lines = (
                        this_line + [other_line.strip() for other_line in file]
                    )

                    return all_lines

                return this_line

    raise LineNotFoundError(
        f"No line starting with '{prefix}' found in the file."
    )


def extract_time_concentrations(file_path):
    """
    Extracts time points and their corresponding concentrations from a
    Kinetiscope data file.

    Args:
    file_path (str): Path to the text file containing the data.

    Returns:
    dict: A dictionary where keys are time points (as floats) and values
          are dictionaries associating chemical names with their
          corresponding concentrations (as floats).

    Raises:
    LineNotFoundError: If no line starting with "1" is found.
    """

    def extract_kinetiscope_names(file_path):
        """
        Opens the provided file, searches for the line starting with "Data",
        and extracts the chemical names, skipping the "Time" and "Data"
        columns.

        Args:
        file_path (str): Path to the text file containing the data.

        Returns:
        list: A list of chemical names (excluding "Time" and "Data").

        Raises:
        LineNotFoundError: If no line starting with "Data" is found.
        """
        data_line_string = find_lines_starting_with(file_path, "Data")[-1]
        data_line_list = data_line_string.split()

        return data_line_list[2:]

    def remove_unnecessary_data(unmodified_time_lines):
        """
        Removes the last two lines from the provided list of time lines
        and splits the remaining lines into lists, removing the first
        element of each.

        Args:
        unmodified_time_lines (list of str): A list of time lines to be
                                               processed.

        Returns:
        list of list: A list of lists, where each inner list contains
                       the split elements of the original lines, excluding
                       the first element.
        """
        end_of_file_removed = unmodified_time_lines[:-2]

        return [string.split()[1:] for string in end_of_file_removed]

    def create_concentration_dict(chemical_names, concentrations):
        """
        Creates a dictionary associating each chemical name with its
        corresponding concentration.

        Args:
        chemical_names (list of str): A list of chemical names.
        concentrations (list of float): A list of concentrations corresponding
                                         to the chemical names.

        Returns:
        dict: A dictionary where keys are chemical names and values are
              their corresponding concentrations (as floats).
        """
        return dict(zip(chemical_names, concentrations))

    chemical_names = extract_kinetiscope_names(file_path)

    time_concentration_dict = {}

    unmodified_time_lines = (
        find_lines_starting_with(file_path, "1", return_all=True)
    )

    corrected_time_lines = remove_unnecessary_data(unmodified_time_lines)

    for time_point in corrected_time_lines:
        time = time_point[0]
        concentrations = time_point[1:]
        concentrations = [
            float(concentration) for concentration in concentrations
        ]

        concentration_dict = (
            create_concentration_dict(chemical_names, concentrations)
        )

        time_concentration_dict[float(time)] = concentration_dict

    return time_concentration_dict


def plot_concentrations(time_concentration_dict, shorthand_to_text_dict,
                        save_path=None, dpi=300):
    """
    Plots the concentrations of chemical species over time using Matplotlib.

    Args:
    time_concentration_dict (dict): A dictionary where keys are time points and
                                    values are dictionaries of chemical species
                                    and their concentrations.
    shorthand_to_text_dict (dict): A dictionary to map shorthand species names
                                   to more descriptive text for the legend.
    save_path (str, optional): Path to save the plot. Defaults to None.
    dpi (int, optional): Resolution of the saved plot. Defaults to 300.
    """
    # Increase figure size
    plt.figure(figsize=(10, 6))  # Width of 10, height of 6

    # Define colors for the top 5 species
    top_colors = ['#2F9696', '#2F7596', '#2F5496', '#2F3496', '#4B2F96']

    # Get the last time point and its concentrations
    last_time_point = max(time_concentration_dict.keys())
    last_concentrations = time_concentration_dict[last_time_point]

    # Sort the species by their concentrations at the last time point
    top_species = sorted(
        last_concentrations, key=last_concentrations.get, reverse=True
    )[:5]

    # Prepare for plotting
    times = list(time_concentration_dict.keys())

    # Plot each species
    for chemical_name in next(iter(time_concentration_dict.values())).keys():
        if chemical_name == "absorption":  # Skip plotting "absorption"
            continue

        # Extract concentration data
        concentrations = [
            time_concentration_dict[time][chemical_name] for time in times
        ]

        # Use the specified color for top species, otherwise default color
        if chemical_name in top_species:
            color = top_colors[top_species.index(chemical_name)]
            # Use shorthand_to_text_dict to map names for the legend
            label = shorthand_to_text_dict.get(chemical_name, chemical_name)
            plt.plot(times, concentrations, label=label, color=color,
                     linewidth=2)  # Bold line for top species
        else:
            plt.plot(times, concentrations, color='lightgray', linewidth=1)

    # Add labels and title
    plt.xlabel('Time (s)')
    plt.ylabel('Events per photon absorbed')
    plt.title('Proton Transfer Reaction Dynamics')

    # Show legend for only top species, position it below the plot
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)

    plt.grid(True)  # Add grid for better readability
    plt.tight_layout()  # Ensure the layout doesn't overlap

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')

    # Show the plot
    plt.show()


def convert_mols_to_particles(mol_dict):
    """Converts a dictionary of moles to particles."""
    particles_per_mol = 2.765e-17

    particle_dict = {
        name: int((1 / particles_per_mol) * moles)
        for name, moles in mol_dict.items()
    }

    return particle_dict


def build_time_particle_dict(time_concentration_dict):
    """Converts concentrations to particle counts for each time point."""
    return {
        time_point: convert_mols_to_particles(conc_dict)
        for time_point, conc_dict in time_concentration_dict.items()
    }


def process_absorption(time_particle_dict):
    """Creates a new dictionary with modified values based on 'absorption'.
       For time points where 'absorption' is 0, the dictionary is not added."""

    new_dict = {}

    for time_point, particle_dict in time_particle_dict.items():

        # For the first time point, just remove 'absorption' and return dict

        if float(time_point) == 0.0:

            new_conc_dict = {key: value for key, value in particle_dict.items()
                             if key != 'absorption'}

            new_dict[time_point] = new_conc_dict

        else:

            absorption_value = particle_dict.pop('absorption', None)

            if absorption_value is None:
                raise ValueError(
                    f"'absorption' key not found for time point {time_point}"
                )

            # Skip adding the time point if absorption is 0 so no div by 0

            if absorption_value == 0:

                continue

            # Create the new concentration dictionary with updated values
            new_conc_dict = {
                key: round(particle_number / absorption_value, 2)
                for key, particle_number in particle_dict.items()
            }

            new_dict[time_point] = new_conc_dict

    return new_dict


if __name__ == "__main__":

    kinetiscope_files_dir = (
        r"G:\My Drive\Kinetiscope\production_simulations_092124"
    )

    correct_path_change_dir(kinetiscope_files_dir)

    conc_time_filepath = (
        "proton_transfers_mols_vs_time_10^6_with_excitations_WithAbs.txt"
    )

    time_concentration_dict = extract_time_concentrations(conc_time_filepath)

    time_particle_dict = build_time_particle_dict(time_concentration_dict)

    shorthand_to_text_dict = {
        "proton_transfer": "all proton transfers",
        "c_n_proton_transfer": "cation transfers proton to neutral",
        "ra_n_proton_transfer": "neutral transfers proton to radical anion",
        "a_c_proton_transfer": "cation transfers proton to anion",
        "a_n_proton_transfer": "neutral transfers proton to anion"
     }

    events_per_photon_dict = process_absorption(time_particle_dict)
    save_path = "excited_proton_transfers_per_photon.png"

    plot_concentrations(
        events_per_photon_dict, shorthand_to_text_dict, save_path
    )
