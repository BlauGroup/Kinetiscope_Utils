# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 10:51:48 2024

@author: jacob
"""
import sys
import matplotlib.pyplot as plt
sys.path.append('../common')
from utilities import correct_path_change_dir


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


def plot_concentrations(time_concentration_dict):
    """
    Plots the concentrations of chemical species over time using Matplotlib.

    Args:
    time_concentration_dict (dict): A dictionary where keys are time points and 
                                     values are dictionaries of chemical species 
                                     and their concentrations.
    """
    # Define colors for the top 5 species
    top_colors = ['#57EAEB', '#57BCEB', '#2F5496', '#57EBB9', '#575EEB']
    
    # Get the last time point and its concentrations
    last_time_point = max(time_concentration_dict.keys())
    last_concentrations = time_concentration_dict[last_time_point]

    # Sort the species by their concentrations at the last time point
    top_species = sorted(last_concentrations, key=last_concentrations.get, reverse=True)[:5]
    
    # Prepare for plotting
    times = list(time_concentration_dict.keys())

    # Plot each species
    for i, chemical_name in enumerate(next(iter(time_concentration_dict.values())).keys()):
        # Extract concentration data
        concentrations = [
            time_concentration_dict[time][chemical_name] for time in times
        ]

        # Use the specified color for top species, otherwise default color
        color = top_colors[top_species.index(chemical_name)] if chemical_name in top_species else 'lightgray'
        
        # Plot the concentration vs. time with a bold line for top species
        linewidth = 2 if chemical_name in top_species else 1  # Bold for top species
        plt.plot(times, concentrations, label=chemical_name, color=color, linewidth=linewidth)

    # Add labels and title
    plt.xlabel('Time (s)')
    plt.ylabel('Concentration (M)')
    plt.title('Proton Transfer Reaction Dynamics')
    # plt.legend()  # Show legend
    plt.grid(True)  # Add grid for better readability

    # Show the plot
    plt.show()


kinetiscope_files_dir = (
    r"G:\My Drive\Kinetiscope\production_simulations_092124"
)

correct_path_change_dir(kinetiscope_files_dir)
conc_time_filepath = "proton_transfers_conc_vs_time_10^6_with_excitations.txt"
time_concentration_dict = extract_time_concentrations(conc_time_filepath)
plot_concentrations(time_concentration_dict)
