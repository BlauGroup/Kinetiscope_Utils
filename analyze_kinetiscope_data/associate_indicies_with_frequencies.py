# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 11:26:28 2024

@author: jacob
"""

def extract_numbers_from_line(line):
    """
    Extracts integer numbers from a line of text.

    Parameters
    ----------
    line : str
        The line of text containing numbers.

    Returns
    -------
    list
        A list of integers extracted from the line.
    """
    
    return [int(string) for string in line.split() if string.isdigit()]

def find_selection_freqs(all_lines, last_point_index):
    """
    Extracts selection frequencies from a specific line in a file.

    Parameters
    ----------
    all_lines : list
        A list containing each line in the file as a list item.
    last_point_index : int
        The index of the line containing selection frequencies.

    Returns
    -------
    list
        A list of selection frequencies.
    """
    
    last_line = all_lines[last_point_index]
    point_and_data_removed = last_line.split()[2:]
    string_to_process = " ".join(point_and_data_removed)
    
    return extract_numbers_from_line(string_to_process)

def find_reaction_indicies(all_lines, reaction_line_index):
    """
    Extracts reaction indices from a specific line in a file.

    Parameters
    ----------
    all_lines : list
        A list containing each line in the file as a list item.
    reaction_line_index : int
        The index of the line containing reaction indices.

    Returns
    -------
    list
        A list of reaction indices.
    """
    
    reaction_line = all_lines[reaction_line_index]
    reaction_line = reaction_line.replace("Step", "")
    
    return extract_numbers_from_line(reaction_line)
    
def convert_to_indicies(reaction_line, last_point_line):
    """
    Helper function to convert indicies in the kinetiscope file to zero-based
    indicies

    Parameters
    ----------
    reaction_line : int
        the line number in the file that has the step number, where each "step"
        represents a chemical reaction
    last_point_line : int
        the line of the last time point in the simulation

    Returns
    -------
    int
        reaction_line corrected to zero indexing
    int
        last_point_line corrected to zero indexing

    """
    return reaction_line-1, last_point_line-1  

def find_numbers_and_frequencies(all_lines, reaction_line, last_point_line):
    """
    A wrapper function that creates two lists, reaction_inicies and
    selection_freqs, which will be associated with each other in a dictionary.

    Parameters
    ----------
    all_lines : string
        a list, generated using the .readlines() method, containing all lines
        in a file
    reaction_line : int
       denotes the line of a file where we'll find reaction 
        indicies
    last_point_line : int
        denotes the line where the last recorded data point in the simulation
        is

    Returns
    -------
    reaction_indicies : list
        list of integers of reaction indicies
    selection_freqs : list
        list of selection frequencies

    """
    
    reaction_line_index, last_point_index = (
        convert_to_indicies(reaction_line, last_point_line)
    )
    
    reaction_indicies = find_reaction_indicies(all_lines, reaction_line_index)
    selection_freqs = find_selection_freqs(all_lines, last_point_index)
    
    return reaction_indicies, selection_freqs

def check_for_mismatched_lengths(reaction_indices, selection_freqs):
    """
    Checks if the lengths of the reaction indices and selection frequencies lists match.

    Parameters
    ----------
    reaction_indices : list
        List of reaction indices.
    selection_freqs : list
        List of selection frequencies.

    Raises
    ------
    ValueError
        If the lengths of the reaction indices and selection frequencies lists do not match.
    """
    if len(reaction_indices) != len(selection_freqs):
        print("Mismatch in the number of step numbers and selection frequencies.")
        print(f"Reaction indices: {reaction_indices}")
        print(f"Selection frequencies: {selection_freqs}")
        raise ValueError("Mismatch in the number of step numbers and selection frequencies.")

def check_line_ranges(reaction_line, last_point_line, num_lines):
    """
    Checks if the given line numbers are within the valid range of the file.

    Parameters
    ----------
    reaction_line : int
        The line number for reaction indices.
    last_point_line : int
        The line number for selection frequencies.
    num_lines : int
        The total number of lines in the file.

    Raises
    ------
    IndexError
        If either line number is out of range.
    """
    if reaction_line > num_lines:
        raise IndexError(f"Reaction line number {reaction_line} is out of range. The file has {num_lines} lines.")
    if last_point_line > num_lines:
        raise IndexError(f"Last point line number {last_point_line} is out of range. The file has {num_lines} lines.")
        
def build_index_freq_dict(select_freq_file, reaction_line, last_point_line):
    """
    Builds a dictionary mapping reaction indices to their corresponding 
    selection frequencies from a file.

    Parameters
    ----------
    select_freq_file : str
        The path to the file containing the selection frequencies and reaction indices.
    reaction_line : int
        The line number in the file where reaction indices are listed.
    last_point_line : int
        The line number in the file where the selection frequencies are listed.

    Returns
    -------
    dict
        A dictionary where keys are reaction indices and values are their 
        corresponding selection frequencies.

    Raises
    ------
    IndexError
        If reaction_line or last_point_line is out of range for the file.
    ValueError
        If there is a mismatch between the number of reaction indices and selection frequencies.
    """
    
    with open(select_freq_file, 'r') as file:
        all_lines = file.readlines()

    num_lines = len(all_lines)  
    check_line_ranges(reaction_line, last_point_line, num_lines)
    reaction_indicies, selection_freqs = find_numbers_and_frequencies(all_lines, reaction_line, last_point_line) 
    check_for_mismatched_lengths(reaction_indicies, selection_freqs)
    index_freq_dict = dict(zip(reaction_indicies, selection_freqs))

    return index_freq_dict