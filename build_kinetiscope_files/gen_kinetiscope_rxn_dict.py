# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:54:29 2024

@author: JRMilton
"""
__version__ = '1.1.0'

import os
import sys
import csv
from monty.serialization import loadfn, dumpfn
from kinetiscope_reaction_writing_utilities import (
    ReactionDataStorage,
    reclassify_crosslinking_reactions,
    check_species_lengths,
    check_reaction_count
)
from select_ionization_builder import select_ionization_builder
from select_chemical_builder import select_chemical_builder
from correct_names_remove_duplicates import (
    validate_and_correct_reaction_name,
    shorten_PCET
)

"""
This module uses a json associating mpculeids from HiPRGen with human-readable
chemical names, as well as classified HiPRGen_reaction objects (Rxn_classes)
generated in classify_HiPRGen_reactions to create Kinetiscope_reaction objects
associated with each HiPRGen reaction.

Tags based on reaction classifications--e.g. ion-molecule, proton transfer--
will be added as reactants and/or products so we can monitor the production of
marker species in kinetiscope. The module, after creating the reaction objects,
orders them in a user-defined manner, which will be reflected upon import into
kinetiscope, and finally saves information related to the reactions to a csv
file which can be imported into kinetiscope.
"""

class DuplicateFileError(Exception):
    """
    Exception raised for errors related to duplicate files.

    This exception is used when an operation encounters a file that already
    exists when it is expected to be unique.

    Attributes:
    - filename (str): The name of the file that caused the exception.

    Methods:
    - __init__(filename): Initializes the exception with the file name.
    """

    def __init__(self, filename):
        """
        Initializes the DuplicateFileError with a specific filename.

        Parameters:
        - filename (str): The name of the file that already exists.

        The error message is constructed to indicate that the specified file
        already exists. This message is passed to the base Exception class.
        """
        super().__init__(f"The file '{filename}' already exists.")
        self.filename = filename

        
class DuplicateReactionError(Exception):
    """
    Exception raised for errors related to duplicate reactions.

    This exception is used when a reaction is found to be duplicated more
    times than expected.

    Attributes:
    - name (str): The name of the duplicate reaction.
    - count (int): The number of times the reaction was found to be duplicated.

    Methods:
    - __init__(name, count): Initializes the exception with the reaction name
      and count of occurrences.
    """

    def __init__(self, name, count):
        """
        Initializes the DuplicateReactionError with the name of the reaction and
        the count of its occurrences.

        Parameters:
        - name (str): The name of the reaction that is duplicated.
        - count (int): The number of times the reaction appears.

        The error message is constructed to indicate the reaction name and its
        occurrence count. This message is passed to the base Exception class.
        """
        self.name = name
        self.count = count
        super().__init__(f"Duplicate reaction found: '{name}' appears {count} times.")

def remove_duplicate_reactions(kinetiscope_reaction_list):
    """
    Remove duplicate Kinetiscope reactions from a list.

    This function filters out duplicate reactions based on their 
    `kinetiscope_name`. Only the first occurrence of each unique name 
    is kept in the returned list.

    Parameters:
    -----------
    kinetiscope_reaction_list : list
        A list of Kinetiscope reaction objects, each of which should have 
        a `kinetiscope_name` attribute.

    Returns:
    --------
    list
        A new list containing only unique Kinetiscope reactions, with 
        duplicates removed.
    """

    new_list = []
    name_set = set()
    
    for reaction in kinetiscope_reaction_list:
        if reaction.kinetiscope_name not in name_set:
            
            new_list.append(reaction)
            name_set.add(reaction.kinetiscope_name)
            
    return new_list

def count_kinetiscope_names(kinetiscope_reaction_list):
    """
    Count occurrences of each Kinetiscope reaction name in a list.
    
    This function tallies how many times each Kinetiscope reaction name 
    appears in the provided list. If duplicates are found, a 
    DuplicateReactionError is raised for each name with more than one 
    occurrence.
    
    Parameters:
    -----------
    kinetiscope_reaction_list : list
        A list of Kinetiscope reaction objects, each of which should have 
        a `kinetiscope_name` attribute.
    
    Raises:
    -------
    DuplicateReactionError
        If any Kinetiscope reaction name appears more than once in the list.
    
    Returns:
    --------
    None
    """
    
    name_counts = {}


    for reaction in kinetiscope_reaction_list:

        kinetiscope_name = reaction.kinetiscope_name
        
        if kinetiscope_name in name_counts:
            
            name_counts[kinetiscope_name] += 1
            
        else:
            
            name_counts[kinetiscope_name] = 1

    # Check for duplicates and raise an error if any are found
    
    for name, count in name_counts.items():
        
        if count > 1:
            
            raise DuplicateReactionError(name, count)
            
def collect_lists_from_nested_dict(d):
    """
    Collect all list items from a nested dictionary.

    This function traverses a nested dictionary and collects all items
    from lists contained within the dictionary. It recurses through nested
    dictionaries and aggregates items from all lists found.

    Parameters:
    -----------
    d : dict
        The nested dictionary from which lists are to be collected.

    Returns:
    --------
    list
        A list containing all items collected from lists within the nested
        dictionary.
    """
    
    collected_items = []
    recurse_through_dict(d, collected_items)
    
    return collected_items

def recurse_through_dict(d, collected_items):
    """
    Recursively traverse a nested dictionary to collect items from lists.

    This helper function is used by `collect_lists_from_nested_dict` to
    perform a recursive traversal of the dictionary. It extends the list
    of collected items with elements from any lists found in the dictionary.

    Parameters:
    -----------
    d : dict
        The nested dictionary to traverse.
    collected_items : list
        The list to which items from nested lists are added.
    """
    
    for value in d.values():
        
        if isinstance(value, dict):
            
            recurse_through_dict(value, collected_items)
            
        elif isinstance(value, list):
            
            collected_items.extend(value)

def get_supercategory_index(reaction, supercategory_order):
    """
    Get the index of the supercategory for a reaction.

    This function determines the index of the supercategory of a reaction
    based on its marker species and the provided supercategory order.

    Parameters:
    -----------
    reaction : object
        The Kinetiscope reaction object to check.
    supercategory_order : list
        A list defining the order of supercategories.

    Returns:
    --------
    int
        The index of the supercategory, or a large number if not found.
    """
    
    for index, supercategory in enumerate(supercategory_order):
        if supercategory in reaction.marker_species:
            
            return index
        
    return len(supercategory_order)

def get_subcategory_index(
        reaction, supercategory_order, 
        supercategories_with_subcategories, subcategory_order
):
    """
    Get the index of the subcategory for a reaction.

    This function determines the index of the subcategory of a reaction
    based on its marker species and the provided subcategory order, only if
    the reaction's supercategory has associated subcategories.

    Parameters:
    -----------
    reaction : object
        The Kinetiscope reaction object to check.
    supercategory_order : list
        A list defining the order of supercategories.
    supercategories_with_subcategories : list
        A list of supercategories that have associated subcategories.
    subcategory_order : list
        A list defining the order of subcategories.

    Returns:
    --------
    int
        The index of the subcategory, or a large number if not found.
    """
    for supercategory in supercategory_order:
        
        supercategory_present = supercategory in reaction.marker_species
        
        is_supercategory_with_subcategories = (
        supercategory in supercategories_with_subcategories
        )
        
        if supercategory_present and is_supercategory_with_subcategories:
            
            for index, subcategory in enumerate(subcategory_order):
                if subcategory in reaction.marker_species:
                    
                    return index
                
    return len(subcategory_order)

def get_sort_key(
        reaction, supercategory_order, 
        supercategories_with_subcategories, subcategory_order
):
    """
    Generate a sort key for a Kinetiscope reaction based on its supercategory 
    and subcategory.
    
    This function calculates indices for both supercategory and subcategory 
    of the reaction. These indices are used to sort reactions by their 
    supercategory and subcategory in the specified orders.
    
    Parameters:
    -----------
    reaction : object
        The Kinetiscope reaction object for which the sort key is generated.
    supercategory_order : list
        A list specifying the order of supercategories.
    supercategories_with_subcategories : list
        A list of supercategories that include subcategories.
    subcategory_order : list
        A list specifying the order of subcategories.
    
    Returns:
    --------
    tuple
        A tuple containing the supercategory index and subcategory index 
        used for sorting the reaction.
    """
    supercategory_index = (
        get_supercategory_index(reaction, supercategory_order)
    )
    
    subcategory_index = get_subcategory_index(
        reaction, 
        supercategory_order, 
        supercategories_with_subcategories, 
        subcategory_order
    )
    
    return supercategory_index, subcategory_index

def order_kinetiscope_reactions(
    kinetiscope_reactions, supercategory_order, 
    supercategories_with_subcategories, subcategory_order
):
    """
    Sort a list of Kinetiscope reactions based on supercategory and 
    subcategory.

    This function orders a list of Kinetiscope reactions first by 
    supercategory and then by subcategory. The sorting is determined 
    using the provided orders for supercategories and subcategories.

    Parameters:
    -----------
    kinetiscope_reactions : list
        A list of Kinetiscope reaction objects to be sorted.
    supercategory_order : list
        A list specifying the order of supercategories.
    supercategories_with_subcategories : list
        A list of supercategories that include subcategories.
    subcategory_order : list
        A list specifying the order of subcategories.

    Returns:
    --------
    list
        The list of Kinetiscope reactions, sorted by supercategory 
        and subcategory.
    """
    # Sort the reactions using the key provided by get_sort_key
    ordered_reactions = sorted(
        kinetiscope_reactions,
        key=lambda reaction: get_sort_key(
            reaction, 
            supercategory_order, 
            supercategories_with_subcategories, 
            subcategory_order
        )
    )
    
    return ordered_reactions

def create_list_for_csv(ordered_reactions):
    """
    Each reaction, when we import it to kinetiscope, has a lot of information
    and/or numbers associated with it, so we create a dictionary associated
    with each reaction to write those numbers to the incipient excel file.

    Parameters
    ----------
    ordered_reactions : list
        list of reactions, which we have already ordered

    Returns
    -------
    dict_list : list
        list of dictionaries associated with each reaction. Has the same order
        as ordered_reactions

    """
    
    dict_list = []
    
    for reaction in ordered_reactions:
        
        csv_dict = {}
        
        rate_coefficient_format = (
            3 if "absorption" in reaction.marker_species else 0
        )
        
        csv_dict['# equation'] = reaction.kinetiscope_name
        csv_dict['fwd_A'] = 1
        csv_dict['fwd_temp_coeff'] = 0
        csv_dict['fwd_Ea'] = 0
        
        csv_dict['fwd_k'] = (
            reaction.rate_coefficient if "absorption" not in reaction.marker_species else 1
        )
        
        csv_dict['rev_A'] = 1
        csv_dict['rev_M'] = 0
        csv_dict['rev_Ea'] = 0
        csv_dict['rev_k'] = 1
        csv_dict['fwd_k0'] = 1
        csv_dict['rev_k0'] = 1
        csv_dict['alpha_val'] = 0.5
        csv_dict['equil_potential'] = 0
        csv_dict['num_electrons'] = 0
        
        csv_dict['fwd_prog_k'] = (
            reaction.rate_coefficient if "absorption" in reaction.marker_species else 1
        )
        
        csv_dict['rev_prog_k'] = 1
        csv_dict['non_stoichiometric'] = 0
        csv_dict['rate_constant_format'] = rate_coefficient_format
        
        dict_list.append(csv_dict)
    
    return dict_list

def check_and_raise_if_duplicate(filename):
    """
    Check if a file with the given name exists in the current directory.
    Raises DuplicateFileError if the file already exists.

    Parameters:
    -----------
    filename : str
        The name of the file to check.

    Raises:
    -------
    DuplicateFileError
        If the file exists in the current directory.
    """
    if os.path.isfile(filename):
        raise DuplicateFileError(filename)
        
def write_reactions_to_json(dict_list, new_filename):
    """
    Writes a list of reaction dictionaries to a JSON file.

    The function first checks if a file with the given filename already exists
    using the `check_and_raise_if_duplicate` function. If the file exists, it
    raises a `DuplicateFileError` and exits the program. Otherwise, it writes
    the reaction data to the JSON file, creating a new file with the specified
    `new_filename`.

    Parameters:
    - dict_list (list of dict): A list of dictionaries, where each dictionary
      represents a reaction and contains the data to be written to the JSON 
      file.
    - new_filename (str): The name of the file to write the JSON data to.

    Raises:
    - DuplicateFileError: If a file with the same name already exists.

    Example:
    ```
    dict_list = [{'equation': 'H2 + O2 -> H2O', 'fwd_A': 1, ...}, ...]
    new_filename = 'reactions.json'
    write_reactions_to_json(dict_list, new_filename)
    """
    
    try:
        # Check if the file already exists
        check_and_raise_if_duplicate(new_filename)
        print(f"The file '{new_filename}' does not exist. Safe to proceed.")
        
    except DuplicateFileError as e:
        print(e)
        sys.exit(1)

    dumpfn(dict_list, new_filename)

def write_reactions_to_csv(dict_list, new_filename):
    """
    Writes a list of reaction dictionaries to a CSV file.
    
    The function first checks if a file with the given filename already exists
    using the `check_and_raise_if_duplicate` function. If the file exists, it
    raises a `DuplicateFileError` and exits the program. Otherwise, it writes
    the reaction data to the CSV file, creating a new file with the specified
    `new_filename`.
    
    Parameters:
    - dict_list (list of dict): A list of dictionaries, where each dictionary
      represents a reaction and contains the data to be written to the CSV file.
    - new_filename (str): The name of the file to write the CSV data to.

    Raises:
    - DuplicateFileError: If a file with the same name already exists.
    
    Example:
    ```
    dict_list = [{'# equation': 'H2 + O2 -> H2O', 'fwd_A': 1, ...}, ...]
    new_filename = 'reactions.csv'
    write_reactions_to_csv(dict_list, new_filename)
    ```
    """
    
    try:
        # Check if the file already exists
        check_and_raise_if_duplicate(new_filename)
        print(f"The file '{new_filename}' does not exist. Safe to proceed.")
        
    except DuplicateFileError as e:
        print(e)
        sys.exit(1)
    
    # Write the list of dictionaries to a CSV file
    with open(new_filename, 'w', newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=dict_list[0].keys())
        writer.writeheader()
        
        for reaction in dict_list:
            writer.writerow(reaction)

kinetiscope_reaction_list = []
os.chdir("G:/My Drive/Kinetiscope/production_simulations_092124")
# test_rxns = "HiPRGen_rxns_to_name.json"
# HiPRGen_reaction_list = loadfn(test_rxns)
full_rxns = "HiPRGen_rxns_to_name_full_092124.json"
HiPRGen_reaction_list = loadfn(full_rxns)
# name_mpculeid_file = "name_test_mpculeid_080624.json"
name_mpculeid_file = "name_full_mpculeid_092124.json"
name_mpculeid_dict = loadfn(name_mpculeid_file)
# for name in name_mpculeid_dict.keys():
#     print(name)

absorption_rate_constants = {
    "4864aee73a83d357c31fadd50b81e3cd-C10H20O2-0-1":1.4E-01,
    "00a7dcc352b0d613f58e850935bf5609-C10H14O1-0-1":1.0E-01,
    "bfed458e642b8daa1eab6bc02d5e5682-C18H15S1-1-1":1.8E-01,
    "9a8a88b8b92c714d7f65b8526ffabc7a-C4F9O3S1-m1-1":5.8E-01,
    "17f31f89123edbaa0e3b9c7eb49d26f3-C8H4N1O2-m1-1":1.5E-01
}
    
marker_species_dict = {
    "radical_cation":"rc",
    "radical_anion":"ra",
    "neutral_radical":"nr",
    "cation":"c",
    "anion":"a",
    "neutral":"n",
    "proton_coupled_electron_transfer":"PCET"
}

excitation_set = set()

#we store this in a ReactionDataStorage object, described in 
# kinetiscope_reaction_writing_utilities so that we can pass all of this stuff
#as a single arguement to functions

reaction_writing_data = ReactionDataStorage(
    name_mpculeid_dict,
    marker_species_dict,
    excitation_set,
    absorption_rate_constants
)

HiPRGen_ionization_reactions = HiPRGen_reaction_list["ionization"].values()

#generate kinetiscope reactions for ionization reacctions

for rxn_list in HiPRGen_ionization_reactions:
    
    for HiPRGen_rxn in rxn_list:
        
        ionization_reaction_list = \
            select_ionization_builder(HiPRGen_rxn, reaction_writing_data) 
            
        #ionization_reaction has >=1 elements
        
        kinetiscope_reaction_list.extend(ionization_reaction_list)

chemical_reaction_list = (
    collect_lists_from_nested_dict(HiPRGen_reaction_list["chemical"])
)

chemical_reaction_list = (
    collect_lists_from_nested_dict(HiPRGen_reaction_list["chemical"])
)

#we need the kinetiscope names to do this, but just do it for the HiPRGen
#reactions so we can write the names with these tags later

reclassify_crosslinking_reactions(chemical_reaction_list, reaction_writing_data)

# repeat for chemical reactions

for HiPRGen_rxn in chemical_reaction_list:

    chemical_reaction_list, reaction_writing_data = (
        select_chemical_builder(HiPRGen_rxn, reaction_writing_data)
    )
    
    #chemical_reaction_list has >=1 elements
    
    kinetiscope_reaction_list.extend(chemical_reaction_list)

#correct names of reactions when they contain duplicate species and a marker
#species that is way too long

for index, reaction in enumerate(kinetiscope_reaction_list):
    
    reaction_has_duplicate_species, corrected_reaction = (
        validate_and_correct_reaction_name(reaction)
    )
    
    if reaction_has_duplicate_species:
        
        kinetiscope_reaction_list[index] = corrected_reaction
        
    if "proton_coupled_electron_transfer" in reaction.marker_species:
        
        kinetiscope_reaction_list[index] = shorten_PCET(reaction)
 
#make sure each reaction in the list is unique
    
try:
    
    kinetiscope_reaction_list = (
        remove_duplicate_reactions(kinetiscope_reaction_list)
    )
    
    count_kinetiscope_names(kinetiscope_reaction_list)
    
except DuplicateReactionError as e:
    print(f"Error: {e}")


barrierless_list = [
    "absorption", "electron_ionization", "recombination", "attachment",
    "excitation", "dexcitation", "fragmentation", "isomerization", "ion-ion",
    "proton_transfer", "H_atom_abstraction", "hydride_abstraction",
    "electron_transfer", "crosslinking", "combination"]

print(f"length of reaction list before pruning: {len(kinetiscope_reaction_list)}")
barrierless_reaction_list = []
for reaction in kinetiscope_reaction_list:

    reaction_is_barrierless = (
        any(elem in reaction.marker_species for elem in barrierless_list)
    )

    if reaction_is_barrierless:
        barrierless_reaction_list.append(reaction)
    else:
        try:
            assert "reaction" in reaction.marker_species or "proton_coupled_electron_transfer" in reaction.marker_species
        except:
            print(reaction.marker_species)
            sys.exit()
    
print(f"length of reaction list after pruning: {len(barrierless_reaction_list)}")
sys.exit()

supercategory_order = [
    "absorption", "electron_ionization", "recombination", "attachment", 
    "excitation", "dexcitation", "crosslinking", "fragmentation", 
    "isomerization", "ion-ion", "ion-molecule", "neutral"
]

supercategories_with_subcategories = (
    set(["ion-ion", "ion-molecule", "neutral"])
)

subcategory_order = [
    "proton_transfer", "H_atom_abstraction", "hydride_abstraction", 
    "proton_coupled_electron_transfer", "electron_transfer", "reaction"
    ]

# order reactions based on all of these categories, to later parse through our
# excel file

ordered_reactions = order_kinetiscope_reactions(
    kinetiscope_reaction_list, 
    supercategory_order, 
    supercategories_with_subcategories, 
    subcategory_order
)

for reaction in ordered_reactions:

    check_species_lengths(reaction.kinetiscope_name)
    check_reaction_count(reaction.kinetiscope_name)

parent_filename = "barrierless_only_101424"
json_filename = parent_filename + ".json"

write_reactions_to_json(ordered_reactions, json_filename)

print('Writing reactions to csv file...')

# create a dict associated with each reaction for easy saving

list_for_csv = create_list_for_csv(ordered_reactions)

new_csv_filename = parent_filename + ".csv"

#write those reactions to a csv file, which can be imported into kinetiscope

write_reactions_to_csv(list_for_csv, new_csv_filename)
  
print('Done!') 