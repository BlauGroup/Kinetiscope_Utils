# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 10:42:47 2024

@author: jacob
"""
import re

def shorten_PCET(reaction):
    """
    The full name, proton_coupled_electron_transfer, may be too long for
    Kinetiscope. This function just alters any name containing that phrase with
    the shorthand "PCET."

    Parameters
    ----------
    reaction : kinetiscope reaction object
        the reaction whose name we're modifying

    Returns
    -------
    reaction : kinetiscope reaction object
        the reaction with an updated name

    """
    reaction.kinetiscope_name = (
        reaction.kinetiscope_name.replace("proton_coupled_electron_transfer", "PCET")
    )

    return reaction

def replace_species_with_2_species(species_list):
    """
    Replace a single species with a notation indicating two species.
    
    This function takes a list containing a single species and returns a
    string that denotes two of that species. The notation used is "2 " 
    followed by the species name.
    
    Parameters:
    -----------
    species_list : list
        A list containing a single species. The function will use this 
        species to create a notation indicating two of that species.
    
    Returns:
    --------
    str
        A string indicating two of the species, formatted as "2 <species>".
    """
    
    species = species_list[0]
    two_species = "2 " + species
    
    return two_species

def star_test(species_list):
    """
    Check if any species in the list contains a star character.

    This function examines a list of species to determine if any of the
    species contain the star character ("*"). It returns `True` if at least
    one species contains the star character and `False` otherwise.

    Parameters:
    -----------
    species_list : list
        A list of species to be checked for the presence of the star
        character.

    Returns:
    --------
    bool
        `True` if any species in the list contains the star character;
        otherwise, `False`.
    """

    for species in species_list:

        if "*" in species:

            return True

    return False

def replace_second_star_reactant(reactant_list):
    """
    Replace the second reactant with its star character removed. 
    
    This function takes a list of reactants and removes the star character
    ("*") from the second reactant in the list. The function modifies the 
    list in place and returns the updated list.
    
    Parameters:
    -----------
    reactant_list : list
        A list of reactants where the second reactant will have its star 
        character removed.
    
    Returns:
    --------
    list
        The updated list of reactants with the star character removed from 
        the second reactant.
    """

    reactant_list[1] = reactant_list[1].replace("*", "")

    return reactant_list

def correct_name_reactants(kinetiscope_reaction, reactant_list):
    """
    Modify the Kinetiscope reaction name by updating either the reactants or
    products.
    
    This function updates the `kinetiscope_name` of a Kinetiscope reaction 
    object by replacing the reactants or products side of the reaction name 
    with a new list of species. The reaction name is normalized by removing any
    extra whitespace around the "=>" separator. The function joins the species
    list with " + " and updates the appropriate side of the reaction name based
    on the provided `side_identifier`.
    
    Parameters:
    -----------
    kinetiscope_reaction : object
        An object representing a Kinetiscope reaction, which has an attribute
        `kinetiscope_name` 
        that will be modified.
    species_list : list
        A list of species to be used for updating either the reactants or 
        products side of the reaction name.
    side_identifier : str
        A string indicating which side of the reaction to update. Should be 
        "reactants" to update the reactants side or any other value to update 
        the products side.
    
    Returns:
    --------
    str
        The modified Kinetiscope reaction name with the updated side (reactants
        or products) and normalized whitespace.
    """
    
    if star_test(reactant_list):
        
        reactant_list = replace_second_star_reactant(reactant_list)
        
    else:
        
        reactant_list = [replace_species_with_2_species(reactant_list)]
    
    kinetiscope_reaction.kinetiscope_name = (
        modify_kinetiscope_name(kinetiscope_reaction, reactant_list, "reactants")
    )
    
    return kinetiscope_reaction

def modify_kinetiscope_name(kinetiscope_reaction, species_list, side_identifier):
    
    # Normalize whitespace around "=>"
    
    sides = re.split(r"\s*=>\s*", kinetiscope_reaction.kinetiscope_name)
    
    if side_identifier == "reactants":
        
        sides[0] = " + ".join(species_list)
        
    else:
        
        sides[1] = " + ".join(species_list)
    
    return " => ".join(sides)

def correct_name_products(kinetiscope_reaction, product_list):
    """
    Update the Kinetiscope reaction name with corrected product names.
    
    This function takes a Kinetiscope reaction object and a list of products, 
    processes the product list to ensure correct naming, and updates the 
    Kinetiscope reaction's name with the corrected product names. The function
    replaces the product list with a version where each species is replaced by 
    a designation indicating two species if needed, and then modifies the 
    reaction name accordingly.
    
    Parameters:
    -----------
    kinetiscope_reaction : object
        An object representing a Kinetiscope reaction, which should have an 
        attribute `kinetiscope_name` 
        to be updated.
    product_list : list
        A list of products to be processed and used for updating the reaction
        name.
    
    Returns:
    --------
    object
        The updated Kinetiscope reaction object with the corrected product 
        names applied to its `kinetiscope_name`.
    """
    
    product_list = [replace_species_with_2_species(product_list)]
    
    kinetiscope_reaction.kinetiscope_name = (
        modify_kinetiscope_name(kinetiscope_reaction, product_list, "products")
    )
    
    return kinetiscope_reaction

def find_reactants_and_products(kinetiscope_reaction):
    """
    Extract and clean reactants and products from a Kinetiscope reaction name.
    
    This function takes a Kinetiscope reaction object, extracts the reactants
    and products from its name, and returns cleaned lists of reactants and 
    products. The lists are cleaned by removing any "+" characters and limiting
    the product list to the first three elements.
    
    Parameters:
    -----------
    kinetiscope_reaction : object
        An object representing a Kinetiscope reaction, which should have an 
        attribute `kinetiscope_name` containing the reaction name in the format
        "reactants => products".
    
    Returns:
    --------
    tuple
        A tuple containing two lists:
        - `corrected_reactant_list`: A list of reactants with "+" characters
        removed.
        - `corrected_product_list`: A list of products with "+" characters 
        removed and limited to the first three elements.
    """

    kinetiscope_reaction_name = kinetiscope_reaction.kinetiscope_name
    
    reactant_str, product_str = kinetiscope_reaction_name.split("=>")
    
    reactant_list = reactant_str.split()
    product_list = product_str.split()[:3] #ignores marker species
    
    corrected_reactant_list = [s for s in reactant_list if s != "+"]
    corrected_product_list = [s for s in product_list if s != "+"]
        
    return corrected_reactant_list, corrected_product_list

def correct_if_duplicates(kinetiscope_reaction, species_list, correction_function):
    """
    Apply a correction to the Kinetiscope reaction if duplicate species are 
    detected.
    
    This function checks if the given list of species contains duplicates. 
    If duplicates are found, it applies a specified correction function to the 
    Kinetiscope reaction. If no duplicates are found, the original reaction is
    returned unchanged.
    
    Parameters:
    -----------
    kinetiscope_reaction : object
        An object representing a Kinetiscope reaction, which will be modified 
        if duplicates are detected.
    species_list : list
        A list of species to be checked for duplicates.
    correction_function : function
        A function that applies the necessary correction to the Kinetiscope 
        reaction if duplicates are found. It takes the Kinetiscope reaction and
        the list of species as arguments.
    
    Returns:
    --------
    object
        The corrected Kinetiscope reaction if duplicates are found; otherwise, 
        the original reaction is returned unchanged.
    """

    if duplicate_test(species_list):
        
        return correction_function(kinetiscope_reaction, species_list)
    
    return kinetiscope_reaction
    
def correct_kinetiscope_name(kinetiscope_reaction, product_list, reactant_list):
    """
    Correct the Kinetiscope reaction name based on duplicates in reactants and 
    products.
    
    This function checks for duplicate species in both the reactants and
    products of a Kinetiscope reaction. If duplicates are found, it applies the
    appropriate correction to the reaction name.
    
    Parameters:
    -----------
    kinetiscope_reaction : object
        An object representing a Kinetiscope reaction, which should have 
        methods or attributes 
        that allow for correcting the reaction name.
    product_list : list
        A list of products to be checked for duplicates.
    reactant_list : list
        A list of reactants to be checked for duplicates.
    
    Returns:
    --------
    object
        The corrected Kinetiscope reaction object if duplicates are found; 
        otherwise, the original reaction is returned unchanged.
    """

    kinetiscope_reaction = (
        correct_if_duplicates(kinetiscope_reaction, reactant_list, correct_name_reactants)
    )
    
    kinetiscope_reaction = (
        correct_if_duplicates(kinetiscope_reaction, product_list, correct_name_products)
    )
    
    return kinetiscope_reaction

def are_species_duplicates(species_list):
    """
    Determine if there are duplicate species in a list.
    
    This function checks if the list of species contains exactly two elements 
    and if both elements are identical, indicating a duplicate.
    
    Parameters:
    -----------
    species_list : list
        A list of species to be checked for duplicates. The list should contain
        exactly two elements for this function to determine if they are 
        duplicates. Note that sorting shouldn't matter for duplicates here.
    
    Returns:
    --------
    bool
        `True` if both elements in the list are identical (i.e., duplicates), 
        otherwise `False`.
    """
    
    return species_list[0] == species_list[1]

def duplicate_test(species_list):
    """
    Test for duplicate species in a list.
    
    This function checks if the list of species contains exactly two elements 
    and if those two elements are duplicates. It combines the checks for the 
    length of the list and the duplication of its elements.
    
    Parameters:
    -----------
    species_list : list
        A list of species to be tested for duplicates. The list should contain 
        exactly two elements for the function to determine if they are 
        duplicates.
    
    Returns:
    --------
    bool
        `True` if the list contains exactly two elements and they are 
        duplicates, otherwise `False`.
    """
    
    return len(species_list) == 2 and are_species_duplicates(species_list)

def validate_and_correct_reaction_name(kinetiscope_reaction):
    """
    Check for duplicate species in a Kinetiscope reaction and correct the 
    reaction name if needed.
    
    This function analyzes a Kinetiscope reaction to determine if there are 
    duplicate species among the reactants or products. If duplicates are found,
    it corrects the Kinetiscope name accordingly. If no duplicates are found, 
    the original Kinetiscope reaction is returned unchanged.
    
    Parameters:
    -----------
    kinetiscope_reaction : object
        An object representing a Kinetiscope reaction, which should have 
        methods or attributes that allow for extracting reactants and products.
    
    Returns:
    --------
    tuple
        A tuple containing:
        - A boolean indicating whether duplicates were found in either the 
        reactants or products.
        - The corrected Kinetiscope reaction name if duplicates were found; 
        otherwise, the original reaction is returned unchanged.
    """

    reactant_list, product_list = (
        find_reactants_and_products(kinetiscope_reaction)
    )
    
    reactants_are_duplicates = duplicate_test(reactant_list)
    products_are_duplicates = duplicate_test(product_list)
    
    #True if either reactants or products have duplicates, False otherwise
    
    reaction_has_duplicates = (
        reactants_are_duplicates or products_are_duplicates
    )
    
    reaction = (
        correct_kinetiscope_name(kinetiscope_reaction, product_list, reactant_list)
    )
    
    return reaction_has_duplicates, reaction
