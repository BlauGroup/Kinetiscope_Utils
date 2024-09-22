# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 14:24:24 2024

@author: jacob
"""

import re
import sys
sys.path.append('../common')
from Rxn_classes import Kinetiscope_Reaction

class ReactionDataStorage:
    def __init__(self, name_mpculeid_dict,marker_species_dict,excitation_set,absorption_dict=None):
        """
        A lot of the functions used to tag reactions uses convoluted data
        structures. This just collects them into a class so we can only call
        what we need in a given function.

        Parameters
        ----------
        name_mpculeid_dict : dict
            a dictionary with chemical names as keys and mpculeids as values
        marker_species_dict : dict
            a dictionary associating the long forms of names--e.g. radical
            cation--with their shorthand--e.g. rc
        excitation_set : set
            a set we use to check to see if an excitation reaction has already
            been created for a given species
        absorption_dict : dict, optional
            a dictionary with mpculeids as the keys pseudo-1st-order absorption
            rate constants as values. By default None.

        """

        self.mpculeid_dict = self.build_mpculeid_dict(name_mpculeid_dict)
        self.marker_species_dict = marker_species_dict
        self.excitation_set = excitation_set
        self.rate_constant_dict = self.create_rate_constant_dict(absorption_dict)
    
    def build_mpculeid_dict(self, name_mpculeid_dict):
        """
        This function just creates a new dictionary with mpculeids as the keys
        and names as the values

        Parameters
        ----------
        name_mpculeid_dict : dict
            a dictionary with chemical names as keys and mpculeids as values

        Returns
        -------
        dict
            a dictionary with mpculeids as keys and names as values

        """
        
        return {mpculeid: name for name, mpculeid in name_mpculeid_dict.items()}
    
    def create_rate_constant_dict(self, absorption_dict=None):
        """
        Creates and returns a nested dictionary to hold rate constants for 
        various reaction types.

        The dictionary is organized into two main categories: 'ionization' and 
        'chemical'. Each category contains specific subcategories for different
        reaction processes.

        Parameters
        ----------
        absorption_dict : dict, optional
            A dictionary mapping molecular IDs to absorption rate constants. 
            If provided, these values will be added to the 'absorption' 
            subcategory.

        Returns
        -------
        dict
            A dictionary with the following structure:
            - 'ionization':
                - 'absorption': {<molecular_id>: <rate_constant>, ...}
                - 'electron_ionization': {
                    'eV_80': 1.4e+15,
                    'eV_55': 1.2e+15,
                    'eV_30': 8.5e+14,
                    'LEE': 1.7e+14
                }
                - 'recombination': {'all': 4.1e+13}
                - 'attachment': {'positive': 4.1e+13, 'neutral': 3.1e+12}
            - 'chemical':
                - '1st_order': 6.2e+12
                - '2nd_order': 1.3e+12
                
        all values are the bottom of the dict are floats written in scientific
        notation. 
        
        """
        def add_absorption_values(rate_constant_dict, absorption_dict):
            """
            Adds a dictionary of absorption rate constants to the 
            rate_constant_dict.

            Parameters
            ----------
            rate_constant_dict : dict
                The rate constant dictionary created by 
                create_rate_constant_dict.
            absorption_dict : dict
                A dictionary mapping molecular IDs to absorption rate constants
                Example: {"4864aee73a83d357c31fadd50b81e3cd-C10H20O2-0-1": 1.4e-01}

            """
            rate_constant_dict["ionization"]["absorption"].update(absorption_dict)
            
        rate_constant_dict = {
            "ionization": {
                "absorption": {},
                "electron_ionization": {
                    "eV_80": 1.4e+15,
                    "eV_55": 1.2e+15,
                    "eV_30": 8.5e+14,
                    "LEE": 1.7e+14,
                },
                "recombination": {"all": 4.1e+13},
                "attachment": {"positive": 4.1e+13, "neutral": 3.1e+12},
            },
            "chemical": {
                "1st_order": 6.2e+12,
                "2nd_order": 1.3e+12,
            }
        }

        if absorption_dict:
            
            add_absorption_values(rate_constant_dict, absorption_dict)

        return rate_constant_dict

def determine_ordinal_number_order(HiPRGen_rxn):
    """
    Determine the order of a chemical reaction based on the number of reactants.
    
    This function assesses the number of reactants in a chemical reaction and 
    determines if the reaction is first-order or second-order. It returns 
    "1st_order" if there is exactly one reactant and "2nd_order" otherwise.
    
    Parameters:
    -----------
    HiPRGen_rxn : object
        An object representing a chemical reaction, which should have a 
        `reactants` attribute that is a list of reactants.
    
    Returns:
    --------
    str
        A string indicating the reaction order: "1st_order" if there is one 
        reactant, "2nd_order" otherwise.
    """
  
    if len(HiPRGen_rxn.reactants) == 1:
        return "1st_order"
    
    return "2nd_order"
        
def add_reactants(base_name, reactants_to_add):
    """
    Add reactants to a base name in a chemical reaction string.
    
    This function takes a base name and a list of reactants, and concatenates 
    them into a single string. The reactants are joined with the base name 
    using " + " as a separator, with the reactants appearing before the base 
    name, as is typical in chemical reaction notation.
    
    Parameters:
    -----------
    base_name : str
        The name or identifier of the chemical species to which reactants 
        will be added.
    reactants_to_add : list
        A list of reactants to be prepended to the base name.
    
    Returns:
    --------
    str
        A string representing the reactants followed by the base name, 
        separated by " + ".
    """

    return " + ".join(reactants_to_add + [base_name])

def add_products(base_name, added_products):
    """
    Add products to a base name in a chemical reaction string.
    
    This function takes a base name and a list of products, and concatenates 
    them into a single string. The base name is joined with the products using 
    " + " as a separator, which is typical in chemical reaction notation.
    
    Parameters:
    -----------
    base_name : str
        The starting name or identifier of the chemical species.
    added_products : list
        A list of products to be appended to the base name.
    
    Returns:
    --------
    str
        A string representing the base name followed by the added products, 
        separated by " + ".
    """

    return " + ".join([base_name] + added_products)
    
def add_reacting_species(initial_name, reactants_to_add, added_products):
    """
    Add reactants and products to a chemical species name.
    
    This function takes an initial name and potentially adds reactants and/or 
    products to it. If `reactants_to_add` is provided, they are added to the 
    `initial_name`. Then, regardless of whether reactants were added, products
    are added to the resulting name.
    
    Parameters:
    -----------
    initial_name : str
        The starting name of the chemical species.
    reactants_to_add : list
        A list of reactants to be added to the initial name. If the list is 
        empty, no reactants 
        will be added.
    added_products : list
        A list of products to be added to the name after reactants (if any) 
        have been added.
    
    Returns:
    --------
    str
        The name after potentially adding both reactants and products.
    """

    current_name = initial_name

    if reactants_to_add:
        
        current_name = add_reactants(initial_name, reactants_to_add)
    
    return add_products(current_name, added_products)

def replace_mpculeids_with_names(HiPRGen_name, mpculeid_dict):
    """
    Replace molecular IDs in a reaction string with their corresponding names.
    
    This function takes a reaction string generated by HiPRGen, which may 
    include molecular IDs (mpculeids) and other symbols like "+" or "=>". It 
    replaces each mpculeid with the corresponding Kinetiscope name from a 
    provided dictionary, while leaving other symbols unchanged.
    
    Parameters:
    -----------
    HiPRGen_name : str
        A string representing a reaction, potentially including mpculeids and 
        symbols.
    mpculeid_dict : dict
        A dictionary mapping mpculeids (keys) to their corresponding
        Kinetiscope names (values).
    
    Returns:
    --------
    str
        The reaction string with mpculeids replaced by their corresponding 
        Kinetiscope names.
    """
    
    #this will include mpculeids, which we will replacy, as well as nonsense
    #like "+" or "=>", which we will not change
    
    original_reaction_list = HiPRGen_name.split()
    reaction_list_copy = original_reaction_list[:]
    
    for index, string in enumerate(reaction_list_copy):
        
        kinetiscope_name = mpculeid_dict.get(string, None)
        string_is_mpculeid = kinetiscope_name
        
        if string_is_mpculeid:
            
            original_reaction_list[index] = kinetiscope_name
    
    return " ".join(original_reaction_list)
    
def write_kinetiscope_name(HiPRGen_name, mpculeid_dict, added_reactants, added_products):
    """
    Takes the name of a chemical reaction from HiPRGen and converts it to a
    more human-readable, and more importantly, kinetiscope-readable form. Also
    adds marker species to that reaction.

    Parameters
    ----------
    HiPRGen_name : str
        The name written for this particular reaction in HiPRGen, with
        mpculeids representing chemical species
    mpculeid_dict : dict
        a dictionary, described in the ReactionDataStorage class, whose keys
        are mpculeids and whose values are chemical names generated via
        name_molecules
    added_reactants : list
        a list of reactants to added to the name of the reaction. May represent
        marker species or real chemical species, or may also be none.
    added_products : list
        a list of products to added to the name of the reaction. May represent
        marker species or real chemical species, or may also be none.

    Returns
    -------
    kinetiscope_name : str
        contains one or more chemical reactants or products, whose names are
        generated via name_molecules, and connected with "=>"

    """
    
    reaction_with_names = (
        replace_mpculeids_with_names(HiPRGen_name, mpculeid_dict)
    )
    
    kinetiscope_name = (
        add_reacting_species(reaction_with_names, added_reactants, added_products)
    )
    
    return kinetiscope_name 

def build_rxn_object(HiPRGen_rxn, kinetiscope_name, rate_constant, order, marker_species):
    """
    A simple function just for taking in attributes and building a
    Kinetiscope_reaction object, as defined in Rxn_classes, from them.

    Parameters
    ----------
    HiPRGen_rxn : HiPRGen_rxn obj
        data related to a reaction from HiPRGen, outlined in Rxn_classes
    kinetiscope_name : str
        contains one or more chemical reactants or products, whose names are
        generated via name_molecules, and connected with "=>"
    rate_constant : float
        written in scientific notation
    order : int
        order of the chemical reaction
    marker_species : list
        a list of marker species associated with a given reaction

    Returns
    -------
    Kinetiscope_reaction object
        the Kinetiscope_reaction object associated with this chemical reaction

    """
    
    return Kinetiscope_Reaction(HiPRGen_rxn, kinetiscope_name, rate_constant, order, marker_species)

def replace_tag_with_shorthand(marker_species, shorthand_dict):
    """
    Kinetiscope's reaction names can only have a maximum of 32 characters. This
    function replaces some long names (e.g. radical cation) with shorthand to
    make sure each individual species has a name with less than 32 characters.

    Parameters
    ----------
    marker_species : list
        the list of current marker species
    shorthand_dict : dict
        a dict where each key is the longhand form of a name we're trying to
        replace and whose values are shorthand

    Returns
    -------
    marker_species : list
        the list of marker species, with updated shortened names

    """
    
    #tag is a very specific tag for a chemical reaction, generated in
    #gen_HiPRGen_rxn_dict_direct
    
    tag = marker_species[-1]
    
    for longhand, shorthand in shorthand_dict.items():
        if longhand in tag:
            tag = tag.replace(longhand, shorthand)
            
    marker_species[-1] = tag
    
    return marker_species

def find_ionization_rate_constant(reaction_writing_data, reaction_dict, ionization_type):
    """
    Each Kineticope_Reaction object needs a rate constant assigned to it for
    our simulations. This function finds the rate constant associated with this
    particular ionization reaction and returns it for assignment.

    Parameters
    ----------
    reaction_writing_data: ReactionDataStorage object
        a class that stores data related to chemical reactions necessary to
        write their kinetiscope names, defined in 
        kinetiscope_reaction_writing_utilities
    reaction_dict : dict
        a dictionary, built using the build_reaction_parameter_dict function
        in kinetiscope_reaction_writing_utilities, which contains information
        associated with this particular kinetiscope reaction
    ionization_type : string
        A string, describing the ionization superclass this reaction belongs to

    Returns
    -------
    float
        a floating point, written in scientific notation, representing the rate
        constant associated with this particular reaction. Units for first-order
        rate constants are s^-1 while for second-order they are M^-1 s^-1.

    """
    
    rate_constant_key = reaction_dict["rate_constant_key"]
    all_rate_constants = reaction_writing_data.rate_constant_dict
    ionization_rate_constants = all_rate_constants["ionization"][ionization_type]
        
    return ionization_rate_constants.get(rate_constant_key, None)

def write_ionization_name(HiPRGen_reaction, reaction_writing_data, reaction_dict, ionization_type, products_to_add):
    """
    Generates a name for an ionization reaction to be assigned to a 
    Kinetiscope_Reaction object.

    Parameters
    ----------
    HiPRGen_reaction : HiPRGen reaction object
        the reaction we're building kinetiscope reaction objects from, defined 
        in Rxn_classes
    reaction_writing_data: ReactionDataStorage object
        a class that stores data related to chemical reactions necessary to
        write their kinetiscope names, defined in 
        kinetiscope_reaction_writing_utilities
    reaction_dict : dict
        a dictionary, built using the build_reaction_parameter_dict function
        in kinetiscope_reaction_writing_utilities, which contains information
        associated with this particular kinetiscope reaction
    ionization_type : string
        A string, describing the ionization superclass this reaction belongs to,
        to be added as a product of the ionization reaction
    products_to_add : list
        a list of products of this reaction to which ionization_type will be
        appended

    Returns
    -------
    string
        a name for this kinetiscope reaction, written using 
        write_kinetiscope_name defined in 
        kinetiscope_reaction_writing_utilities

    """
    
    mpculeid_dict = reaction_writing_data.mpculeid_dict
    reactants_to_add = reaction_dict["reactants_to_add"]    
    HiPRGen_name = HiPRGen_reaction.name
    
    return write_kinetiscope_name(HiPRGen_name, mpculeid_dict, reactants_to_add, products_to_add)
   
def build_ionization_reaction(HiPRGen_reaction, reaction_writing_data, reaction_dict):
    """
    Ionization reactions can have reactants or products to add to the 
    Kinetiscope name depending on the reaction type. This function generates
    a name for this ionization reaction, finds the appropriate rate constant
    for the reaction, and finally builds and returns a Kinetiscope_Reaction
    object associated with this ionization reaction.

    Parameters
    ----------
    HiPRGen_reaction : HiPRGen reaction object
        the reaction we're building kinetiscope reaction objects from, defined 
        in Rxn_classes
    reaction_writing_data: ReactionDataStorage object
        a class that stores data related to chemical reactions necessary to
        write their kinetiscope names, defined in 
        kinetiscope_reaction_writing_utilities
    reaction_dict : dict
        a dictionary, which contains information associated with this
        particular kinetiscope reaction. This is build before this function is
        called so that the parameters passed aren't unwieldy

    Returns
    -------
    ionization_reaction: Kinetiscope_Reaction object
        the newly constructed Kinetiscope_Reaction object associated with this 
        ionization reaction

    """
    
    ionization_type = reaction_dict["superclass"]
    products_to_add = reaction_dict["products_to_add"]
    products_to_add.append(ionization_type)
    reaction_order = reaction_dict["reaction_order"]
    
    kinetiscope_name = (
        write_ionization_name(HiPRGen_reaction, reaction_writing_data, reaction_dict, ionization_type, products_to_add)
    )
    
    rate_constant = (
        find_ionization_rate_constant(reaction_writing_data, reaction_dict, ionization_type)
    )
    
    ionization_reaction = (
    build_rxn_object(HiPRGen_reaction, kinetiscope_name, rate_constant, reaction_order, products_to_add)
    )
    
    return ionization_reaction

def reclassify_reaction(HiPRGen_reaction):
    """
    Reclassifies the reaction by updating its classification tag.
    
    Parameters:
    - HiPRGen_reaction (object): The reaction object whose classification will
    be updated.
    
    Process:
    1. Retrieves the current classification tag.
    2. Replaces "combination" with "crosslinking" in the tag.
    3. Updates the reaction's classification list and tag.
    
    Returns:
    - None
    """

    current_classifications = HiPRGen_reaction.classification_list
    current_tag = current_classifications[-1]
    new_tag = current_tag.replace("combination", "crosslinking")
    current_classifications[-1] = new_tag
    HiPRGen_reaction.tag = new_tag
    
def reaction_is_crosslinking(product_name):
    """
    Checks if the input string indicates a crosslinking reaction by verifying:
    - The presence of both "PHS" and "pMMA" with any number following them.
    - OR the presence of "PHS" or "pMMA" with a number greater than or equal to
    2.

    Parameters:
    - product_name (str): The string to check for crosslinking indicators.

    Returns:
    - bool: True if the string meets either of the conditions, False otherwise.
    """

    # Pattern for both "PHS" and "pMMA" with any number
    both_backbones_pattern = r"PHSb\d*.*PtBMAb\d*|PtBMAb\d*.*PHSb\d*"

    # Pattern for "PHS" or "pMMA" with number >= 2
    two_same_backbone_pattern = r"(PHSb[2-9]\d*|PtBMAb[2-9]\d*)"

    contains_both_backbones = re.search(both_backbones_pattern, product_name)

    contains_two_same = re.search(two_same_backbone_pattern, product_name)

    return contains_both_backbones or contains_two_same

def reclassify_if_crosslinking(HiPRGen_reaction, reaction_writing_data):
    """
    Reclassifies the reaction if it is identified as crosslinking.
    
    Parameters:
    - HiPRGen_reaction (object): The reaction object with product details.
    - reaction_writing_data (object): Contains the dictionary mapping product
    IDs to names.
    
    Process:
    1. Retrieves the product name from the dictionary.
    2. Raises a KeyError if the product ID is missing.
    3. Reclassifies the reaction if it is crosslinking.
    
    Raises:
    - KeyError: If the mpculeid ID is not found in the dictionary.
    
    Returns:
    - None
    """
    
    product_mpculeid =  HiPRGen_reaction.products[0] #combination reactions 
                                                     #have only one product
    
    product_name = (
        reaction_writing_data.mpculeid_dict.get(product_mpculeid, None)
    )
    
    if not product_name:
        
        raise KeyError(f"mpculeid {product_mpculeid} not in dictionary")
        
    if reaction_is_crosslinking(product_name):

        reclassify_reaction(HiPRGen_reaction)

def reclassify_crosslinking_reactions(
        HiPRGen_reaction_list, reaction_writing_data
    ):
    """
    Reclassifies reactions in a list if they are identified as crosslinking.
    
    Parameters:
    - HiPRGen_reaction_list (list): List of reaction objects to be checked and 
    potentially reclassified.
    - reaction_writing_data (object): Contains data for checking crosslinking, 
    including a dictionary of molecule IDs.
    
    Process:
    1. Iterates through each reaction in the list.
    2. Checks if "combination" is in the reaction's classifications.
    3. Reclassifies the reaction if it is crosslinking.
    
    Returns:
    - None
    """

    for index, HiPRGen_reaction in enumerate(HiPRGen_reaction_list):

        current_classifications = HiPRGen_reaction.classification_list

        if "combination" in current_classifications:

            reclassify_if_crosslinking(HiPRGen_reaction, reaction_writing_data)
    
def check_species_lengths(kinetiscope_reaction):
    """
    Processes the input string to ensure that each element of the resulting
    list is less than or equal to 32 characters long. Raises an error if 
    anyelement exceeds this length.
    
    Parameters:
    - input_string (str): The string to be processed. Whitespace is removed,
    and the string is split into a list of non-whitespace elements.
    
    Raises:
    - ValueError: If any element in the list is longer than 32 characters.
    
    Returns:
    - None
    """
    
    # Remove any extra whitespace and split the string into a list of words
    species_list = kinetiscope_reaction.split()
    
    # Check the length of each word and raise an error if it's too long
    for species_name in species_list:
        if len(species_name) > 32:
            raise ValueError(f"{species_name} has {len(species_name)} chars.")

def check_reaction_count(reaction_string):
    """
    Checks if the number of reactants and products in the reaction string is
    less than 9. Raises an error if there are more than 8 reactants or products.
    
    Parameters:
    - reaction_string (str): The reaction string to be checked. It should be 
    in the format 'reactant1 + reactant2 + ... => product1 + product2 + ...'.
    
    Raises:
    - ValueError: If there are more than 8 reactants or products.
    
    Returns:
    - None
    """
    
    try:
        
        reactants, products = reaction_string.split(" => ")
        
    except ValueError:
        
        raise ValueError("Reaction string must contain exactly one ' => ' separator.")
    
    # Split reactants and products into lists
    reactants_list = reactants.split(" + ")
    products_list = products.split(" + ")
    
    # Check if the number of reactants or products exceeds 8
    if len(reactants_list) > 8:
        raise ValueError("Number of reactants exceeds 8.")
    if len(products_list) > 8:
        raise ValueError("Number of products exceeds 8.")