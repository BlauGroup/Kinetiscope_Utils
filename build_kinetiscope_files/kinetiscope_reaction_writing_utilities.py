# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 14:24:24 2024

@author: jacob
"""

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

# def build_reaction_parameter_dict(rate_constant_key, superclass, reaction_order, reactants_to_add=None, products_to_add=None):
#     reaction_parameter_dict = {
#         "rate_constant_key":rate_constant_key,
#         "superclass":superclass,
#         "reaction_order":reaction_order,
#         "reactants_to_add":reactants_to_add,
#         "products_to_add":products_to_add}
#     return reaction_parameter_dict
        
def add_reactants(base_name, added_reactants):
    return " + ".join(added_reactants + [base_name])

def add_products(base_name, added_products):
    return " + ".join([base_name] + added_products)
    
def add_reacting_species(initial_name, added_reactants, added_products):
    name_after_reactants_alias = initial_name

    if added_reactants:
        name_after_reactants_alias = add_reactants(initial_name, added_reactants)
    
    return add_products(name_after_reactants_alias, added_products)

def replace_mpculeids_with_names(H_name, mpculeid_dict):
    reaction_as_list = H_name.split()
    
    for index, obj in enumerate(reaction_as_list[:]):
        kinetiscope_name = mpculeid_dict.get(obj, None)
        if kinetiscope_name:
            reaction_as_list[index] = kinetiscope_name
    
    return " ".join(reaction_as_list)
    
def write_kinetiscope_name(H_name, mpculeid_dict, added_reactants, added_products):
    reaction_with_names = replace_mpculeids_with_names(H_name, mpculeid_dict)
    kinetiscope_name = \
        add_reacting_species(reaction_with_names, added_reactants, added_products)
    
    return kinetiscope_name 

def build_rxn_object(HiPRGen_rxn, kinetiscope_name, rate_constant, order, marker_species):
    return Kinetiscope_Reaction(HiPRGen_rxn, kinetiscope_name, rate_constant, order, marker_species)

def replace_tag_with_shorthand(marker_species, shorthand_dict):
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
    H_name = HiPRGen_reaction.name
    
    return write_kinetiscope_name(H_name, mpculeid_dict, reactants_to_add, products_to_add)
   
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