# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 14:24:24 2024

@author: jacob
"""

from Rxn_classes import Kinetiscope_Reaction

class ReactionDataStorage:
    def __init__(self, name_mpculeid_dict,marker_species_dict,excitation_set,absorption_dict=None):

        self.mpculeid_dict = self.build_mpculeid_dict(name_mpculeid_dict)
        self.marker_species_dict = marker_species_dict
        self.excitation_set = excitation_set
        self.rate_constant_dict = self.create_rate_constant_dict(absorption_dict)
    
    def build_mpculeid_dict(self, name_mpculeid_dict):
        return {mpculeid: name for name, mpculeid in name_mpculeid_dict.items()}
    
    def create_rate_constant_dict(self, absorption_dict=None):
        """
        Creates and returns a nested dictionary to hold rate constants for various reaction types.

        The dictionary is organized into two main categories: 'ionization' and 'chemical'.
        Each category contains specific subcategories for different reaction processes.

        Parameters
        ----------
        absorption_dict : dict, optional
            A dictionary mapping molecular IDs to absorption rate constants. 
            If provided, these values will be added to the 'absorption' subcategory.

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
            Adds a dictionary of absorption rate constants to the rate_constant_dict.

            Parameters
            ----------
            rate_constant_dict : dict
                The rate constant dictionary created by create_rate_constant_dict.
            absorption_dict : dict
                A dictionary mapping molecular IDs to absorption rate constants.
                Example: {"4864aee73a83d357c31fadd50b81e3cd-C10H20O2-0-1": 1.4e-01}

            """
            if absorption_dict:
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

def build_reaction_parameter_dict(rate_constant_key, superclass, reaction_order, reactants_to_add=None, products_to_add=None):
    reaction_parameter_dict = {
        "rate_constant_key":rate_constant_key,
        "superclass":superclass,
        "reaction_order":reaction_order,
        "reactants_to_add":reactants_to_add,
        "products_to_add":products_to_add}
    return reaction_parameter_dict
        
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

