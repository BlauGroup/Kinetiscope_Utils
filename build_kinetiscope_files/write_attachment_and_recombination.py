# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 11:16:15 2024

@author: jacob
"""

from kinetiscope_reaction_writing_utilities import build_ionization_reaction
import sys
sys.path.append('../classify_HiPRGen_reactions')
from reaction_classification_utilities import find_mpculeid_charge

def test_for_positive_attachment(HiPRGen_reaction):
    """
    Helper function to determine whether a partiuclar attachment reaction,
    the only type it is called on in set_ionization_type_and_rate_constant, 
    is electron attachment to a cation

    Parameters
    ----------
    HiPRGen_reaction : HiPRGen reaction object
        the reaction we're building kinetiscope reaction objects from, defined 
        in Rxn_classes

    Returns
    -------
    Bool
        True if reaction is electron attachment to a cation, False otherwise

    """
    
    reactant_mpculeid = HiPRGen_reaction.reactants[0]
    reactant_charge = find_mpculeid_charge(reactant_mpculeid)
    
    return  reactant_charge > 0 

def set_ionization_type_and_rate_constant_key(HiPRGen_reaction):
    """
    Electron attachment to positive ions should have the same rate constant as
    recombination reactions, owing to the attraction between the negative
    electrons and positive ions. This function sets the ionization type and
    determines the appropriate key to find the rate constant for the reaction
    it's called for.

    Parameters
    ----------
    HiPRGen_reaction : HiPRGen reaction object
        the reaction we're building kinetiscope reaction objects from, defined 
        in Rxn_classes

    Returns
    -------
    ionization_type : str
        string describing the ionization type
    rate_constant_key : str
        string denoting the key we'll look for in the rate_constant_dict to
        find the rate constant for this particular reaction

    """
    if "recombination" in HiPRGen_reaction.tag:
        
        ionization_type, rate_constant_key = "recombination", "all"
        
    else:
        
        ionization_type = "attachment"
        attachment_is_positive = test_for_positive_attachment(HiPRGen_reaction)
        
        rate_constant_key = "positive" if attachment_is_positive else "neutral"
    
    return ionization_type, rate_constant_key

def build_recombination_or_attachment_dict(HiPRGen_reaction):
    """
    For building recombination or attachment reactions, this function gathers
    the necessary parameters into a dictionary for clarity and ease of use.
    
    Parameters
    ----------
    HiPRGen_reaction : HiPRGen reaction object
        The reaction we're building Kinetiscope reaction objects from, defined
        in Rxn_classes.
    
    Returns
    -------
    reaction_dict : dict
        A dictionary containing information used to construct the Kinetiscope_reaction
        object, including the ionization type, rate constant key, reaction order,
        reactants, and products.
    """
    
    ionization_type, rate_constant_key = (
        set_ionization_type_and_rate_constant_key(HiPRGen_reaction) 
    )
    
    product_list = HiPRGen_reaction.classification_list[:-1]
    reactant_list = ["TE"]
    
    reaction_dict = { 
        "rate_constant_key":rate_constant_key,
        "superclass":ionization_type,
        "reaction_order":2,
        "reactants_to_add":reactant_list,
        "products_to_add":product_list,
    }
    
    return reaction_dict

def build_attachment_or_recombination(HiPRGen_reaction, reaction_writing_data):
    """
    Electron attachment and recombination reactions are pretty similar in how
    they're written, so this function just lets us construct both after
    classifying the reaction we're currently building.

    Parameters
    ----------
    HiPRGen_reaction : HiPRGen reaction object
        the reaction we're building kinetiscope reaction objects from, defined 
        in Rxn_classes
    reaction_writing_data: ReactionDataStorage obj
        contains info related to a given kinetiscope reaction. Defined in 
        kinetiscope_reaction_writing_utilities
    Returns
    -------
    list
        the constructed Kinetiscope_reaction obj, contained as the only element
        of a list

    """
    
    reaction_dict = build_recombination_or_attachment_dict(HiPRGen_reaction)
    
    reaction = (
        build_ionization_reaction(HiPRGen_reaction, reaction_writing_data,reaction_dict)
    )
    
    return [reaction]