# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:54:29 2024

@author: JRMilton
"""

from reaction_classification_utilities import (
    find_reactant_and_product_charges)

def reaction_is_ionization(reaction):
    """
    "Phase 1" of our EUV exposure model contains ionization reactions that
    we want to tag later when we're building our Kinetiscope reaction objects.
    This function tags these reactions.
    

    Parameters
    ----------
    reaction : HiPRGen_Reaction object
        the reaction we're testing

    Returns
    -------
    bool
        True if reaction is ionization, False otherwise.
    string
        None if reaction IS NOT ionization, string determined by 
        determine_broad_ionization_type if rxn is ionization
    """
    
    num_reactants = len(reaction.reactants)
    num_products = len(reaction.products)
    
    if num_reactants+num_products == 2:
        
        reactant_charge, product_charge = find_reactant_and_product_charges(reaction)
        
        if reactant_charge != product_charge: 
            return True
        
    return False

def determine_broad_ionization_tag(reaction):
    """
    Calculates the change in charge over the course of an ionization reaction
    and determines if the ionization reaction is attachment, recomibnation,
    or positive ionization

    Parameters
    ----------
    reactant_charge : int
        charge of the reactant species
    product_charge : int
        charge of the product species

    Returns
    -------
    str
        a tag denoting, broadly, what type of ionziation reaction this is.
        We will distinguish between attachment and recombination later.
    """
    
    reactant_charge, product_charge = find_reactant_and_product_charges(reaction)
    delta_charge = product_charge-reactant_charge
    
    if delta_charge < 0 and reactant_charge > 0:
        
        return "attachment_or_recombination"
    
    elif delta_charge < 0:
        
        return "electron_attachment"
    
    else:
        
        return "positive_ionization"

# def process_ionization_reactions(new_rxn, rxns_for_simulation, rxns_already_added):
#     """
#     This function updates ionization_rxn_list and ionization_hashes by adding
#     info for new_rxn if deemed necessary.

#     Parameters
#     ----------
#     new_rxn : HiPRGen rxn object
#         reaction we're classifying
#     ionization_rxn_list : list
#         list containing all HiPRGen rxns we've already made
#     ionization_hashes : set
#         set of hashes for ionization reactions we've already created

#     Returns
#     -------
#     ionization_rxn_list : list
#         potentially updated list
#     ionization_hashes : hash
#         potentially updated set
#     """
    
#     if new_rxn.reaction_hash not in rxns_already_added:
        
#         new_rxn.tag = determine_broad_ionization_tag(new_rxn)  
#         rxns_for_simulation.append(new_rxn)
#         rxns_already_added.add(new_rxn.reaction_hash)
    
#     return rxns_for_simulation, rxns_already_added

def is_reversed_reaction(reaction, rxn):
    """
    If a given reaction is recombination, that means that it has a
    positive_ionization reaction that its reverse. This function tests that.

    Parameters
    ----------
    reaction : HiPRGen rxn object
        the potential recombination reaction
    rxn : HiPRGen rxn object
        a positive_ionization reaction that may be associated with a
        we're testing reaction against
    Returns
    -------
    bool
        True if the tested reactions are reverse of eachother, False otherwise.
    """
    
    return reaction.reactants == rxn.products and reaction.products == rxn.reactants

def narrow_down_ionization_type(reaction, rxns_for_simulation):
    """
    If in our list of reactions, there exists a positive_reaction that is the
    reverse of our tested reaction, it is a recombination reaction. Otherwise,
    the reaction is electron attachment. Thus this function distinguishes
    between those two possibilities.

    Parameters
    ----------
    reaction : HiPRGen rxn object
        reaction we're testing to narrow down its ionization type
    ionization_rxn_list : list
        list of all reaction objects we've generated

    Returns
    -------
    str
        "recombination" if reaction is recombination, "electron_attachment"
        otherwise.
    """
    
    if any(is_reversed_reaction(reaction, rxn) for rxn in rxns_for_simulation):
        
        return "electron_cation_recombination"
    
    else:
        
        return "electron_attachment"