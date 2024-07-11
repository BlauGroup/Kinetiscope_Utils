# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 09:33:22 2024

@author: jacob
"""

def find_mpculeid_charge(mpculeid):
    """
    This function simply returns the charge of a given species as an interger
    based on its mpculeid.

    Parameters
    ----------
    mpculeid : string
        a string of the form: "graph_hash-formula-charge-spin" 

    Returns
    -------
    charge : int
        charge of the species
    """
    
    charge_str = mpculeid.split("-")[2]
    if "m" in charge_str: #m stands for minus in the string
        charge = -int(charge_str.replace("m", ""))
        
    else:
        charge = int(charge_str)
    return charge

def find_reactant_and_product_charges(reaction):
    reactant_mpculeid = reaction.reactants[0]
    product_mpculeid = reaction.products[0]
    reactant_charge = find_mpculeid_charge(reactant_mpculeid)
    product_charge = find_mpculeid_charge(product_mpculeid)
    
    return reactant_charge, product_charge

def add_reaction_if_new(new_rxn, tag, rxns_for_simulation, rxns_already_added):
    
    if new_rxn.reaction_hash not in rxns_already_added:
        
        new_rxn.tag = tag
        rxns_for_simulation.append(new_rxn)
        rxns_already_added.add(new_rxn.reaction_hash)
        
    return rxns_for_simulation, rxns_already_added