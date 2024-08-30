# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 09:33:22 2024

@author: jacob
"""

def find_mpculeid_formula(mpculeid):
    """
    This function returns the formula of an mpculeid.

    Parameters
    ----------
    mpculeid : string
        a string of the form: "graph_hash-formula-charge-spin" 

    Returns
    -------
    string
        chemical formula of the mpculeid
    """
    
    return mpculeid.split("-")[1]

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

def find_mpculeid_spin(mpculeid):
    """
    This function returns the spin of an mpculeid.

    Parameters
    ----------
    mpculeid : string
        a string of the form: "graph_hash-formula-charge-spin" 

    Returns
    -------
    int
        spin associated with this mpculeid
    """
    
    return int(mpculeid.split("-")[3])

def find_reactant_and_product_charges(reaction):
    """
    For ionization reactions, we only have one product and one reactant, so
    this function returns the charge of the reactant and product.

    Parameters
    ----------
    reaction : HiPRGen reaction object
        the reaction we're calculating the charges for

    Returns
    -------
    reactant_charge : int
        charge of the reactant
    product_charge : int
        charge of the product

    """
    
    reactant_mpculeid = reaction.reactants[0]
    product_mpculeid = reaction.products[0]
    reactant_charge = find_mpculeid_charge(reactant_mpculeid)
    product_charge = find_mpculeid_charge(product_mpculeid)
    
    return reactant_charge, product_charge

def add_reaction_if_new(new_rxn, tag, rxns_for_simulation, rxns_already_added):
    """
    Some reactions are in both P1 and P2, and some are already pulled from the
    P2 reaction tally before we look for pathways forming network products.
    This function ensures that a new reaction has not already been added to 
    our list.

    Parameters
    ----------
    new_rxn : HiPRGen reaction object
        the potentially new reaction we're testing
    tag : string
        a tag describing what kind of reaction this is
    rxns_for_simulation : list
        list of reactions we'll be simulating
    rxns_already_added : set
        set of reaction names for reactions already added

    Returns
    -------
    rxns_for_simulation : list
        modified if new_rxn was added, unmodified otherwise
    rxns_already_added : set
        modified if new_rxn was added, unmodified otherwise

    """
    
    if new_rxn.name not in rxns_already_added:
        
        new_rxn.tag = tag
        rxns_for_simulation.append(new_rxn)
        rxns_already_added.add(new_rxn.name)
        
    return rxns_for_simulation, rxns_already_added

def find_reactant_charges(reaction):
    reactant_1_charge = find_mpculeid_charge(reaction.reactants[0])
    reactant_2_charge = find_mpculeid_charge(reaction.reactants[1])
    return reactant_1_charge, reactant_2_charge

def reaction_is_neutral(reaction):
    reactant_1_charge, reactant_2_charge = find_reactant_charges(reaction)
    return reactant_1_charge == reactant_2_charge == 0

def reaction_is_ion_molecule(reaction):
    reactant_1_charge, reactant_2_charge = find_reactant_charges(reaction)
    charges_are_different = reactant_1_charge != reactant_2_charge
    one_charge_zero = reactant_1_charge == 0 or reactant_2_charge == 0
    return charges_are_different and one_charge_zero

def determine_bimolecular_subtype(reaction):
    if reaction_is_neutral(reaction):
        return "neutral"
    elif reaction_is_ion_molecule(reaction):
        return "ion-molecule"
    else:
        return "ion-ion"

def write_bimolecular_classification(reaction):
    subtype = determine_bimolecular_subtype(reaction)
    if "combination" in reaction.tag:
        subclass = "combination"
        return ["chemical", "bimolecular", subclass, subtype, reaction.tag]
    
    subclass = "biproduct"
    subcategories = [
        "proton_transfer",
        "hydride_abstraction",
        "H_atom_abstraction",
        "proton_coupled_electron_transfer",
        "electron_transfer",
        "reaction"
    ]
    
    for subcategory in subcategories:
        if subcategory in reaction.tag:
            return ["chemical", "bimolecular", subclass, subtype, subcategory, reaction.tag]

def write_unimolecular_classification(reaction):
    subclass = "fragmentation" if "fragmentation" in reaction.tag else "isomerization"
    return ["chemical", "unimolecular", subclass, reaction.tag]

def write_chemical_classification(reaction):
    number_of_reactants = len(reaction.reactants)
    if number_of_reactants == 1:
        return write_unimolecular_classification(reaction)
    return write_bimolecular_classification(reaction)

def write_ionization_classification(reaction):
    return ["ionization", reaction.tag]

def write_reaction_classification(reaction):
    ionization_classifications = set(["positive_ionization", "electron_attachment",
                                      "electron_cation_recombination"])
    if reaction.tag in ionization_classifications:
        return write_ionization_classification(reaction)
    return write_chemical_classification(reaction)