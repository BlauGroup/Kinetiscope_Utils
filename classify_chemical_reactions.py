 # -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:54:29 2024

@author: JRMilton
"""
import re
from reaction_classification_utilities import (
find_mpculeid_charge,
find_mpculeid_formula,
find_mpculeid_spin
)
from Rxn_classes import HiPRGen_Reaction
import sys
import copy

def narrow_H_rxn_type(reactant_gaining_H, product_with_H):
    """
    This function is called when we have determined that a given reaction
    involves the transfer of a hydrogen. It narrows down whether the reaction
    is proton transfer, hydride abstraction, H_atom abstraction, or 
    proton_coupled_electron_transfer based on charges or spins. 

    Parameters
    ----------
    reactant_gaining_H : string
       mpculeid of the form "graph_hash-formula-charge-spin"
    product_with_H : string
        mpculeid of the form "graph_hash-formula-charge-spin"

    Returns
    -------
    str
        narrowed description of this hydrogen-transfer reaction

    """
    
    reactant_charge = find_mpculeid_charge(reactant_gaining_H)
    product_charge = find_mpculeid_charge(product_with_H)
    reactant_gaining_H_spin = find_mpculeid_spin(reactant_gaining_H)
    product_with_H_spin = find_mpculeid_spin(product_with_H)
    
    delta_charge = product_charge - reactant_charge
    
    if delta_charge == 1:
        
        return "proton_transfer"
    
    elif delta_charge == -1:
        
        return "hydride_abstraction"
    
    is_radical_reaction = reactant_gaining_H_spin != product_with_H_spin
    
    if delta_charge == 0 and is_radical_reaction:
        
        return "H_atom_abstraction"
    
    else:
        
        return "proton_coupled_electron_transfer"

def count_number_hydrogens(formula):
    """
    Counts the number of hydrogen atoms in a chemical formula.

    Parameters
    ----------
    formula : str
        The chemical formula string.

    Returns:
    -------
    int
        The number of hydrogen atoms in the formula.
    """
    
    match = re.search(r'H(\d+)', formula)
    
    if match:
        hydrogen_count = int(match.group(1))
        
    else:
        hydrogen_count = 0
        
    return hydrogen_count

def add_hydrogen_to_formula(formula):
    """
    Some formulas we need to modify don't contain any hydrogens. Thus, we need
    to add hydrogen to these formulas, but because pymatgen generates formulas
    in alphabetical order, we need to insert hydrogen where it belongs for
    comparison. This function does that.

    Parameters
    ----------
    formula : str
        String of the form atom1#atom2# etc (e.g., 'C2H6O').

    Returns
    -------
    modified_formula : str
        formula with one hydrogen added

    """
    
    elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula) #list of tuples
                                                          #(atom, #)
    elements.append(('H', '1'))
    elements.sort(key=lambda x: x[0]) #sorts list by first element of tuples
    modified_formula = ''.join(f"{el}{count}" for el, count in elements)
    
    return modified_formula
    
def add_one_more_hydrogen(formula):
    """
    Increments the number of hydrogens in the given formula by 1.
    If there are no hydrogens in the formula, adds 'H1' in alphabetical order.

    Parameters
    ----------
    formula : str
        String of the form atom1#atom2# etc (e.g., 'C2H6O').

    Returns
    -------
    modified_formula : str
        The formula with the number of hydrogens increased by one.
    """
    
    number_hydrogens = count_number_hydrogens(formula)
    
    if number_hydrogens == 0:
        
        modified_formula = add_hydrogen_to_formula(formula)
        
    else:
        
        modified_hydrogen_count = number_hydrogens + 1
        modified_formula = re.sub(r'H\d*', f'H{modified_hydrogen_count}', formula)
    
    return modified_formula

def classify_H(rxn):
    """
    Determines whether a given reaction involves the transfer of a hydrogen atom
    by comparing the reactants of the reaction to the products. If a product has 
    one more hydrogen atom than a reactant, the reaction is an H atom transfer.
    Then, we narrow down the type via narrow_H_rxn_type and return the specific type.
    
    Parameters
    ----------
    rxn : HiPRGen rxn object
        The reaction we're testing.

    Returns
    -------
    str or None
        Narrowed H reaction type if the reaction is an H transfer, None otherwise.
    """
    
    for reactant_mpculeid in rxn.reactants:
        
        reactant_formula = find_mpculeid_formula(reactant_mpculeid)
        reactant_with_one_more_hydrogen = add_one_more_hydrogen(reactant_formula)
        
        for product_mpculeid in rxn.products:
            
            product_formula = find_mpculeid_formula(product_mpculeid)
            
            if product_formula == reactant_with_one_more_hydrogen:
                return narrow_H_rxn_type(reactant_mpculeid, product_mpculeid)
                
    return ""

def generate_formula_charge_dict(mpculeid):
    formula = find_mpculeid_formula(mpculeid)
    charge = find_mpculeid_charge(mpculeid)
    return {"formula":formula, "charge":charge}

def generate_mpculeid_formula_charge_dict(rxn):
    mpculeid_formula_charge_dict = {"reactants":{}, "products":{}}
    
    for reactant_mpculeid in rxn.reactants:
        mpculeid_formula_charge_dict["reactants"][reactant_mpculeid] = \
            generate_formula_charge_dict(reactant_mpculeid)
            
    for product_mpculeid in rxn.products:
        mpculeid_formula_charge_dict["products"][product_mpculeid] = \
            generate_formula_charge_dict(product_mpculeid)
    
    return mpculeid_formula_charge_dict

def get_sorted_formulas(mpculeid_formula_charge_dict, category):
    """
    Extracts and sorts formulas from the specified category (reactants or products)
    of the mpculeid_formula_charge_dict.
    """
    formulas = [entry["formula"] for entry in mpculeid_formula_charge_dict[category].values()]
    formulas.sort()
    return formulas

def charges_differ_by_one(reactant_charge, product_charge):
    """
    Checks if the charges of the reactant and product differ by one.
    """
    return abs(reactant_charge - product_charge) == 1

def is_matching_formula(reactant_formula, product_formula):
    """
    Checks if the reactant and product formulas match.
    """
    return reactant_formula == product_formula

def find_matching_product(reactant, products):
    """
    Finds a product with a matching formula and a charge differing by one.
    """
    for product_mpculeid, product_info in products.items():
        if is_matching_formula(reactant["formula"], product_info["formula"]) and \
           charges_differ_by_one(reactant["charge"], product_info["charge"]):
            return product_mpculeid
    return None

def reaction_is_electron_transfer(rxn):
    """
    Electron transfer reactions are those where the charges of the reactants
    are opposites, and the formulas of the products do not change, but their
    charges do. 

    Parameters
    ----------
    rxn : HiPRGen reaction object
        The reaction we're testing to see if it is electron transfer

    Returns
    -------
    bool
        True if the reaction is electron transfer, False otherwise
    """
    mpculeid_formula_charge_dict = generate_mpculeid_formula_charge_dict(rxn)
    
    sorted_reactant_formulas = get_sorted_formulas(mpculeid_formula_charge_dict, "reactants")
    sorted_product_formulas = get_sorted_formulas(mpculeid_formula_charge_dict, "products")
    
    if sorted_reactant_formulas != sorted_product_formulas:
        return ""
    
    test_dict = copy.deepcopy(mpculeid_formula_charge_dict)
    num_matches = 0
    
    for reactant_mpculeid, reactant_info in test_dict["reactants"].items():
        matching_product_mpculeid = find_matching_product(reactant_info, test_dict["products"])
        if matching_product_mpculeid:
            num_matches += 1
            test_dict["products"].pop(matching_product_mpculeid)
    
    if num_matches == 2:
        return "electron_transfer"
    
    return ""

def determine_charge_name(mpculeid):
    charge = find_mpculeid_charge(mpculeid)
    
    if charge > 0:
        return "cation"
    elif charge < 0:
        return "anion"
    else:
        return "neutral"
    
def determine_species_subclass(mpculeid):
    charge_name = determine_charge_name(mpculeid)
    spin = find_mpculeid_spin(mpculeid)
    
    if spin == 2 and charge_name == "neutral":
        
        return charge_name + "_radical"
    
    elif spin == 2:
        
        return "radical_" + charge_name
    
    return charge_name
    
def handle_unimolecular_reactions(rxn):
    single_reactant_mpculeid = rxn.reactants[0]
    species_subclass = determine_species_subclass(single_reactant_mpculeid)
    
    if len(rxn.products) == 1:
        return species_subclass + "_isomerization"
    
    return species_subclass + "_fragmentation"

def generate_ordered_name(reactant_1_subclass, reactant_2_subclass):
    
    priority_dict = {
            "radical_anion":5,
            "radical_cation":4,
            "anion":3,
            "cation":2,
            "neutral_radical":1,
            "neutral":0
            }
    
    reactant_1_priority = priority_dict.get(reactant_1_subclass)
    reactant_2_priority = priority_dict.get(reactant_2_subclass)
    
    if reactant_1_priority > reactant_2_priority:
        ordered_name = reactant_1_subclass + "_" + reactant_2_subclass
    else:
        ordered_name = reactant_2_subclass + "_" + reactant_1_subclass
    
    return ordered_name  

def determine_bimolecular_reactant_subclasses(rxn):
    reactant_1 = rxn.reactants[0]
    reactant_2 = rxn.reactants[1]
    
    reactant_1_subclass = determine_species_subclass(reactant_1)
    reactant_2_subclass = determine_species_subclass(reactant_2)
    
    ordered_name = \
        generate_ordered_name(reactant_1_subclass, reactant_2_subclass)

    return ordered_name

def handle_bimolecular_reactions(rxn):
    bimolecular_subclass = determine_bimolecular_reactant_subclasses(rxn)
    
    if len(rxn.products) == 1:
        return bimolecular_subclass + "_combination"
    
    H_name = classify_H(rxn)
    
    if H_name:
        return bimolecular_subclass + "_" + H_name
    
    electron_transfer_name = reaction_is_electron_transfer(rxn)
    
    if electron_transfer_name:
        return bimolecular_subclass + "_" + electron_transfer_name
    
    return bimolecular_subclass + "_reaction"

def determine_chemical_reaction_tag(rxn):
    """
    Tags chemical reactions by calling classification functions in a specific order.

    Parameters
    ----------
    rxn : HiPRGen rxn object
        reaction we're classifying
        
    Returns
    -------
    list
        The tag for the classified reaction
    """
    if len(rxn.reactants) == 1:
        
        return handle_unimolecular_reactions(rxn)
    
    return handle_bimolecular_reactions(rxn)