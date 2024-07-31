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
                
    return None

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
        return False
    
    test_dict = copy.deepcopy(mpculeid_formula_charge_dict)
    num_matches = 0
    
    for reactant_mpculeid, reactant_info in test_dict["reactants"].items():
        matching_product_mpculeid = find_matching_product(reactant_info, test_dict["products"])
        if matching_product_mpculeid:
            num_matches += 1
            test_dict["products"].pop(matching_product_mpculeid)
    
    return num_matches == 2

def classify_microcategory(rxn, supercateogry, subcategory):
    return classify_H(rxn), reaction_is_electron_transfer(rxn)
    
    return reaction_is_electron_transfer(rxn)
def find_reactant_charges(rxn):
    reactant_1_charge = find_mpculeid_charge(rxn.reactants[0])
    reactant_2_charge = find_mpculeid_charge(rxn.reactants[1])
    
    return reactant_1_charge, reactant_2_charge

def reaction_is_neutral(rxn):
    reactant_1_charge, reactant_2_charge = find_reactant_charges(rxn) 
       
    return reactant_1_charge == reactant_2_charge and reactant_1_charge == 0

def reaction_is_ion_molecule(rxn):
    reactant_1_charge, reactant_2_charge = find_reactant_charges(rxn)
    charges_are_different = reactant_1_charge != reactant_2_charge
    one_charge_zero = reactant_1_charge == 0 or reactant_2_charge == 0
    
    return charges_are_different and one_charge_zero
    
def classify_subcategory(rxn, supercategory):
    """
    For reactions with two reactants, determines if the reactions are
    oppositely charged. If they are, returns the ion-ion classification.

    Parameters
    ----------
    rxn : HiPRGen rxn object
        the reaction being tested

    Returns
    -------
    str
        "ion-ion" if reaction is ion-ion, None otherwise

    """
    
    if reaction_is_neutral(rxn):
        subcategory = "neutral"
    
    elif reaction_is_ion_molecule(rxn):
        subcategory = "ion-molecule"
        
    else:
        subcategory = "ion-ion"
        
    if supercategory != "bimolecular":
        return subcategory
    
    if classify_microcategory(rxn, subcategory):
        return classify_microcategory(rxn, subcategory)
    
    return supercategory
    
    # # if reaction_is_electron_transfer(rxn):
    # return classify_microcategory(rxn, "ion-ion")
    #     return "electron transfer"
    
    # return "ion-ion"
    # reactant_1 = rxn.reactants[0]
    # reactant_2 = rxn.reactants[1]
    
    # reactant_1_charge = find_mpculeid_charge(reactant_1)
    # reactant_2_charge = find_mpculeid_charge(reactant_2)

    # if reactant_1_charge == reactant_2_charge and reactant_1_charge == 0: #already filtered same charge ion reactions, this just removed
    #                                                                       #reactions where the charges are both 0
    #     return None

    # charges_are_opposites = -reactant_1_charge == reactant_2_charge
    
    # if charges_are_opposites:
        
    #     if len(rxn.products) == 2:
            
    #     return "ion-ion"
    
    # return None

def classify_ion_molecule(rxn):
    """
    Classifies a reaction as "ion-molecule" if, quite simply, one reactant
    is charged and the other isn't.

    Parameters
    ----------
    rxn : HiPRGen rxn object
        reaction we're classifying

    Returns
    -------
    str
        "ion-molecule" if reaction is ion-molecule, None otherwise

    """
    
    reactant_1 = rxn.reactants[0]
    reactant_2 = rxn.reactants[1]
    
    reactant_1_charge = find_mpculeid_charge(reactant_1)
    reactant_2_charge = find_mpculeid_charge(reactant_2)

    if reactant_1_charge != reactant_2_charge:
        return "ion-molecule"
   
    return None

def classify_neutral_radical(rxn):
    """
    Classifies a reaction as "neutral radical" if both reactants have charges
    of zero and one of the reactants has a spin of 2.

    Parameters
    ----------
    rxn : HiPRGen rxn object
        reaction we're classifying

    Returns
    -------
    str
        "neutral radical" if the reaction is neutral radical, None otherwise.

    """
    if not reaction_is_neutral(rxn):
        print(rxn)
        
        return None
    
    'fired'
    for reactant_mpculeid in rxn.reactants:
        
        spin = find_mpculeid_spin(reactant_mpculeid)
        
        if spin == 2:
            
            # reactant_1_charge = find_mpculeid_charge(rxn.reactants[0])
            # reactant_2_charge = find_mpculeid_charge(rxn.reactants[1])
            
            # reactant_charges_are_zero = \
            #     reactant_1_charge == reactant_2_charge and reactant_1_charge == 0
            
            # if reactant_charges_are_zero:
                
            return "neutral_radical"
                
    return "neutral"

def handle_combination_reactions(rxn):
    """
    Helper function for handling the case when a reaction has two reactants
    combining to form one product.

    Parameters
    ----------
    rxn : HiPRGen rxn object
        reaction we're classifying

    Returns
    -------
    string
        the classification of this particular type of combination reaction

    """
        
    possible_classifications = [
        classify_ion_ion(rxn),
        classify_ion_molecule(rxn),
        classify_neutral_radical(rxn)
     ]
    
    for classification in possible_classifications:
        if classification:
            return classification + "_combination"
        
    return "neutral_recombination"

def handle_fragmentation_reactions(rxn):
    """
    Helper function for classifying fragmentation reactions, where one reactant
    breaks apart into two products.

    Parameters
    ----------
    rxn : HiPRGen rxn object
        reaction we're classifying

    Returns
    -------
    string
        the classification of this particular type of fragmentation reaction

    """
    
    reverse_reaction_dict = \
        {"reactants":rxn.products, "products":rxn.reactants}
        
    reverse_reaction = HiPRGen_Reaction(reverse_reaction_dict)
    
    possible_classifications = [
        classify_ion_ion(reverse_reaction),
        classify_ion_molecule(reverse_reaction),
        classify_neutral_radical(reverse_reaction)
      ]
    
    for classification in possible_classifications:
        if classification:
            return classification + "_fragmentation"
    
    return "misc_fragmentation"

def handle_isomerization_reactions(rxn):
    """
    Helper function dealing with unimolecular isomerization reactions,
    classifying them based on whether they occur for a neutral species, an ion,
    or a radical

    Parameters
    ----------
    rxn : HiPRGen rxn object
        reaction we're classifying

    Returns
    -------
    string
        the classification of this particular type of isomerization reaction

    """
    reactant = rxn.reactants[0]
    product = rxn.products[0]
    
    reactant_charge = find_mpculeid_charge(reactant)
    product_charge = find_mpculeid_charge(product)
    reactant_spin = find_mpculeid_spin(reactant)
    product_spin = find_mpculeid_spin(product)
    
    if reactant_charge != product_charge:
        
        return "charge isomerization"
    
    if reactant_spin != 1 or product_spin != 1:
        
        return "radical isomerization"
    
    return "neutral isomerization"   
    
def handle_bimolecular_biproduct_reactions(rxn):
    """
    Helper function that classifies reactions with two reactants and two
    products

    Parameters
    ----------
    rxn : HiPRGen rxn object
        reaction we're classifying

    Returns
    -------
    string
        the classification of this reaction
    """
    
    possible_classifications = [
        classify_H(rxn),
        classify_ion_ion(rxn),
        classify_ion_molecule(rxn),
        classify_neutral_radical(rxn)
    ]

    for classification in possible_classifications:
        if classification:
            return classification

    return "misc_chemical"
    
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
    
    if len(rxn.reactants) == 2 and len(rxn.products) == 1:
        
        return handle_combination_reactions(rxn)
    
    if len(rxn.reactants) == 1 and len(rxn.products) == 2:
        
        return handle_fragmentation_reactions(rxn)
    
    if len(rxn.reactants) == 1 and len(rxn.products) == 1:
        
        return handle_isomerization_reactions(rxn)
          
    return handle_bimolecular_biproduct_reactions(rxn)