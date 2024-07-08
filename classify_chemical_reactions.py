# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:54:29 2024

@author: JRMilton
"""
import re
from classify_ionization_reactions import find_charge
from Rxn_classes import HiPRGen_Reaction

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
    
    reactant_charge = find_charge(reactant_gaining_H)
    product_charge = find_charge(product_with_H)
    
    delta_charge = product_charge - reactant_charge
    
    if delta_charge == 1:
        return "proton_transfer"
    
    elif delta_charge == -1:
        return "hydride_abstraction"
    
    is_radical_reaction = \
        reactant_gaining_H.split("-")[3] != product_with_H.split("-")[3]
    
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
        
        reactant_formula = reactant_mpculeid.split("-")[1]
        reactant_with_one_more_hydrogen = add_one_more_hydrogen(reactant_formula)
        
        for product_mpculeid in rxn.products:
            
            product_formula = product_mpculeid.split("-")[1]
            
            if product_formula == reactant_with_one_more_hydrogen:
                return narrow_H_rxn_type(reactant_mpculeid, product_mpculeid)
                
    return None
            
def classify_ion_ion(rxn):
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
    
    reactant_1 = rxn.reactants[0]
    reactant_2 = rxn.reactants[1]
    
    reactant_1_charge = find_charge(reactant_1)
    reactant_2_charge = find_charge(reactant_2)

    if reactant_1_charge == reactant_2_charge and reactant_1_charge == 0: #already filtered same charge ion reactions, this just removed
                                                                          #reactions where the charges are both 0
        return None

    charges_are_opposites = -reactant_1_charge == reactant_2_charge
    
    if charges_are_opposites:
        return "ion-ion"
    
    return None

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
    
    reactant_1_charge = find_charge(reactant_1)
    reactant_2_charge = find_charge(reactant_2)

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

    for reactant_mpculeid in rxn.reactants:
        spin = int(reactant_mpculeid.split('-')[3])
        if spin == 2:
            
            reactant_1_charge = find_charge(rxn.reactants[0])
            reactant_2_charge = find_charge(rxn.reactants[1])
            
            reactant_charges_are_zero = reactant_1_charge == reactant_2_charge and reactant_1_charge == 0
            
            if reactant_charges_are_zero:
                return "neutral_radical"
                
    return None

def tag_chemical_reaction(rxn):
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
    
    if len(rxn.reactants) == 2 and len(rxn.products) == 1: #combination reactions
        
        classifications = [
            classify_ion_ion(rxn),
            classify_ion_molecule(rxn),
            classify_neutral_radical(rxn)
         ]
        
        for classification in classifications:
            if classification:
                return classification + "_combination"
            
        return "misc_neutral_recombination"
    
    if len(rxn.reactants) == 1 and len(rxn.products) == 2: #fragmentation reactions
        
        reverse_reaction_dict = \
            {"reactants":rxn.products, "products":rxn.reactants}
            
        reverse_reaction = HiPRGen_Reaction(reverse_reaction_dict)
        
        classifications = [
            classify_ion_ion(reverse_reaction),
            classify_ion_molecule(reverse_reaction),
            classify_neutral_radical(reverse_reaction)
          ]
        
        for classification in classifications:
            if classification:
                return classification + "_fragmentation"
        
        return "misc_fragmentation"
    
    if len(rxn.reactants) == 2: #2 reactant, 2 product reactions
        
        classifications = [
            classify_H(rxn),
            classify_ion_ion(rxn),
            classify_ion_molecule(rxn),
            classify_neutral_radical(rxn)
        ]
    
        for classification in classifications:
            if classification:
                return classification

    return "misc_chemical"  # tag if no classification matched