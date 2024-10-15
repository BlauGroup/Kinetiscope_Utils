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
import copy
import sys


def generate_mpculeid_formula_charge_dict(rxn):
    """
    Generates a dictionary mapping molecular IDs (mpculeid) to their formulas
    and charges for both reactants and products.

    This function creates a dictionary where the keys are 'reactants' and
    'products', and each value is another dictionary that maps molecular
    IDs to their associated formula and charge.

    Parameters
    ----------
    rxn : object
        The reaction object containing reactants and products.

    Returns
    -------
    dict
        A dictionary containing the molecular IDs of reactants and products
        mapped to their formulas and charges.
    """

    def generate_formula_charge_dict(mpculeid):
        """
        Generates a dictionary containing the formula and charge of a molecule
        based on its molecular ID.

        Parameters
        ----------
        mpculeid : str
            The molecular ID of the molecule.

        Returns
        -------
        dict
            A dictionary containing the formula and charge of the molecule.
        """
        formula = find_mpculeid_formula(mpculeid)
        charge = find_mpculeid_charge(mpculeid)
        return {"formula": formula, "charge": charge}

    mpculeid_formula_charge_dict = {
        "reactants": {},
        "products": {}
    }

    for reactant_mpculeid in rxn.reactants:
        mpculeid_formula_charge_dict["reactants"][
            reactant_mpculeid
        ] = generate_formula_charge_dict(reactant_mpculeid)

    for product_mpculeid in rxn.products:
        mpculeid_formula_charge_dict["products"][
            product_mpculeid
        ] = generate_formula_charge_dict(product_mpculeid)

    return mpculeid_formula_charge_dict


def find_matching_product(reactant, products):
    """
    Finds a product with a matching formula and a charge differing by one.
    """

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

    for product_mpculeid, product_info in products.items():

        matching_formula = is_matching_formula(
            reactant["formula"], product_info["formula"]
        )

        charges_differ = charges_differ_by_one(
            reactant["charge"], product_info["charge"]
        )

        if matching_formula and charges_differ:

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

    def get_sorted_formulas(mpculeid_formula_charge_dict, category):
        """
        Extracts and sorts formulas from the specified category (reactants or
        products) of the mpculeid_formula_charge_dict.
        """

        formulas = []

        for entry in mpculeid_formula_charge_dict[category].values():

            formulas.append(entry["formula"])

        formulas.sort()

        return formulas

    mpculeid_formula_charge_dict = generate_mpculeid_formula_charge_dict(rxn)

    sorted_reactant_formulas = get_sorted_formulas(
        mpculeid_formula_charge_dict, "reactants"
    )

    sorted_product_formulas = get_sorted_formulas(
        mpculeid_formula_charge_dict, "products"
    )

    if sorted_reactant_formulas != sorted_product_formulas:
        return ""

    test_dict = copy.deepcopy(mpculeid_formula_charge_dict)
    num_matches = 0

    for reactant_mpculeid, reactant_info in test_dict["reactants"].items():

        matching_product_mpculeid = find_matching_product(
            reactant_info, test_dict["products"]
        )

        if matching_product_mpculeid:
            num_matches += 1
            test_dict["products"].pop(matching_product_mpculeid)

    if num_matches == 2:
        return "electron_transfer"

    return ""


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


def add_one_more_hydrogen(formula):
    """
    Increments the number of hydrogens in the given formula by 1.
    If there are no hydrogens in the formula, adds 'H1' in alphabetical
    order.

    Parameters
    ----------
    formula : str
        String of the form atom1#atom2# etc (e.g., 'C2H6O').

    Returns
    -------
    modified_formula : str
        The formula with the number of hydrogens increased by one.
    """
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
        Some formulas we need to modify don't contain any hydrogens. Thus, we
        need to add hydrogen to these formulas, but because pymatgen generates
        formulas in alphabetical order, we need to insert hydrogen where it
        belongs forcomparison. This function does that.

        Parameters
        ----------
        formula : str
            String of the form atom1#atom2# etc (e.g., 'C2H6O').

        Returns
        -------
        modified_formula : str
            formula with one hydrogen added

        """
        # returns  list of tuples (atom, #)
        elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula)

        elements.append(('H', '1'))
        elements.sort(key=lambda x: x[0])  # sorts list by first ele of tuples
        modified_formula = ''.join(f"{el}{count}" for el, count in elements)

        return modified_formula

    number_hydrogens = count_number_hydrogens(formula)

    if number_hydrogens == 0:

        modified_formula = add_hydrogen_to_formula(formula)

    else:

        modified_hydrogen_count = number_hydrogens + 1

        modified_formula = re.sub(
            r'H\d*', f'H{modified_hydrogen_count}', formula
        )

    return modified_formula


def classify_H(rxn):
    """
    Determines whether a given reaction involves the transfer of a hydrogen
    atom by comparing the reactants of the reaction to the products. If a
    product has one more hydrogen atom than a reactant, the reaction is an H
    atom transfer.Then, we narrow down the type via narrow_H_rxn_type and
    return the specific type.

    Parameters
    ----------
    rxn : HiPRGen rxn object
        The reaction we're testing.

    Returns
    -------
    str or None
        Narrowed H reaction type if the reaction is an H transfer,
        None otherwise.
    """
    def handle_multiple_classifications(possible_classifications):
        # tested and possible_classifications is always len(2)
        if possible_classifications[0] == possible_classifications[1]:
            print("same classification")
            return possible_classifications[0]

        classification_hierarchy = [
            "proton_transfer",
            "H_atom_abstraction",
            "hydride_abstraction",
            "proton_coupled_electron_transfer"]

        for classification in classification_hierarchy:
            if classification in possible_classifications:
                print("different classifications")
                return classification

    possible_classifications = []

    for reactant_mpculeid in rxn.reactants:

        reactant_formula = find_mpculeid_formula(reactant_mpculeid)

        reactant_with_one_more_hydrogen = add_one_more_hydrogen(
            reactant_formula
        )

        for product_mpculeid in rxn.products:

            product_formula = find_mpculeid_formula(product_mpculeid)

            if product_formula == reactant_with_one_more_hydrogen:

                possible_classification = narrow_H_rxn_type(
                    reactant_mpculeid, product_mpculeid
                )

                possible_classifications.append(possible_classification)

    if not possible_classifications:

        return ""

    elif len(possible_classifications) == 1:
        return possible_classifications[0]

    else:
        return handle_multiple_classifications(possible_classifications)


def handle_biproduct_reactions(rxn, bireactant_subclass):
    """
    Handles biproduct reactions, checking for hydrogen classification and
    electron transfer, and returns the appropriate subclass.
    """

    H_name = classify_H(rxn)

    if H_name:

        return bireactant_subclass + "_" + H_name

    electron_transfer_name = reaction_is_electron_transfer(rxn)

    if electron_transfer_name:

        return bireactant_subclass + "_" + electron_transfer_name

    return bireactant_subclass + "_reaction"


def determine_bireactant_subclass(rxn):
    """
    Determines the subclass for a bireactant reaction based on the species
    subclasses of the reactants.

    This function assigns a subclass to a bireactant reaction by first
    determining the species subclasses of both reactants and then generating
    an ordered name for the reaction based on the priority of the subclasses.

    Parameters
    ----------
    rxn : object
        The reaction object with two reactants.

    Returns
    -------
    str
        The ordered subclass name of the bireactant reaction.
    """

    def generate_ordered_name(reactant_1_subclass, reactant_2_subclass):
        """
        Generates an ordered name for the reaction based on the priority of the
        reactant subclasses.

        The function uses a priority dictionary to ensure that the higher
        priority subclass comes first in the name.

        Parameters
        ----------
        reactant_1_subclass : str
            The subclass of the first reactant.
        reactant_2_subclass : str
            The subclass of the second reactant.

        Returns
        -------
        str
            The ordered name of the reaction, with the higher priority subclass
            listed first.
        """

        priority_dict = {
            "radical_anion": 5,
            "radical_cation": 4,
            "anion": 3,
            "cation": 2,
            "neutral_radical": 1,
            "neutral": 0
        }

        reactant_1_priority = priority_dict.get(reactant_1_subclass)
        reactant_2_priority = priority_dict.get(reactant_2_subclass)

        if reactant_1_priority > reactant_2_priority:

            ordered_name = reactant_1_subclass + "_" + reactant_2_subclass

        else:

            ordered_name = reactant_2_subclass + "_" + reactant_1_subclass

        return ordered_name

    reactant_1 = rxn.reactants[0]
    reactant_2 = rxn.reactants[1]

    reactant_1_subclass = determine_species_subclass(reactant_1)
    reactant_2_subclass = determine_species_subclass(reactant_2)

    ordered_name = generate_ordered_name(
        reactant_1_subclass, reactant_2_subclass
    )

    return ordered_name


def handle_bireactant_reactions(rxn):
    """
    Handles classification of bimolecular reactions and returns the appropriate
    subclass.
    """
    def is_combination_reaction(rxn):
        """
        Check if the reaction is a combination reaction based on the number of
        products
        """
        return len(rxn.products) == 1

    bireactant_subclass = determine_bireactant_subclass(rxn)

    if is_combination_reaction(rxn):

        return bireactant_subclass + "_combination"

    return handle_biproduct_reactions(rxn, bireactant_subclass)


def determine_species_subclass(mpculeid):
    """
    Determines the species subclass based on the charge and spin of the
    molecule.

    This function first determines the charge of the molecule and assigns
    it a charge-related name ('cation', 'anion', or 'neutral'). It then
    checks the spin of the molecule to determine whether it is a radical,
    and appends the appropriate prefix.

    Parameters
    ----------
    mpculeid : str
        The identifier for the molecule.

    Returns
    -------
    str
        The subclass of the species, which could be a combination of its
        charge and whether it is a radical.
    """

    def determine_charge_name(mpculeid):
        """
        Determines the charge name of the molecule ('cation', 'anion', or
        'neutral') based on the charge.

        Parameters
        ----------
        mpculeid : str
            The identifier for the molecule.

        Returns
        -------
        str
            The charge name ('cation', 'anion', or 'neutral').
        """

        charge = find_mpculeid_charge(mpculeid)

        if charge > 0:
            return "cation"
        elif charge < 0:
            return "anion"
        else:
            return "neutral"

    charge_name = determine_charge_name(mpculeid)
    spin = find_mpculeid_spin(mpculeid)

    if spin == 2 and charge_name == "neutral":
        return charge_name + "_radical"

    elif spin == 2:
        return "radical_" + charge_name

    return charge_name


def handle_unireactant_reactions(rxn):
    """
    Handles classification of unireactant reactions and returns the appropriate
    subclass.

    This function takes a unireactant reaction and determines whether it is
    an isomerization or fragmentation based on the number of products.

    Parameters
    ----------
    rxn : object
        A reaction object with reactants and products.

    Returns
    -------
    str
        The classified tag of the unireactant reaction.
    """

    single_reactant_mpculeid = rxn.reactants[0]
    species_subclass = determine_species_subclass(single_reactant_mpculeid)

    if len(rxn.products) == 1:

        return species_subclass + "_isomerization"

    return species_subclass + "_fragmentation"


def determine_chemical_reaction_tag(rxn):
    """
    Tags chemical reactions by calling classification functions in a specific
    order.

    Parameters
    ----------
    rxn : HiPRGen rxn object
        Reaction we're classifying.

    Returns
    -------
    list
        The tag for the classified reaction.
    """

    def is_unireactant_reaction(rxn):
        """
        Check if the reaction is a unireactant reaction based on the number
        of reactants.
        """
        return len(rxn.reactants) == 1

    if is_unireactant_reaction(rxn):

        return handle_unireactant_reactions(rxn)

    return handle_bireactant_reactions(rxn)
