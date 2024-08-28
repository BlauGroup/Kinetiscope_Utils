# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 12:29:36 2024

@author: jacob
"""
from kinetiscope_reaction_writing_utilities import (
    determine_ordinal_number_order,
    replace_tag_with_shorthand,
    build_rxn_object,
    write_kinetiscope_name)

def write_name_with_excited_reactant(kinetiscope_name, reactant):
    """
    Modify the given kinetiscope name by appending an asterisk to the specified
    reactant.
    
    Parameters
    ----------
    kinetiscope_name : str
        The name of the kinetiscope reaction, whose form is defined in
        Rxn_classes
    reactant : str
        The reactant name to be marked with an asterisk.
    
    Returns
    -------
    str
        The modified kinetiscope name with the specified reactant marked by an
        asterisk.
    """
    
    reactant = reactant.strip()
    name_list = kinetiscope_name.split()
    new_name_list = []

    for name in name_list:
        
        if name == reactant:
            
            new_name_list.append(name + "*")
            
        else:
            
            new_name_list.append(name)
            
    return " ".join(new_name_list)

def add_reaction_with_excited_reactant(name_without_excitation, HiPRGen_rxn, reaction_list, reactant, reaction_writing_data):
    """
    Add a reaction with an excited reactant to the list of reactions.
    
    Parameters
    ----------
    name_without_excitation : str
        The name of the reaction without the excited reactant.
    HiPRGen_rxn : HiPRGen reaction object
        The reaction object from which to build the Kinetiscope reaction.
    reaction_list : list
        The list of reactions to which the new reaction will be appended.
    reactant : str
        The reactant to be excited in the reaction name.
    reaction_writing_data : ReactionDataStorage obj
        Contains information related to the reaction, including rate constants
        and marker species.
    
    Returns
    -------
    list
        The updated list of reactions, including the newly added reaction with
        the excited reactant.
    """
    
    ordinal_number_order = determine_ordinal_number_order(HiPRGen_rxn)
    order = 2 if ordinal_number_order == "2nd_order" else 1 
    
    rate_constant = (
        reaction_writing_data.rate_constant_dict["chemical"].get(ordinal_number_order, None)
        )
    
    # Generate the name for the reaction with the excited reactant
    
    name_with_excited_reactant = (
        write_name_with_excited_reactant(name_without_excitation, reactant)
    )
    
    # Replace tags with shorthand for marker species
    
    marker_species = (
        replace_tag_with_shorthand(HiPRGen_rxn.classification_list, reaction_writing_data.marker_species_dict)
    )
    
    # Build the reaction object using the generated name and rate constant
    
    excited_reaction = (
        build_rxn_object(HiPRGen_rxn, name_with_excited_reactant, rate_constant, order, marker_species)
        )
    
    reaction_list.append(excited_reaction)
    
    return reaction_list
    
def build_dexcitation_reaction(HiPRGen_rxn, dexcitation_reaction_name, reaction_list):
    """
    Builds a deexcitation reaction and appends it to the provided reaction 
    list.
    
    This function creates a deexcitation reaction using the specified name and 
    standard parameters. The newly created reaction is then appended to the 
    given reaction list, which is returned with the updated contents.
    
    Parameters
    ----------
    HiPRGen_rxn : HiPRGen reaction object
        The reaction object from which to build the Kinetiscope reaction.
    dexcitation_reaction_name : str
        The name to be assigned to the deexcitation reaction.
    reaction_list : list
        The list of reactions to which the new deexcitation reaction will be 
        appended.
    
    Returns
    -------
    list
        The updated reaction list with the appended deexcitation reaction.
    """
    
    order = 1
    rate_constant = 2.0E+06
    marker_species = ["dexcitation"]
    
    dexcitation_reaction = (
        build_rxn_object(HiPRGen_rxn, dexcitation_reaction_name, rate_constant, order, marker_species)
    )
    
    reaction_list.append(dexcitation_reaction)
    
    return reaction_list

def write_dexcitation_reaction_name(excitation_reaction_name, reactant_name):
    """
    Generates the name for a deexcitation reaction based on the provided e
    xcitation reaction name and reactant name.
    
    This function constructs the name for a deexcitation reaction by processing
    the given excitation reaction name. It identifies the relevant components
    related to the specified reactant and arranges them in reverse order,
    forming the deexcitation reaction name.
    
    Parameters
    ----------
    excitation_reaction_name : str
        The name of the excitation reaction from which the deexcitation 
        reaction name will be derived.
    reactant_name : str
        The name of the reactant involved in the deexcitation reaction.
    
    Returns
    -------
    str
        The generated name for the deexcitation reaction.
    """
    
    reactant_name = reactant_name.strip(" ") #removes whitespace from string
    excitation_species_list = excitation_reaction_name.split()
    
    dexcitation_list = (
        [name for name in excitation_species_list if reactant_name in name])
    
    dexcitation_list.reverse()
    
    return " => ".join(dexcitation_list)

def add_dexcitation_reaction(excitation_name, reactant, HiPRGen_rxn, reaction_writing_data, reaction_list):
    """
    Adds a deexcitation reaction to the reaction list based on the given excitation reaction name and reactant.
    
    This function generates the name for the deexcitation reaction using the provided excitation reaction name and reactant. It then builds the deexcitation reaction object and appends it to the reaction list.
    
    Parameters
    ----------
    excitation_name : str
        The name of the excitation reaction for which the deexcitation reaction is to be created.
    reactant : str
        The name of the reactant involved in the deexcitation reaction.
    HiPRGen_rxn : HiPRGen reaction object
        The reaction object used to build the deexcitation reaction, defined in Rxn_classes.
    reaction_writing_data : ReactionDataStorage object
        Contains information related to the reaction being processed, defined in kinetiscope_reaction_writing_utilities.
    reaction_list : list
        The list to which the constructed deexcitation reaction will be appended.
    
    Returns
    -------
    list
        The updated list of reactions including the newly added deexcitation reaction.
    """
    
    dexcitation_name = write_dexcitation_reaction_name(excitation_name, reactant)
    
    return build_dexcitation_reaction(HiPRGen_rxn, dexcitation_name, reaction_list)

def update_reaction_writing_data(reaction_writing_data, excitation_name):
    """
    Updates the reaction writing data by adding the specified excitation 
    reaction name to the excitation set.
    
    This function adds a new excitation reaction name to the excitation set 
    within the reaction writing data. It then returns the updated reaction 
    writing data.
    
    Parameters
    ----------
    reaction_writing_data : ReactionDataStorage object
        Contains information related to the reaction being processed,
        including the set of excitation reactions.
    excitation_name : str
        The name of the excitation reaction to be added to the excitation set.
    
    Returns
    -------
    ReactionDataStorage
        The updated reaction writing data object with the new excitation 
        reaction name added to the excitation set.
    """
    
    excitation_set = reaction_writing_data.excitation_set
    excitation_set.add(excitation_name)
    reaction_writing_data.excitation_set = excitation_set
    
    return reaction_writing_data

def add_excitation_reaction(HiPRGen_rxn, reaction_writing_data, excitation_name, reaction_list):
    """
    Adds an excitation reaction to the reaction list and updates the reaction 
    writing data.
    
    This function updates the reaction writing data by adding the specified 
    excitation reaction name to the excitation set, constructs an excitation 
    reaction object, and appends it to the provided reaction list.
    
    Parameters
    ----------
    HiPRGen_rxn : HiPRGen reaction object
        The reaction object used to build the Kinetiscope reaction.
    reaction_writing_data : ReactionDataStorage object
        Contains information related to the reaction being processed, 
        including the rate constant dictionary and excitation set.
    excitation_name : str
        The name of the excitation reaction to be added to the reaction list
        and the excitation set.
    reaction_list : list
        A list of reaction objects to which the new excitation reaction will 
        be appended.
    
    Returns
    -------
    tuple
        A tuple containing:
        - Updated list with the new excitation reaction appended.
        - Updated ReactionDataStorage object with the new excitation reaction 
        name added to the excitation set.
    """

    reaction_writing_data = (
        update_reaction_writing_data(reaction_writing_data, excitation_name)
    )
    
    order = 2
    
    rate_constant = (
        reaction_writing_data.rate_constant_dict["chemical"].get("2nd_order", None)
    )
    
    marker_species = ["excitation"]
    
    excitation_reaction = (
        build_rxn_object(HiPRGen_rxn, excitation_name, rate_constant, order,  marker_species)
    )
    
    reaction_list.append(excitation_reaction)
    
    return reaction_list, reaction_writing_data

def add_excitation_and_dexcitation(reactant, HiPRGen_rxn, reaction_writing_data, reaction_list, excitation_name):
    """
    Adds both an excitation and a deexcitation reaction to the reaction list 
    and updates the reaction writing data.
    
    This function first adds an excitation reaction to the reaction list and 
    updates the reaction writing data. Then, it adds a corresponding 
    deexcitation reaction to the list.
    
    Parameters
    ----------
    reactant : str
        The name of the reactant involved in the reactions.
    HiPRGen_rxn : HiPRGen reaction object
        The reaction object used to build the Kinetiscope reactions.
    reaction_writing_data : ReactionDataStorage object
        Contains information related to the reaction being processed, 
        including the rate constant dictionary and excitation set.
    reaction_list : list
        A list of reaction objects to which the new excitation and deexcitation
        reactions will be appended.
    excitation_name : str
        The name of the excitation reaction to be added to the reaction list 
        and used to generate the deexcitation reaction name.
    
    Returns
    -------
    tuple
        A tuple containing:
        - Updated list with the new excitation and deexcitation reactions 
        appended.
        - Updated ReactionDataStorage object with the new excitation reaction
        name added to the excitation set.
    """

    reaction_list, reaction_writing_data = (
        add_excitation_reaction(HiPRGen_rxn, reaction_writing_data, excitation_name, reaction_list)
    )
    
    reaction_list = (
        add_dexcitation_reaction(excitation_name, reactant, HiPRGen_rxn, reaction_writing_data, reaction_list)
    )
    
    return reaction_list, reaction_writing_data
 
def write_excitation_reaction_name(reactant):
    """
    Generates the name of an excitation reaction based on the reactant.
    
    This function formats a string to represent the excitation reaction where 
    the reactant is excited. The format follows the pattern: 
    `reactant + LEE => excited_reactant + TE + excitation`.
    
    Parameters
    ----------
    reactant : str
        The name of the reactant that will be excited in the reaction.
    
    Returns
    -------
    str
        The formatted name of the excitation reaction.
    """

    excited_reactant = f"{reactant}* "
    
    formated_excitation_name = (
        f"{reactant} + LEE => {excited_reactant} + TE + excitation"
    )
    
    return  formated_excitation_name
   
def add_excitation_dexcitation_if_new(reactant, HiPRGen_rxn, reaction_writing_data, reaction_list):
    """
    Adds an excitation and deexcitation reaction to the reaction list if the 
    excitation reaction is new.
    
    This function first generates the name for the excitation reaction based on
    the provided reactant. It then checks if this excitation reaction is 
    already in the set of recorded excitation reactions. If it is new, the 
    function adds both the excitation and deexcitation reactions to the 
    reaction list and updates the reaction writing data.
    
    Parameters
    ----------
    reactant : str
        The name of the reactant for which the excitation and deexcitation 
        reactions are being processed.
    HiPRGen_rxn : HiPRGen reaction object
        The HiPRGen reaction object used to build Kinetiscope reactions.
    reaction_writing_data : ReactionDataStorage obj
        Contains information related to Kinetiscope reactions, including the
        set of recorded excitation reactions.
    reaction_list : list
        The list of reactions to which new reactions will be added.
    
    Returns
    -------
    tuple
        A tuple containing the updated reaction list and the updated reaction
        writing data.
    """
    
    excitation_name = write_excitation_reaction_name(reactant)
    
    excitation_is_new = (
        excitation_name not in reaction_writing_data.excitation_set
    )
    
    if excitation_is_new:
        
       reaction_list, reaction_writing_data = (
           add_excitation_and_dexcitation(reactant, HiPRGen_rxn, reaction_writing_data, reaction_list, excitation_name)
       )
       
    return reaction_list, reaction_writing_data
 
def find_reactants_to_excite(name_without_excitation):
    """
    Extracts the list of reactants from a reaction name without excitation.
    
    This function takes a reaction name string (without excitation) and 
    extracts the list of reactants from it. The reactants are separated by 
    spaces and any '+' symbols are removed from the list.
    
    Parameters
    ----------
    name_without_excitation : str
        The name of the reaction without excitation, where reactants are listed
        before the '=>' symbol.
    
    Returns
    -------
    list
        A list of reactant names extracted from the provided reaction name 
        string.
    """

    reactants_str = name_without_excitation.split(" => ")[0]
    reactant_list = reactants_str.split()
    
    if "+" in reactant_list:
        reactant_list.remove("+")
        
    return reactant_list

def write_name_without_excitation(HiPRGen_rxn, reaction_writing_data):
    """
    Generates the Kinetiscope reaction name from the HiPRGen reaction object.
    
    This function creates a Kinetiscope-compatible reaction name by replacing 
    tags with shorthand notations from the HiPRGen reaction object. 
    It constructs the name based on the classification list and uses the
    provided shorthand dictionary to format the reaction name.
    
    Parameters
    ----------
    HiPRGen_rxn : HiPRGen reaction object
        The reaction object containing classification information and the 
        original name.
    reaction_writing_data : ReactionDataStorage object
        Contains the shorthand dictionary for marker species and other 
        reaction-related data.
    
    Returns
    -------
    str
        The formatted Kinetiscope reaction name with shorthand notation for the
        species.
    """

    species_list = HiPRGen_rxn.classification_list[:]
    shorthand_dict = reaction_writing_data.marker_species_dict
    
    shorthand_species = (
        replace_tag_with_shorthand(species_list, shorthand_dict)
    )
    
    added_reactants = None  #currently no marker species added as reactants
    HiPRGen_name = HiPRGen_rxn.name
    
    kinetiscope_name = write_kinetiscope_name(
        HiPRGen_name, shorthand_dict, added_reactants, shorthand_species
    )
        
    return kinetiscope_name

def write_phase_1_chemical_reactions(HiPRGen_rxn, reaction_writing_data):
    """
    Processes Phase 1 chemical reactions by generating excitation and 
    deexcitation reactions.
    
    This function handles the construction of excitation and deexcitation 
    reactions for a given HiPRGen reaction. It creates a list of reactions 
    based on the reactants that need excitation, adds new reactions if they 
    are not already present, and updates the reaction writing data.
    
    Parameters
    ----------
    HiPRGen_rxn : HiPRGen reaction object
        The reaction object from which to derive the Kinetiscope reactions.
    reaction_writing_data : ReactionDataStorage object
        Contains information about existing reactions and rate constants.
    
    Returns
    -------
    tuple
        A tuple containing:
        - list : The updated list of reactions with new excitation and 
        deexcitation reactions added.
        - ReactionDataStorage : The updated reaction writing data with new 
        excitation information.
    """
    
    reaction_list = []
    
    name_without_excitation = (
        write_name_without_excitation(HiPRGen_rxn, reaction_writing_data)
    )
    
    reactants_to_excite = find_reactants_to_excite(name_without_excitation)
    
    for reactant in reactants_to_excite:
        
        reaction_list, reaction_writing_data = (
            add_excitation_dexcitation_if_new(reactant, HiPRGen_rxn, reaction_writing_data, reaction_list)
        )
        
        reaction_list = (
            add_reaction_with_excited_reactant(name_without_excitation, HiPRGen_rxn, reaction_list, reactant, reaction_writing_data)
        )
        
    return reaction_list, reaction_writing_data