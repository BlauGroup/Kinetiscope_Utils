# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 14:24:02 2024

@author: jacob
"""
from reaction_classification_utilities import find_mpculeid_charge
from kinetiscope_reaction_writing_utilities import (
write_kinetiscope_name,
build_rxn_object,
build_reaction_parameter_dict)

def set_ionization_type_and_rate_constant_key(H_reaction):
    if "recombination" in H_reaction.tag:
        
        ionization_type = "recombination"
        rate_constant_key = "all"
        
    else:
        
        ionization_type = "attachment"
        
        if reaction_is_positive(H_reaction):
            
            rate_constant_key = "positive"
            
        else:
            
            rate_constant_key = "neutral"
    
    return ionization_type, rate_constant_key

def build_attachment_or_recombination(H_reaction, mpculeid_name_dict, rate_constant_dict):
    ionization_type, rate_constant_key = set_ionization_type_and_rate_constant_key(H_reaction) 
    product_list = H_reaction.classification_list[:-1]
    reactant_list = ["TE"]   
    reaction = build_ionization_reaction(H_reaction, mpculeid_name_dict, rate_constant_key, rate_constant_dict, ionization_type, 2, reactant_list, product_list)
    return [reaction]

def reaction_is_positive(H_reaction):
    reactant_mpculeid = H_reaction.reactants[0]
    reactant_charge = find_mpculeid_charge(reactant_mpculeid)
    reactant_is_positive = reactant_charge > 0
    return  reactant_is_positive  

def build_ionization_reaction(H_reaction, reaction_writing_data, rate_constant_key, ionization_type, reaction_order, reactants_to_add, products_to_add):
    """
    Ionization reactions can have reactants or products to add to the 
    Kinetiscope name depending on the reaction type. This function generates
    a name for this ionization reaction, finds the appropriate rate constant
    for the reaction, and finally builds and returns a Kinetiscope_Reaction
    object associated with this ionization reaction.

    Parameters
    ----------
    HiPRGen_reaction : HiPRGen reaction
        the reaction we're building kinetiscope reaction objects from
    mpculeid_dict : dict
        a dict of the form: {mpculeid:chemical name}
    rate_constant_dict : dict
        a nested dict created using the create_rate_constant_dict function from
        kinetiscope_reaction_writing_utilities
    rate_constant_key : string
        the key we'll look up in rate_constant_dict to find the associated
        rate constant
    ionization_type : string
        string related to finding where the rate_constant_key is nested in the
        rate_constant_dict
    reaction_order : int
        order of the chemical reaction
    reactants_to_add : list
        reactants to add to the kinetiscope name of the reaction
    products_to_add : list
        products to add to the kinetiscope name of the reaction

    Returns
    -------
    Kinetiscope_Reaction object
        the newly constructed Kinetiscope_Reaction object associated with this 
        reaction

    """
    
    products_to_add.append(ionization_type)
    kinetiscope_name = write_kinetiscope_name(H_reaction.name, reaction_writing_data, reactants_to_add, products_to_add)
    rate_constant =  reaction_writing_data.rate_constant_dict["ionization"][ionization_type].get(rate_constant_key, None)
    
    return build_rxn_object(H_reaction, kinetiscope_name, rate_constant, reaction_order, products_to_add)

def build_electron_ionization_reactions(H_reaction, reaction_writing_data):
    """
    Builds three electron ionization reactions for each HiPRGen positive
    ionization. 

    Parameters
    ----------
    HiPRGen_reaction : HiPRGen reaction
        the reaction we're building kinetiscope reaction objects from
    mpculeid_dict : dict
        a dict of the form: {mpculeid:chemical name}
    rate_constant_dict : dict
        a nested dict created using the create_rate_constant_dict function from
        kinetiscope_reaction_writing_utilities

    Returns
    -------
    electron_ionization_reaction_list : list
        list of Kinetiscope_Reaction objects representing electron ionization
        reactions for this species

    """
    
    electron_ionization_reaction_list = []
    electron_collision_dict = {
        "eV_80": ["eV_55", "LEE"],
        "eV_55": ["eV_30", "LEE"],
        "eV_30": ["2 LEE"]
    } 
    reaction_order = 2
    
    for reactant_electron, product_electron_list in electron_collision_dict.items():
        
        marker_species_without_tag = H_reaction.classification_list[:-1]
        full_product_list = marker_species_without_tag + product_electron_list
        
        electron_ionization_reaction = (
            build_ionization_reaction(H_reaction, reaction_writing_data, reactant_electron, "electron_ionization", reaction_order, [reactant_electron], full_product_list)
        )
        
        electron_ionization_reaction_list.append(electron_ionization_reaction)
        
    return electron_ionization_reaction_list


    
def build_absorption_reaction(H_reaction, reaction_writing_data):
    """
    Absorption of photons at the kinetic energies we're working with lead to
    the emission of photoelectrons with kinetic energies ~80 eV. This writes
    reactions modelling that process in a way that is readable by kinetiscope.

    Parameters
    ----------
    HiPRGen_reaction : HiPRGen reaction
        the reaction we're building kinetiscope reaction objects from
    mpculeid_dict : dict
        a dict of the form: {mpculeid:chemical name}
    rate_constant_dict : dict
        a nested dict created using the create_rate_constant_dict function from
        kinetiscope_reaction_writing_utilities

    Returns
    -------
    absorption_reaction : Kinetiscope_Reaction object
        a reaction representing absorption of a photon leading to the emission
        of a photoelectron.

    """
    def add_photoelectron_as_product(H_reaction):
        classifications_without_tag =  H_reaction.classification_list[:-1]
        kinetiscope_reaction_products = classifications_without_tag
        photoelectron = "eV_80"
        kinetiscope_reaction_products.append(photoelectron)
        return kinetiscope_reaction_products
    
    kinetiscope_reaction_products = add_photoelectron_as_product(H_reaction)
    absorber = H_reaction.reactants[0]
    reaction_order = 1
    
    reaction_dict = (
        build_reaction_parameter_dict(absorber, "absorption", reaction_order,products_to_add=kinetiscope_reaction_products)
    )
    
    absorption_reaction = (
        build_ionization_reaction(H_reaction, reaction_writing_data, reaction_dict)
    )
    
    return absorption_reaction
     
def build_positive_ionization_reactions(H_reaction, reaction_writing_data):
    """
    Each positive ionization reaction from HiPRGen is associated with several
    Kinetiscope reactions. Currently, if a species is a starting species in our
    simulations, it is considered an absorber and has an associated absorption
    reaction in addition to electron ionization reactions. Otherwise, if a
    species is not an absorber, it has only elecrron ionization reactions
    created for it.

    Parameters
    ----------
    HiPRGen_reaction : HiPRGen reaction
        the reaction we're building kinetiscope reaction objects from
    mpculeid_dict : dict
        a dict of the form: {mpculeid:chemical name}
    rate_constant_dict : dict
        a nested dict created using the create_rate_constant_dict function from
        kinetiscope_reaction_writing_utilities

    Returns
    -------
    list
        list of positive ionization reactions associated with H_reaction

    """
    def add_absorption_reaction(H_reaction, reaction_writing_data, positive_ionization_reactions):
        
        absorption_reaction = (
            build_absorption_reaction(H_reaction, reaction_writing_data)
        )
        
        positive_ionization_reactions.append(absorption_reaction)
        
        return positive_ionization_reactions
    
    def potential_absorber_absorbs(potential_absorber, reaction_writing_data):
        rate_constant_dict = reaction_writing_data.rate_constant_dict
        absorbers = rate_constant_dict["ionization"]["absorption"].keys()
        return potential_absorber in absorbers
        
    positive_ionization_reactions = []
    potential_absorber = H_reaction.reactants[0]
    
    reactant_is_absorber = (
        potential_absorber_absorbs(potential_absorber, reaction_writing_data)
    )
    
    if reactant_is_absorber:
        
        positive_ionization_reactions = (
            add_absorption_reaction(H_reaction, reaction_writing_data, positive_ionization_reactions)
        )
        
    electron_ionization_reactions = (
        build_electron_ionization_reactions(H_reaction, reaction_writing_data)
    )
    
    return positive_ionization_reactions + electron_ionization_reactions

def handle_ionization_reactions(HiPRGen_reaction, reaction_writing_data):
    """
    We handle naming positive ionization reactions differently from attachment
    and recombination reactions. This function just calls the appropriate method
    to build the Kinetiscope reactions related to this particular ionization
    reaction.

    Parameters
    ----------
    HiPRGen_reaction : HiPRGen reaction
        the reaction we're building kinetiscope reaction objects from
    reaction_writing_data: ReactionDataStorage object
        a class that stores data related to chemical reactions necessary to
        write their kinetiscope names, defined in 
        kinetiscope_reaction_writing_utilities

    Returns
    -------
    list
        a  list of Kinetiscope_Reaction objects associated with the 
        HiPRGen_reaction we called this function on.

    """
    
    if HiPRGen_reaction.tag == "positive_ionization":
        
        return build_positive_ionization_reactions(HiPRGen_reaction, reaction_writing_data)
    
    return build_attachment_or_recombination(HiPRGen_reaction, reaction_writing_data)