# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 11:15:59 2024

@author: jacob
"""

from kinetiscope_reaction_writing_utilities import build_ionization_reaction

def add_products_to_list(HiPRGen_reaction, product_list):
    """
    Helper function for adding products to the list we'll put into the names
    of reactions.

    Parameters
    ----------
    HiPRGen_reaction : HiPRGen reaction object
        the reaction we're building kinetiscope reaction objects from, defined 
        in Rxn_classes
    product_list : list
        list of products we're going to add

    Returns
    -------
    kinetiscope_reaction_products : list
        classifications_without_tag with the new products added

    """
    classifications_without_tag = HiPRGen_reaction.classification_list[:-1]
    kinetiscope_reaction_products = classifications_without_tag
    
    for product in product_list:
        kinetiscope_reaction_products.append(product)
    
    return kinetiscope_reaction_products

def build_EI_reaction_dict(HiPRGen_reaction, reactant_electron, product_electron_list):
    """
    As for building absorption reactions, we have a lot of data to pass to the
    construction of these electron ionization reactions, so this just puts
    those parameters in a dictionary for readability.

    Parameters
    ----------
    HiPRGen_reaction : HiPRGen reaction object
        the reaction we're building kinetiscope reaction objects from, defined 
        in Rxn_classes
    reactant_electron : str
        the electron (or photoelectron) undergoing collision
    product_electron_list : list
        the electrons that are products of this particular collision--one is
        the reactant_electron with less kinetic energy, while the other is
        a low-energy electron released from the target following collision

    Returns
    -------
    reaction_dict : dict
        a dictionary with information used to construct the Kinetiscope_reaction
        obj

    """
    superclass = "electron_ionization"
    
    full_product_list = (
        add_products_to_list(HiPRGen_reaction, product_electron_list)
    )
    
    reaction_dict = { 
        "rate_constant_key":reactant_electron,
        "superclass":superclass,
        "reaction_order":2,
        "reactants_to_add":[reactant_electron],
        "products_to_add":full_product_list,
    }
    
    return reaction_dict

def build_one_EI_reaction(HiPRGen_reaction, reactant_electron, product_electron_list, reaction_writing_data):
    """
    Creates a single Kinetiscope_Reaction object associated with a given
    electron ionization reaction.

    Parameters
    ----------
    HiPRGen_reaction : HiPRGen reaction object
        the reaction we're building kinetiscope reaction objects from, defined 
        in Rxn_classes
    reactant_electron : string
        string representing the electron whose collision with a particular
        species we're modelling in this reaction
    product_electron_list : list
        list of electrons (represented by string) produced by the collision
    reaction_writing_data: ReactionDataStorage object
        a class that stores data related to chemical reactions necessary to
        write their kinetiscope names, defined in 
        kinetiscope_reaction_writing_utilities

    Returns
    -------
    EI_reaction : Kinetiscope_Reaction object
        a reaction associated with this particular ionization reaction

    """
    
    reaction_dict = (
        build_EI_reaction_dict(HiPRGen_reaction, reactant_electron, product_electron_list)
    )

    return build_ionization_reaction(HiPRGen_reaction, reaction_writing_data, reaction_dict)

def build_all_EI_reactions(HiPRGen_reaction, reaction_writing_data):
    """
    Builds three electron ionization reactions for each HiPRGen positive
    ionization. 

    Parameters
    ----------
    HiPRGen_reaction : HiPRGen reaction object
        the reaction we're building kinetiscope reaction objects from, defined 
        in Rxn_classes
        a class that stores data related to chemical reactions necessary to
        write their kinetiscope names, defined in 
        kinetiscope_reaction_writing_utilities

    Returns
    -------
    EI_reaction_list : list
        list of Kinetiscope_Reaction objects representing electron ionization
        reactions for this species

    """
    
    all_EI_reactions = []
    
    electron_collision_dict = {
        "eV_80": ["eV_55", "LEE"],
        "eV_55": ["eV_30", "LEE"],
        "eV_30": ["2 LEE"]
    } 
    
    for reactant_electron, product_electron_list in electron_collision_dict.items():
        
        EI_reaction = (
            build_one_EI_reaction(HiPRGen_reaction, reactant_electron, product_electron_list, reaction_writing_data)
        )
        
        all_EI_reactions.append(EI_reaction)
        
    return all_EI_reactions

def build_absorption_reaction_dict(HiPRGen_reaction):
    """
    We have to pass a lot of info to build our Kinetiscope_reaction objects,
    so this just constructs it as a dict to make the function call more
    readable

    Parameters
    ----------
    HiPRGen_reaction : HiPRGen_reaction object
        class defined in Rxn_objects. Contains information related to
        a particular reaction found via HiPRGen

    Returns
    -------
    reaction_dict : dict
        a dictionary containing info we'll pass to construct the absorption
        reaction

    """
    
    superclass = "absorption"
    photoelectron = "eV_80" #photoelectron has 80 eV kinetic energy
    product_list = [photoelectron]
    products_to_add = add_products_to_list(HiPRGen_reaction, product_list)
    absorber = HiPRGen_reaction.reactants[0]
    
    reaction_dict = {
        "rate_constant_key":absorber,
        "superclass":superclass,
        "reaction_order":1,
        "reactants_to_add":None,
        "products_to_add":products_to_add
    }
    
    return reaction_dict

def build_absorption_reaction(HiPRGen_reaction, reaction_writing_data):
    """
    Absorption of photons at the kinetic energies we're working with leads to
    the emission of photoelectrons with kinetic energies ~80 eV. This writes
    reactions modelling that process in a way that is readable by kinetiscope.

    Parameters
    ----------
    HiPRGen_reaction : HiPRGen reaction object
        the reaction we're building kinetiscope reaction objects from, defined 
        in Rxn_classes
    reaction_writing_data: ReactionDataStorage object
        a class that stores data related to chemical reactions necessary to
        write their kinetiscope names, defined in 
        kinetiscope_reaction_writing_utilities

    Returns
    -------
    Kinetiscope_Reaction object
        a reaction representing absorption of a photon leading to the emission
        of a photoelectron.

    """
    
    reaction_dict = build_absorption_reaction_dict(HiPRGen_reaction)
    
    return build_ionization_reaction(HiPRGen_reaction, reaction_writing_data, reaction_dict)

def species_is_absorber(potential_absorber, reaction_writing_data):
    """
    A function that checks if a given species is able to absorb in our 
    simulations by checking if it is in the dict.keys() view object 
    representing all absorbers we're modelling.

    Parameters
    ----------
    potential_absorber : string
        an mpculeid of the form graph_hash-formula-charge-spin
    reaction_writing_data: ReactionDataStorage object
        a class that stores data related to chemical reactions necessary to
        write their kinetiscope names, defined in 
        kinetiscope_reaction_writing_utilities

    Returns
    -------
    Bool
        True if the potential_absorber is an absorber, False otherwise

    """
    
    rate_constant_dict = reaction_writing_data.rate_constant_dict
    absorbers = rate_constant_dict["ionization"]["absorption"].keys()
    
    return potential_absorber in absorbers
 
def add_absorption_if_absorber(HiPRGen_reaction, reaction_writing_data):
    """
    Adds an absorption reaction to the positive_ionization_reactions list
    only if the reactant of this HiPRGen reaction is in the list of absorbers.
    Otherwise, just returns an empty list.

    Parameters
    ----------
    HiPRGen_reaction : HiPRGen reaction object
        the reaction we're building kinetiscope reaction objects from, defined 
        in Rxn_classes
     reaction_writing_data: ReactionDataStorage object
         a class that stores data related to chemical reactions necessary to
         write their kinetiscope names, defined in 
         kinetiscope_reaction_writing_utilities

    Returns
    -------
    positive_ionization_reactions : list
        a list containing 1 absorption reaction if the reactant is an absorber,
        an empty list otherwise

    """
    
    positive_ionization_reactions = []
    potential_absorber = HiPRGen_reaction.reactants[0]
    
    
    if species_is_absorber(potential_absorber, reaction_writing_data):
        
        absorption_reaction = (
            build_absorption_reaction(HiPRGen_reaction, reaction_writing_data)
        )
        
        positive_ionization_reactions.append(absorption_reaction)
        
    return positive_ionization_reactions

def build_positive_ionization_reactions(HiPRGen_reaction, reaction_writing_data):
    """
    Each positive ionization reaction from HiPRGen is associated with several
    Kinetiscope reactions. Currently, if a species is a starting species in our
    simulations, it is considered an absorber and has an associated absorption
    reaction in addition to electron ionization (EI) reactions. Otherwise, if a
    species is not an absorber, it has only EI reactions created.

    Parameters
    ----------
    HiPRGen_reaction : HiPRGen reaction object
        the reaction we're building kinetiscope reaction objects from, defined 
        in Rxn_classes
    reaction_writing_data: ReactionDataStorage object
        a class that stores data related to chemical reactions necessary to
        write their kinetiscope names, defined in 
        kinetiscope_reaction_writing_utilities

    Returns
    -------
    list
        list of Kinetiscope positive ionization reactions associated 
        with this HiPRGen_reaction

    """
    
    positive_ionization_reactions = (
        add_absorption_if_absorber(HiPRGen_reaction, reaction_writing_data)
    )  
    
    EI_reactions =(
        build_all_EI_reactions(HiPRGen_reaction, reaction_writing_data)
    )
    
    return positive_ionization_reactions + EI_reactions