# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 14:24:02 2024

@author: jacob
"""
from reaction_classification_utilities import find_mpculeid_charge
from kinetiscope_reaction_writing_utilities import (
write_kinetiscope_name,
build_rxn_object
)

def set_ionization_type_and_rate_constant_key(HiPRGen_reaction):
    if "recombination" in HiPRGen_reaction.tag:
        
        ionization_type = "recombination"
        rate_constant_key = "all"
        
    else:
        
        ionization_type = "attachment"
        
        if reaction_is_positive(HiPRGen_reaction):
            
            rate_constant_key = "positive"
            
        else:
            
            rate_constant_key = "neutral"
    
    return ionization_type, rate_constant_key

def build_attachment_or_recombination(HiPRGen_reaction, mpculeid_name_dict, rate_constant_dict):
    ionization_type, rate_constant_key = set_ionization_type_and_rate_constant_key(HiPRGen_reaction) 
    product_list = HiPRGen_reaction.classification_list[:-1]
    reactant_list = ["TE"]   
    reaction = build_ionization_reaction(HiPRGen_reaction, mpculeid_name_dict, rate_constant_key, rate_constant_dict, ionization_type, 2, reactant_list, product_list)
    return [reaction]

def reaction_is_positive(HiPRGen_reaction):
    reactant_mpculeid = HiPRGen_reaction.reactants[0]
    reactant_charge = find_mpculeid_charge(reactant_mpculeid)
    reactant_is_positive = reactant_charge > 0
    return  reactant_is_positive  

def find_ionization_rate_constant(reaction_writing_data, reaction_dict, ionization_type):
    """
    Each Kineticope_Reaction object needs a rate constant assigned to it for
    our simulations. This function finds the rate constant associated with this
    particular ionization reaction and returns it for assignment.

    Parameters
    ----------
    reaction_writing_data: ReactionDataStorage object
        a class that stores data related to chemical reactions necessary to
        write their kinetiscope names, defined in 
        kinetiscope_reaction_writing_utilities
    reaction_dict : dict
        a dictionary, built using the build_reaction_parameter_dict function
        in kinetiscope_reaction_writing_utilities, which contains information
        associated with this particular kinetiscope reaction
    ionization_type : string
        A string, describing the ionization superclass this reaction belongs to

    Returns
    -------
    float
        a floating point, written in scientific notation, representing the rate
        constant associated with this particular reaction. Units for first-order
        rate constants are s^-1 while for second-order they are M^-1 s^-1.

    """
    
    rate_constant_key = reaction_dict["rate_constant_key"]
    all_rate_constants = reaction_writing_data.rate_constant_dict
    ionization_rate_constants = all_rate_constants["ionization"][ionization_type]
        
    return ionization_rate_constants.get(rate_constant_key, None)

def write_ionization_name(HiPRGen_reaction, reaction_writing_data, reaction_dict, ionization_type, products_to_add):
    """
    Generates a name for an ionization reaction to be assigned to a 
    Kinetiscope_Reaction object.

    Parameters
    ----------
    HiPRGen_reaction : HiPRGen reaction object
        the reaction we're building kinetiscope reaction objects from, defined 
        in Rxn_classes
    reaction_writing_data: ReactionDataStorage object
        a class that stores data related to chemical reactions necessary to
        write their kinetiscope names, defined in 
        kinetiscope_reaction_writing_utilities
    reaction_dict : dict
        a dictionary, built using the build_reaction_parameter_dict function
        in kinetiscope_reaction_writing_utilities, which contains information
        associated with this particular kinetiscope reaction
    ionization_type : string
        A string, describing the ionization superclass this reaction belongs to,
        to be added as a product of the ionization reaction
    products_to_add : list
        a list of products of this reaction to which ionization_type will be
        appended

    Returns
    -------
    string
        a name for this kinetiscope reaction, written using 
        write_kinetiscope_name defined in 
        kinetiscope_reaction_writing_utilities

    """
    
    mpculeid_dict = reaction_writing_data.mpculeid_dict
    reactants_to_add = reaction_dict["reactants_to_add"]    
    H_name = HiPRGen_reaction.name
    
    return write_kinetiscope_name(H_name, mpculeid_dict, reactants_to_add, products_to_add)
   
def build_ionization_reaction(HiPRGen_reaction, reaction_writing_data, reaction_dict):
    """
    Ionization reactions can have reactants or products to add to the 
    Kinetiscope name depending on the reaction type. This function generates
    a name for this ionization reaction, finds the appropriate rate constant
    for the reaction, and finally builds and returns a Kinetiscope_Reaction
    object associated with this ionization reaction.

    Parameters
    ----------
    HiPRGen_reaction : HiPRGen reaction object
        the reaction we're building kinetiscope reaction objects from, defined 
        in Rxn_classes
    reaction_writing_data: ReactionDataStorage object
        a class that stores data related to chemical reactions necessary to
        write their kinetiscope names, defined in 
        kinetiscope_reaction_writing_utilities
    reaction_dict : dict
        a dictionary, which contains information associated with this
        particular kinetiscope reaction. This is build before this function is
        called so that the parameters passed aren't unwieldy

    Returns
    -------
    ionization_reaction: Kinetiscope_Reaction object
        the newly constructed Kinetiscope_Reaction object associated with this 
        ionization reaction

    """
    
    ionization_type = reaction_dict["superclass"]
    products_to_add = reaction_dict["products_to_add"]
    products_to_add.append(ionization_type)
    reaction_order = reaction_dict["reaction_order"]
    
    kinetiscope_name = (
        write_ionization_name(HiPRGen_reaction, reaction_writing_data, reaction_dict, ionization_type, products_to_add)
    )
    
    rate_constant = (
        find_ionization_rate_constant(reaction_writing_data, reaction_dict, ionization_type)
    )
    
    ionization_reaction = (
    build_rxn_object(HiPRGen_reaction, kinetiscope_name, rate_constant, reaction_order, products_to_add)
    )
    
    return ionization_reaction

def add_products_to_list(HiPRGen_reaction, product_list):
    classifications_without_tag = HiPRGen_reaction.classification_list[:-1]
    kinetiscope_reaction_products = classifications_without_tag
    
    for product in product_list:
        kinetiscope_reaction_products.append(product)
    
    return kinetiscope_reaction_products

def build_EI_reaction_dict(HiPRGen_reaction, reactant_electron, product_electron_list):
    superclass = "electron_ionization"
    product_electron_list.append(superclass)
    full_product_list = add_products_to_list(HiPRGen_reaction, product_electron_list)
    
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
        
        EI_reaction = build_one_EI_reaction(HiPRGen_reaction, reactant_electron, product_electron_list, reaction_writing_data)
        all_EI_reactions.append(EI_reaction)
        
    return all_EI_reactions

def build_absorption_reaction_dict(HiPRGen_reaction):
    superclass = "absorption"
    photoelectron = "eV_80" #photoelectron has 80 eV kinetic energy
    product_list = [photoelectron, superclass]
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
    absorption_reaction : Kinetiscope_Reaction object
        a reaction representing absorption of a photon leading to the emission
        of a photoelectron.

    """
    
    # superclass = "absorption"
    # photoelectron = "eV_80"
    # product_list = [superclass, photoelectron]
    # products_to_add = add_products_to_reaction(HiPRGen_reaction, product_list)
    # # kinetiscope_reaction_products = add_photoelectron_as_product(HiPRGen_reaction)
    # absorber = HiPRGen_reaction.reactants[0]
    
    
    # reaction_order = 1
    
    # reaction_dict = (
    # build_reaction_parameter_dict(absorber, "absorption", reaction_order,products_to_add=kinetiscope_reaction_products)
    # )
    reaction_dict = build_absorption_reaction_dict(HiPRGen_reaction)
    # reaction_dict = {
    #     "rate_constant_key":absorber,
    #     "superclass":superclass,
    #     "reaction_order":1,
    #     "reactants_to_add":None
    #     "products_to_add":products_to_add
    # }
    
    absorption_reaction = (
        build_ionization_reaction(HiPRGen_reaction, reaction_writing_data, reaction_dict)
    )
    
    return absorption_reaction  

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
    
    reactant_is_absorber = (
        species_is_absorber(potential_absorber, reaction_writing_data)
    )
    
    if reactant_is_absorber:
        
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