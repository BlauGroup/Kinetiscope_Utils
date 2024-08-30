# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 14:24:02 2024

@author: jacob
"""

from write_positive_ionization_reactions import build_positive_ionization_reactions
from write_attachment_and_recombination import build_attachment_or_recombination

def select_ionization_builder(HiPRGen_reaction, reaction_writing_data):
    """
    We handle naming positive ionization reactions differently from attachment
    and recombination reactions. This function just calls the appropriate script
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