# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 12:22:03 2024

@author: jacob
"""

from write_phase_1_chemical_reactions import write_phase_1_chemical_reactions
from write_phase_2_chemical_reactions import write_phase_2_chemical_reactions

def select_chemical_builder(HiPRGen_rxn, reaction_writing_data):
    
    if HiPRGen_rxn.phase == 1:
        
        return write_phase_1_chemical_reactions(HiPRGen_rxn, reaction_writing_data)
    
    return write_phase_2_chemical_reactions(HiPRGen_rxn, reaction_writing_data)
    