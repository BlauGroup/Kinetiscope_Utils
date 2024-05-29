# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:54:29 2024

@author: JRMilton
"""
import os
import sys
from monty.serialization import loadfn, dumpfn
from Rxn_classes import HiPRGen_Reaction

def find_charge(mpculeid):
    """
    This
    function simply returns the charge of a given species as an interger based
    on its mpculeid.

    Parameters
    ----------
    mpculeid : string
        a string of the form: "graph_hash-formula-charge-spin" 

    Returns
    -------
    charge : int
        charge of the species

    """
    charge_str = mpculeid.split("-")[2]
    if "m" in charge_str: #m stands for minus in the string
        charge = -int(charge_str.replace("m", ""))
        
    else:
        charge = int(charge_str)
    return charge
    
def is_ionization_reaction(reaction):
    """
    "Phase 1" of our EUV exposure model contains ionization reactions that
    we want to tag later when we're building our Kinetiscope reaction objects.
    This function tags these reactions.
    

    Parameters
    ----------
    reaction : HiPRGen_Reaction object
        the reaction we're testing

    Returns
    -------
    bool
        True if reaction is ionization, False otherwise.

    """
    num_reactants = len(reaction.reactants)
    num_products = len(reaction.products)
    if num_reactants+num_products == 2:
        reactant_mpculeid = reaction.reactants[0]
        product_mpculeid = reaction.products[0]
        reactant_charge = find_charge(reactant_mpculeid)
        product_charge = find_charge(product_mpculeid)
        if product_charge != reactant_charge:    
            return True
    return False

#tag P1 rxns

os.chdir("G:/My Drive/CRNs/041323_p1")
P1_rxn_tally = loadfn("reaction_tally.json")
P1_rxn_hashes = set()
HiPRGen_rxn_list = []

for rxn_dict in P1_rxn_tally["reactions"].values():
    new_rxn = HiPRGen_Reaction(rxn_dict,phase=1)
    if is_ionization_reaction(new_rxn):
        new_rxn.tag = "ionization"
        HiPRGen_rxn_list.append(new_rxn)
    P1_rxn_hashes.add(new_rxn.reaction_hash)

#read in HiPRGen rxns we want to write names for

os.chdir("G:/My Drive/CRNs/Kinetiscope_rxn_naming052824")
all_rxn_dicts = loadfn("euvl_TSreactions_041823.json")


for rxn_dict in all_rxn_dicts:
    new_rxn = HiPRGen_Reaction(rxn_dict,tag="chemical")
    new_rxn.phase = 1 if new_rxn.reaction_hash in P1_rxn_hashes else 2
    HiPRGen_rxn_list.append(new_rxn)

# Initialize the dictionary to hold the reactions grouped by their tags
tagged_rxn_dict = {}

for reaction in HiPRGen_rxn_list:
    
    tag = reaction.tag
    
    if tag not in tagged_rxn_dict:
        
        tagged_rxn_dict[tag] = []
    
    tagged_rxn_dict[tag].append(reaction)
 

# Now tagged_rxn_dict contains the grouped reactions by their tags

dumpfn(tagged_rxn_dict,"HiPRGen_rxns_to_name.json")