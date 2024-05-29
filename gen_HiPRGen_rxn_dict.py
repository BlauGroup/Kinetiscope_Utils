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

def determine_ionization_type(reactant_charge, product_charge):
    delta_charge = product_charge-reactant_charge
    if delta_charge < 0:
        # if reactant_charge :
        return "electron_attachment"
        # else:
        #     return "charge_quenching"
    else:
        return "positive_ionization"
    
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
            return determine_ionization_type(reactant_charge, product_charge)
    return False

def is_H_rxn():
    pass

def is_ion_ion():
    pass

def is_aromatic_sub():
    pass

def is_ion_molecule():
    pass

def is_neutral_radical():
    pass

def tag_rxn():
    pass

#tag P1 rxns

os.chdir("G:/My Drive/CRNs/041323_p1")
P1_rxn_tally = loadfn("reaction_tally.json")
P1_rxn_hashes = set()
HiPRGen_rxn_list = []

for rxn_dict in P1_rxn_tally["reactions"].values():
    new_rxn = HiPRGen_Reaction(rxn_dict,phase=1)
    if is_ionization_reaction(new_rxn):
        new_rxn.tag = is_ionization_reaction(new_rxn)
        HiPRGen_rxn_list.append(new_rxn)
    P1_rxn_hashes.add(new_rxn.reaction_hash)


#read in HiPRGen rxns we want to write names for

os.chdir("G:/My Drive/CRNs/Kinetiscope_rxn_naming052824")
all_rxn_dicts = loadfn("euvl_TSreactions_041823.json")


for rxn_dict in all_rxn_dicts:
    new_rxn = HiPRGen_Reaction(rxn_dict,tag="chemical")
    new_rxn.phase = 1 if new_rxn.reaction_hash in P1_rxn_hashes else 2
    HiPRGen_rxn_list.append(new_rxn)

modifications = []
for reaction in HiPRGen_rxn_list:
    if reaction.tag == "electron_attachment":
        for index, rxn in enumerate(HiPRGen_rxn_list):
            if rxn.tag == "positive_ionization":
                if reaction.reactants == rxn.products and reaction.products == rxn.reactants:
                    # Change the tag to recombination
                    reaction.tag = "recombination"
                    # Collect the modification
                    modifications.append((index, {
                        "positive_ionization": rxn,
                        "recombination": reaction
                    }))

# Apply the modifications after iteration
for index, new_value in modifications:
    HiPRGen_rxn_list[index] = new_value
# Initialize the dictionary to hold the reactions grouped by their tags
tagged_rxn_dict = {}

for reaction in HiPRGen_rxn_list:
    
    # if reaction.tag == "electron_attachment":
    #     for index, rxn in enumerate(tagged_rxn_dict["ionization_and_recombination"]):
    #         if reaction.reactants == rxn.products and reaction.products == rxn.reactants:
    #             reaction.tag = "recombination"
    #             tagged_rxn_dict["ionization_and_recombination"][index] = {
    #                 "positive_ionization":rxn,
    #                 "recombination":reaction}
                
    if reaction.tag not in tagged_rxn_dict:
        
        tagged_rxn_dict[reaction.tag] = []
    
    tagged_rxn_dict[reaction.tag].append(reaction)
 
# for reaction in tagged_rxn_dict["electron_attachment"]:
#     for rxn in tagged_rxn_dict["positive_ionization"]:
#         if reaction.reactants == rxn.products and reaction.products == rxn.reactants:
#             reaction.tag = "recombination"

dumpfn(tagged_rxn_dict,"HiPRGen_rxns_to_name.json")