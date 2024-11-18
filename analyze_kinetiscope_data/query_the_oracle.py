# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 14:42:27 2024

@author: jacob
"""
import os
from monty.serialization import loadfn, dumpfn

os.chdir("G:/My Drive/Kinetiscope/production_simulations_092124")
chemical_dict = loadfn("name_full_mpculeid_092124.json")

while True:
    # Query user for a chemical name
    chemical_name = input("Enter a chemical name (or 'end' to stop): ")
    
    # Check if the user wants to end the loop
    if chemical_name.lower() == 'end':
        print("Exiting the script.")
        break
    
    # Check if the chemical name is in the dictionary
    if chemical_name in chemical_dict:
        print(f"The mpculeid for {chemical_name} is {chemical_dict[chemical_name]}.")
    else:
        print(f"Error: {chemical_name} not found in the dictionary.")

# this intakes mpculeid and outputs name

# mpculeid_to_name = {value: key for key, value in chemical_dict.items()}

# while True:
#     # Query user for an mpculeid
#     mpculeid = input("Enter an mpculeid (or 'end' to stop): ")
    
#     # Check if the user wants to end the loop
#     if mpculeid.lower() == 'end':
#         print("Exiting the script.")
#         break
    
#     # Check if the mpculeid is in the reversed dictionary
#     if mpculeid in mpculeid_to_name:
#         print(f"The chemical name for mpculeid {mpculeid} is {mpculeid_to_name[mpculeid]}.")
#     else:
#         print(f"Error: mpculeid {mpculeid} not found in the dictionary.")