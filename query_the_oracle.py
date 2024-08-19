# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 14:42:27 2024

@author: jacob
"""
import os
from monty.serialization import loadfn, dumpfn

os.chdir("G:/My Drive/Kinetiscope/new_kinetiscope_naming_080224")
chemical_dict = loadfn("name_full_mpculeid_corrected_081324.json")

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