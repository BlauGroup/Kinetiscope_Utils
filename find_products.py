# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 13:11:32 2024

@author: jacob
"""

import os
from monty.serialization import loadfn

def remove_counter_and_time(line):
    line_list = line.split()
    return line_list[2:]

def read_lines_and_process(file_path, line_numbers):
    with open(file_path) as f:
        file_list = f.readlines()
    
    result = []
    for line_number in line_numbers:
        line_string = file_list[line_number - 1]  # Adjusting for zero-based indexing
        result.append(remove_counter_and_time(line_string))
    
    return result

os.chdir("G:\My Drive\Kinetiscope\Allabsorbers_020524")

name_list, initial_conc_list, final_conc_list = read_lines_and_process("duplicatesremoved_full_exposure.txt", [9, 12, 980])

final_conc_list = [float(value) for value in final_conc_list] #helps with sorting later

name_info_dict = {}
counter_species = ["IP", "eV80", "eV55", "eV30", "LEE"] #pseudospecies that are not real products

for i, name in enumerate(name_list):
    initial_conc = initial_conc_list[i]
    final_conc = final_conc_list[i]
    if name not in counter_species:
        name_info_dict[name] = {"initial_conc": initial_conc, "final_conc": final_conc}

formed_non_reactants = {name: info for name, info in name_info_dict.items() if float(info["initial_conc"]) == 0.0 and float(info["final_conc"] > 0.0)}
sorted_products = sorted(formed_non_reactants.items(), key=lambda x: x[1]["final_conc"], reverse=True)

top_10_products = [key for key, _ in sorted_products[:10]]
name_index_dict = loadfn("name_index_key.json")

index_list = []

for product in top_10_products:
    if name_index_dict.get(product, False):
        index = name_index_dict[product]
        index_list.append(index)

print(f"list of names of top 10 products: {top_10_products}")
print(f"list of indicies of top 10 products: {index_list}")
print(f"total number of products with nonzero concentration: {len(sorted_products)}")