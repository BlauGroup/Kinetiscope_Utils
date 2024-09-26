# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 13:54:20 2024

@author: jacob
"""

#written by chatgpt

import os

def save_chemical_concentration():
    while True:
        # Ask for the chemical name
        chemical_name = input("Enter the chemical name (or type 'end' to quit): ").strip()

        # Break the loop if the user types 'end'
        if chemical_name.lower() == "end":
            print("Ending the program.")
            break

        # Create the filename with the chemical name
        filename = f"{chemical_name}.txt"

        # Check if the file already exists
        if os.path.exists(filename):
            print(f"Error: A file named '{filename}' already exists.")
            continue

        # Ask for the concentration
        try:
            concentration = float(input(f"Enter the concentration for {chemical_name}: "))
        except ValueError:
            print("Error: Please enter a valid floating point number for the concentration.")
            continue

        # Write the concentration to the file
        with open(filename, 'w') as file:
            file.write(f"1 1 1 {concentration}\n")

        print(f"Saved concentration of {concentration} for {chemical_name} in '{filename}'.")

os.chdir("G:/My Drive/Kinetiscope/production_simulations_092124")
# Run the function
save_chemical_concentration()