# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 16:37:32 2024

@author: jacob
"""

import csv
import os

def remove_duplicates(input_file, output_file):
    unique_rows = set()  # Use a set to store unique rows
    duplicates_removed = 0  # Counter for duplicate rows removed
    with open(input_file, 'r', newline='') as infile:
        reader = csv.reader(infile)
        for row in reader:
            # Convert each row to a tuple to make it hashable
            row_tuple = tuple(row)
            if row_tuple not in unique_rows:
                unique_rows.add(row_tuple)
            else:
                duplicates_removed += 1

    # Write unique rows to the output file
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        for row in unique_rows:
            writer.writerow(row)

    print(f"Number of duplicate rows removed: {duplicates_removed}")

os.chdir("G:\My Drive\Kinetiscope\Allabsorbers_020524")

# Example usage
input_file = 'full_absorbers_021624.csv'
output_file = 'full_absorbers_duplicateremoved_02204.csv'
remove_duplicates(input_file, output_file)
