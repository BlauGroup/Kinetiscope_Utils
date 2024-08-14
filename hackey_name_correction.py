# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 20:51:29 2024

@author: jacob
"""
import os
from monty.serialization import loadfn, dumpfn

new_dict = {}
os.chdir("C:/Users/jacob/Kinetiscope_Utils")
to_fix = loadfn("name_full_mpculeid_080624.json")

count = 0

for name, mpculeid in to_fix.items():
    if len(name) > 32:
        count +=1
        new_name = name.replace("phenyl", "p")
        new_dict[new_name] = mpculeid
    else:
        new_dict[name] = mpculeid

for name in new_dict:
    assert len(name) <= 32

print(f"{count} names successfully changed")
dumpfn(new_dict, "name_full_mpculeid_corrected_081324.json")