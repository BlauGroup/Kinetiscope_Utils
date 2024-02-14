# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 12:21:02 2024

@author: JRMilton
"""
import pickle
from monty.serialization import dumpfn
import networkx as nx

print('Loading mol_entries.pickle...')

with open('mol_entries.pickle', 'rb') as f: #loads pymatgen Molecule objects from pickle file

    mol_entries = pickle.load(f)
    
print('Done!')

mpcule_id_molecule_dict = {}

for entry in mol_entries:
    
    mpcule_id = entry.entry_id
    
    if mpcule_id:
        assert mpcule_id not in mpcule_id_molecule_dict
        undirected_graph = nx.Graph(entry.graph)
        json_graph = nx.node_link_data(undirected_graph)
        mpcule_id_molecule_dict[mpcule_id] = {"molecule" : entry.molecule, "nx_graph" : json_graph}
 
dumpfn(mpcule_id_molecule_dict, "mpcule_id_molecule_association.json")