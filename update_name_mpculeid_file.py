# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 16:28:56 2024

@author: jacob
"""
from name_molecules import build_species_graph
from networkx import weisfeiler_lehman_graph_hash
import os
from monty.serialization import loadfn, dumpfn

def update_mpculeid(old_mpculeid, graph_data):
    mpculeid_as_list = old_mpculeid.split('-')
    species_graph = build_species_graph(graph_data)
    new_graph_hash = weisfeiler_lehman_graph_hash(species_graph, node_attr="specie")
    mpculeid_as_list[0] = new_graph_hash
    return "-".join(mpculeid_as_list)

os.chdir("G:/My Drive/CRNs/fixing_mpculeids_080524")
mpculeid_molecule_dict = loadfn("mpcule_id_molecule_association.json")
name_mpculeid_dict = loadfn("name_mpculeid_association_060424.json")
updated_name_mpculeid_file = {}

for name, old_mpculeid in name_mpculeid_dict.items():
    associated_graph_data = \
        mpculeid_molecule_dict.get(old_mpculeid, None)["nx_graph"]
    new_mpculeid = update_mpculeid(old_mpculeid, associated_graph_data)
    updated_name_mpculeid_file[name] = new_mpculeid

dumpfn(updated_name_mpculeid_file, "name_mpculeid_association_updated_mpculeids_080524.json")