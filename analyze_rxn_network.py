# -*- coding: utf-8 -*-
import networkx as nx
import os

def find_edge_info(rxn_name):
    reactant, product = rxn_name.split(" => ")
    reactant_list = reactant.split(" + ")
    product_list = product.split(" + ")
    return reactant_list, product_list

def write_edge_list(species_node, side_id, reaction):
    if side_id == "reactant":
        return [species_node, reaction["name"]]
    else:
        return [reaction["name"], species_node]

#hypothesis: our network can be represented by a cyclic, undirected graph.
#reaction pathways can be found by searching for trees contained within <-this appears to be incorrect
#this graph.

#represent network

#find the minimum spanning tree connecting product to a given reactant
#check to see if this is a the pathway forming that product

network = nx.MultiGraph()

# name_file_path = "G:\My Drive\Kinetiscope\Allabsorbers_020524"
# os.chdir(name_file_path)
# with open("duplicatesremoved_full_exposure.txt", "r") as f:
#     lines = f.readlines()
#     name_nodes = lines[8].split()[2:]


reactant_nodes = ["A", "B", "C", "D"]
product_nodes = ["G"]
intermediate_nodes = ["E", "F"]
all_nodes = {"reactants": reactant_nodes, "products": product_nodes, 
             "intermediates": intermediate_nodes}

rxn_edges = {
    "A + B => E + F": {"name": "A + B => E + F", "select_freq": 1/1234},
    "C + D => E + F": {"name": "C + D => E + F", "select_freq": 1/123},
    "E + F => G": {"name": "E + F => G", "select_freq": 1/1200}
}

for tag, node_list in all_nodes.items():
    network.add_nodes_from(node_list, tag=tag)

# network.add_nodes_from([rxn["name"] for rxn in rxn_nodes.values()], bipartite=1)

for rxn in rxn_edges.values():
    reactant_list, product_list = find_edge_info(rxn["name"])
    freq = rxn["select_freq"]
    for reactant in reactant_list:
        for product in product_list:
        # start, end = write_edge_list(reactant, "reactant", rxn)
            network.add_edge(reactant, product, weight=freq)
    # for product in product_list:
    #     start, end = write_edge_list(product, "product", rxn)
    #     network.add_edge(start, end, weight=freq)

print(nx.minimum_spanning_tree(network, weight="freq")) #test for cycles
print(sorted(nx.minimum_spanning_tree(network, weight="freq").edges(data=True)))
# for reactant in reactant_nodes:
#     for product in product_nodes:
#         paths = nx.all_simple_paths(network, reactant, product)