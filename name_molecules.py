# -*- coding: utf-8 -*-
"""
Created on Wed May 10 13:26:19 2023

@author: JRMilton
"""

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
# from mc_analysis import ReportGenerator
import glob
# import pickle
from monty.serialization import loadfn, dumpfn
import os
import copy
import networkx as nx
import time
import csv
# import sqlite3
import copy
import operator
import sys

#TODO investigate why TPS_2+ is just TPS_2

"""
This module saves a json associating each molecule involved in relevent reactions
with a unique name based on functional group present within that molecule. It also
writes an excel file for the desired reactions of the form name1 + name2 => name 3
+ name 4, calculates the relative concentrations of species from a given trajectory,
and outputs those concentrations as text files. 
"""

def validate_species_to_name(species_dict):
    """
    Determines whether or not we have the data we need to generate a name for 
    a given species. If we do, returns that data, otherwise, raises a KeyError.

    Parameters
    ----------
    species_dict : dictionary
        dictionary associated with a given species loaded from json

    Raises
    ------
    KeyError
        raises a key error if this dictionary is missing graph or molecule
        data

    Returns
    -------
    species_graph : networkx graph
        networkx undirected graph associated with this species
    species_pymatgen_mol : pymatgen Molecule object
        the Molecule object associated with this species.

    """
    species_node_link_data = species_dict.get("nx_graph", False)
    species_graph = nx.node_link_graph(species_node_link_data)
    species_pymatgen_mol = species_dict.get("molecule", False)
    
    if not species_graph or not species_pymatgen_mol:
        raise KeyError("Missing required keys 'nx_graph' or 'molecule' in species_to_name.")
    
    return species_graph, species_pymatgen_mol

def update_species_name(species_name, func_group_name, func_name_already_added):
    """
    Adds the name of the functional group present within the molecule to its
    species name, accounting for the name of the functional group already 
    having been added to the name. 

    Parameters
    ----------
    species_name : string
        current name generated for the species
    func_group_name : string
        name of the functional group being added to the species name
    func_name_already_added : Boolean
        True if the name of this functional group is already in the species name,
        False otherwise.

    Returns
    -------
    species_name : string
        the updated name of this species

    """
    if func_group_name not in species_name:
        species_name += func_group_name + '_'
    else:
        if func_name_already_added:
            current_number_present = int(species_name[-2])
            new_number = current_number_present + 1
            species_name = species_name[:-2] + str(new_number) + '_'
        else:
            species_name = species_name[:-1] + '2_'
    
    return species_name

def update_remaining_species_graph(remaining_species_graph, func_group_loc):
    """
    

    Parameters
    ----------
    remaining_species_graph : networkx undirected graph
        The graph of a given species after the nodes associated with a
        functional group have been removed.
    func_group_locations : generator
        generator over isomorphisms between a subgraph of two graphs

    Returns
    -------
    remaining_species_graph : networkx undirected graph
        the graph of the species with the functional group removed.

    """
    for atom_index in func_group_loc.keys():
        remaining_species_graph.remove_node(atom_index)
    
    return remaining_species_graph

def update_species_name_with_atom_composition(species_name, name, remaining_species_graph):
    atom_composition = ""
    for node in remaining_species_graph.nodes(data=True):
        element = node[1]['specie']
        if element in atom_composition:
            element_index = atom_composition.index(element)
            old_number = int(atom_composition[element_index + 1])
            new_number = old_number + 1
            atom_composition = atom_composition[:-1] + str(new_number)
        else:
            atom_composition += element + "1"
    species_name += name + atom_composition + '_'
    
    return species_name

def generate_species_name(species_dict, func_group_dict):
    species_graph, species_pymatgen_mol = validate_species_to_name(species_dict)
    remaining_species_graph = copy.deepcopy(species_graph)
    species_name = ""

    for func_group_name, func_group_graph in func_group_dict.items():
        func_group_loc = functional_group_present(remaining_species_graph, func_group_graph)
        while func_group_loc:
            if species_name:
                func_name_already_added = species_name[-2].isnumeric()
            else:
                func_name_already_added = False
            species_name = update_species_name(species_name, func_group_name, func_name_already_added)
            remaining_species_graph = update_remaining_species_graph(remaining_species_graph, func_group_loc)
            func_group_loc = functional_group_present(remaining_species_graph, func_group_graph)
    
    species_name = update_species_name_with_atom_composition(species_name, name, remaining_species_graph) 
        
    charge_suffix = "+" if species_pymatgen_mol.charge == 1 else str(species_pymatgen_mol.charge)
    species_name += charge_suffix

    return species_name

def functional_group_present(mol_graph, func):
    """
    Tests whether or not a given functional group is present in a molecule via
    testing if the Networkx graphical representation of a molecule contains a 
    subgraph that is isomorphic to the functional group.
    
    Parameters
    ----------
    mol_graph : Networkx undirected graph
        Graph of the molecule you want to name
    func : Networkx undirected graph
        Graph of the functional group

    Returns
    -------
    Boolean
        True if a subgraph of the molecule is isomorphic to the functional group
        graph, False otherwise.
    Generator over isomorphisms between a subgraph of G1 and G2.
        This is a generator of mappings between a subgraph of the molecule and
        the functional group--this helps us deterimine where the functional
        group is in the original molecule

    """
    nm = nx.isomorphism.categorical_node_match("specie", None) #ensures isomorphic graphs must have the same atoms
    isomorphism = nx.isomorphism.GraphMatcher(mol_graph, func, node_match = nm)
    return isomorphism.mapping

def stereoisomer_test(stereoisomer_list, name):
    """
    Structural isomers--species with the same formula but different graphs--will return
    the same name when the name_molecule function is called. This function
    numbers each stereoisomer such that the new names are unique.

    Parameters
    ----------
    name_dict : dictionary
        dictionary whose keys are names and values are entry_ids of the species
        associated with that name
    mol_entry_id : string
        the unique entry id associated with a species
    mol_name : string
        the name generated for a species from the name_molecule function, which
        may not yet be unique

    Returns
    -------
    stereos : dictionary
        a dictionary of length 1 or 2, whose keys are now unique names of stereoisomers
        and values are entryids

    """
    new_stereos_list = []
    if len(stereoisomer_list) == 1:
        name_1 = name + '_#1'
        name_2 = name + '_#2'
        new_stereos_list.append(name_1)
        new_stereos_list.append(name_2)
    else:
        current_max_num = str(len(stereoisomer_list))
        new_stereos_list = stereoisomer_list
        new_isomer = name + '#_' + current_max_num
        new_stereos_list.append(new_isomer)
    return new_stereos_list

def update_names(test_name, stereo_dict, name_mpcule_dict, mpculeid):
    current_stereos = stereo_dict.get(test_name)
    new_stereos = stereoisomer_test(current_stereos, test_name)
    stereo_dict[test_name] = new_stereos
    
    old_stereo_mpculeid = name_mpcule_dict.get(test_name, False)
    if old_stereo_mpculeid:
        name_mpcule_dict[old_stereo_mpculeid] = new_stereos[-2]
    name_mpcule_dict[mpculeid] = new_stereos[-1]
    

# def write_reaction(reaction_dict, mpculid_dict): #convert to strings of the appropriate format for kinetiscope
#     """
#     Kinetiscope can read in reactions of the format name1 + name2 => name3
#     + name4. This function uses the names we've generated to write strings for
#     reactions of this format.

#     Parameters
#     ----------
#     reaction_dict : dictionary
#         dictionary whose keys are "reactants" and "products" and whose values
#         are a list of mpculeids associated with a given reaction
#     mol_entry_id : dictionary
#         dictionary associating each mpculeid with its name for fast lookup

#     Returns
#     -------
#     rxn_eqn: string
#         the desired reaction written in the form: name1 + name2 => name3
#         + name4
#     """
#     rxn_eqn = ""
#     if len(reaction_dict["reactants"]) == 1:
#         name = mpculid_dict[reaction_dict["reactants"][0]]
#         rxn_eqn = name + " "
#     else:
#         first_reactant = reaction_dict["reactants"][0]
#         first_name = mpculid_dict[first_reactant]
#         rxn_eqn = first_name + " + "
#         second_reactant = reaction_dict["reactants"][1]
#         second_name = mpculid_dict[second_reactant]
#         if first_name == second_name:
#             # print(first_name)
#             rxn_eqn = "2 " + first_name + " " #if reaction is A + A -> B + C we have to write 2A rather than A + A
#             # print(rxn_eqn)
#         else:
#             rxn_eqn = rxn_eqn + second_name + " "
#     rxn_eqn = rxn_eqn + "=> "
#     if len(reaction_dict["products"]) == 1:
#         name = mpculid_dict[reaction_dict["products"][0]]
#         rxn_eqn = rxn_eqn + name
#     else:
#          first_product = reaction_dict["products"][0]
#          first_name = mpculid_dict[first_product]
#          # rxn_eqn = rxn_eqn + first_name + " + "
#          second_product = reaction_dict["products"][1]
#          second_name = mpculid_dict[second_product]
#          if first_name == second_name:
#              # print(first_name)
#              rxn_eqn += "2 " + first_name #if reaction is A + A -> B + C we have to write 2A rather than A + A
#              # print(rxn_eqn)
#          else:
#              rxn_eqn = rxn_eqn + first_name + " + " + second_name
#     return rxn_eqn
    
# def generate_name_key(network_loader, key_dict, species_report_path):
#     """
#     print all species
#     """
#     report_generator = ReportGenerator(
#         mol_entries,
#         species_report_path,
#         rebuild_mol_pictures=False)

#     report_generator.emit_text("name key")
#     for name, entry_id in key_dict.items():
#         for mol_entry in mol_entries:
#             m_id = mol_entry.entry_id
#             if m_id == entry_id:
#                 report_generator.emit_text(str(entry_id))
#                 report_generator.emit_text(
#                     "formula: " + mol_entry.formula)
#                 report_generator.emit_molecule(mol_entry.index)
#                 report_generator.emit_newline()

#     report_generator.finished()
 
# def export_species_report_to_json(network_loader, path, key_dict):
#     data = {}
#     for i, entry_id in enumerate(key_dict.values()):
#         data[i] = entry_id

#     dumpfn(data, path)

# Change directory to the functional groups folder
func_groups_dir = r"G:\My Drive\Kinetiscope\import_test_021424\func_groups"
os.chdir(func_groups_dir)

# Process functional groups' XYZ files
print('Associating functional groups with their Molecule objects...')
func_group_dict = {}

for filename in glob.glob('*.xyz'):
    func_group = Molecule.from_file(filename)  # Load functional group as a Molecule
    func_group_mol_graph = MoleculeGraph.with_local_env_strategy(func_group, OpenBabelNN(order=False))
    func_group_undirected_graph = nx.Graph(func_group_mol_graph.graph)
    name = filename.replace('.xyz', '')
    func_group_dict[name] = func_group_undirected_graph

print('Functional group association completed.')

# Load molecule data and generate species names
os.chdir(r"G:\My Drive\Kinetiscope\import_test_021424")
mpcule_id_molecule_association_file = "mpcule_id_molecule_association.json"
mpcule_id_molecule_dict = loadfn(mpcule_id_molecule_association_file)

mpcule_name_dict = {}
name_mpcule_dict = {}
stereo_dict = {}

# Process each molecule and generate a species name
for mpcule_id, species_dict in mpcule_id_molecule_dict.items():
    species_name = generate_species_name(species_dict, func_group_dict)
    
    # Check and update stereo_dict if necessary
    if species_name in stereo_dict:
        update_names(species_name, stereo_dict, name_mpcule_dict, mpcule_id)
    else:
        mpcule_name_dict[mpcule_id] = species_name
        name_mpcule_dict[species_name] = mpcule_id
        stereo_dict[species_name] = [species_name]

print('Species names generated.')

# Finalize and output results
reactions_added = set()
reactions = []
# print(mpcule_name_dict)
# print('Naming Molecules...')

# # reactions_added = set()

# print('Adding reactions for phase 2 network products...')                                             
# second_name = 'sink_report'
# second_entries = loadfn(second_name + ".json") 

# for dictionary in second_entries.values(): #iterates through each network product
#     species_index = dictionary["species_index"]
#     # print(species_index)
#     reaction_json = loadfn(str(species_index) + "_pathway.json") #loads each network_product.json file as dictionary containing  
#     pathways = reaction_json["pathways"]                         #two keys: pathways and reactions.
#     pathways.sort(key = operator.itemgetter('frequency'), reverse = True) #sorts by frequency
#     top_pathways = []
#     n = 1
#     ten_saved = False
    
#     while not ten_saved: #makes sure we save top ten reactions even for network products that don't have 10 1-step reactions forming them
#         for dictionary in pathways:
#             if int(dictionary['weight']) == n:
#                 top_pathways.append(dictionary['pathway']) #top_pathways will be a list of ints, with each number corresponding to a reaaction index
#                 if len(top_pathways) == 10:
#                     ten_saved = True
#                     break
#         n += 1
    
#     for pathway in top_pathways:
#         for reaction in pathway:
#             if reaction not in reactions_added:
#                 reactions_added.add(reaction)
#                 rxn = reaction_json["reactions"][str(reaction)]
#                 reactions.append(rxn)
#                 molecule_dict = find_molecules_from_mpculeids(rxn, mol_entries) #want to associate an mpculeid with a name
#                 for l in molecule_dict.values():
#                     for molecule in l:
#                         if molecule.entry_id not in entry_ids: #i.e. we haven't named this molecule yet
#                             name = name_molecule(molecule, func_group_dict)
#                             if stereoisomer_test(name_dict, molecule.entry_id, name): #returns empty dict if no stereoisomers already in dict
#                                 name_dict.update(stereoisomer_test(name_dict, molecule.entry_id, name))
#                                 if len(stereoisomer_test(name_dict, molecule.entry_id, name)) == 2: #If two entries are returned, this means we've named two stereoisomers #1 and #2, and need to remove the duplicate name from the dict
#                                     name_dict.pop(name)
#                             else:
#                                 name_dict[name] = molecule.entry_id
#                             entry_ids.add(molecule.entry_id)
#                 # added.append(num)
#                 # added_hashes.update(reaction_dict.items())
#             # reactants = reaction_json["reactions"][num]["reactants"]
#             # products = reaction_json["reactions"][num]["products"]
#             # participants = [reactants, products] 
#             # reaction_dict = create_reaction_dict(participants)
#             #         if not resonant_reaction(reaction_dict, added_hashes):
#             #             if not reverse_reaction(reaction_dict, added_hashes):
#             #                 if not charge_transfer_reaction(reaction_dict):
                                
#     # n = 1
                
# # print('Done! ', len(mpcule_ids), ' reactions total')

# print('Done!', len(reactions_added), 'added forming network products!')

# print('Opening json...')  
# third_name = 'reaction_tally'
# third_entries = loadfn(third_name + ".json")
# print('Done!') 

# print("Naming reactions that fired more than 500 times...")

# for reaction in third_entries["pathways"].keys():
#     if third_entries["pathways"][reaction] > 500: #only add network products found >500 times
#         if reaction not in reactions_added:
#             reactions_added.add(reaction)
#             rxn = third_entries["reactions"][reaction]
#             reactions.append(rxn)
#             molecule_dict = find_molecules_from_mpculeids(rxn, mol_entries) #want to associate an mpculeid with a name
#             for l in molecule_dict.values():
#                 for molecule in l:
#                     if molecule.entry_id not in entry_ids: #i.e. we haven't named this molecule yet
#                         name = name_molecule(molecule, func_group_dict) 
#                         if stereoisomer_test(name_dict, molecule.entry_id, name): #returns empty dict if no stereoisomers already in dict
#                             name_dict.update(stereoisomer_test(name_dict, molecule.entry_id, name))
#                             if len(stereoisomer_test(name_dict, molecule.entry_id, name)) == 2: #If two entries are returned, this means we've named two stereoisomers #1 and #2, and need to remove the duplicate name from the dict
#                                 name_dict.pop(name)
#                         else:
#                             name_dict[name] = molecule.entry_id
#                         entry_ids.add(molecule.entry_id)
                                
# print('Done!', len(reactions_added), 'added!')
# # print('Done naming molecules!')
# # end = time.time()
# # total = end - start
# # time_min = total/60
# # time_min = round(time_min, 2)
# # print('named ', len(name_dict), 'molecules and took', time_min, ' minutes')

# print("Done naming molecules!")

# n = 0
# longest = 0

# for name in name_dict.keys():
#     if len(name) > 32:
#         n += 1
#         print(name)
#         if len(name) > longest:
#             longest = len(name)
# if longest > 0:  
#     print('Warning, ', n, 'too long names have been generated!' )
#     print('Longest length: ', longest)
    
# dict_set = set(name_dict.values())
# if entry_ids.difference(dict_set):
#     print('Unnamed molecule hashes: ', entry_ids.difference(dict_set))
        
# l_list = []
    
# print('Testing for duplicates...')
# for name, h in name_dict.items(): 
#     for na, ha in name_dict.items():
#         if na == name and ha != h:
#             l_list.append((ha, h))

# if l_list:
#     print('# of duplicate entries found: ', len(l_list))
#     print(l_list)
    
# print('Done!') 

# # name_dict = loadfn('named_molecules.json')
# mpculid_dict = dict([(value, key) for key, value in name_dict.items()]) #wanted to do this while generating the names, but can't because the names can change as we're building the dictionary

# key_dict = {} #consider printing a key instead if it seems necessary
# for name, mpculeid in name_dict.items():
#     for species in mol_entries:
#         if mpculeid == species.entry_id:
#             key_dict[name] = species.ind

# if len(key_dict) != len(name_dict):
#     print('Error: not all names associated with a molecule index')
#     for name in name_dict.keys():
#         if not key_dict.get(name, False):
#             print('Name not in key: ', name)
# else:
#     print('Saving key json...')
#     os.chdir(new_dir)
#     dumpfn(key_dict, 'name_index_key.json')
#     os.chdir(original_directory)
#     print('Done!')

# index_dict = dict([(value, key) for key, value in key_dict.items()]) #want to look up name by index later
# named_reactions = []

# print('Converting reactions containing names...')
# # for reaction in third_entries["pathways"].keys():
# #     if third_entries["pathways"][reaction] > 500: #only add network products found >500 times
# #         rxn = third_entries["reactions"][reaction]
# for rxn in reactions:
#         # for rxn in third_entries["reactions"].keys():
#             # if str(reaction) == rxn:
#     named_reaction = write_reaction(rxn, mpculid_dict) #iterate through the desired reactions, and find the name from the mpcule id
#     named_reactions.append(named_reaction)

# print('Done!')

# print('Writing reactions to csv file...')
# # with open('Kinetiscope_rxn_template.csv', newline = "") as csvfile:
# #     # reader = csv.reader(csvfile)
# #     # fields = list(next(reader))

# dict_list = []

# for reaction in named_reactions:
#     csv_dict = {}
#     csv_dict['# equation'] = reaction
#     csv_dict['fwd_A'] = 1
#     csv_dict['fwd_temp_coeff'] = 0
#     csv_dict['fwd_Ea'] = 0
#     csv_dict['fwd_k'] = 10379761818429.5 #kT/h
#     csv_dict['rev_A'] = 1
#     csv_dict['rev_temp_coeff'] = 0
#     csv_dict['rev_Ea'] = 0
#     csv_dict['rev_k'] = 1
#     csv_dict['fwd_k0'] = 1
#     csv_dict['rev_k0'] = 1
#     csv_dict['alpha_alv'] = 0.5
#     csv_dict['equil_potential'] = 0.5
#     csv_dict['num_electrons'] = 0
#     csv_dict['fwd_prog_k'] = 1
#     csv_dict['rev_prog_k'] = 1
#     csv_dict['non_stoichiometric'] = 0
#     csv_dict['rate_constant_format'] = 0
#     dict_list.append(csv_dict)

# os.chdir(new_dir)
    
# with open("euvl_phase2_reactions.csv", 'w', newline = "") as csvfile:
#     writer = csv.DictWriter(csvfile, fieldnames = dict_list[0].keys())
#     writer.writeheader()
#     for reaction in dict_list:
#         writer.writerow(reaction)
  
# print('Done!')     
#     # writer = csv.writer(csvfile)

# print('Writing concentration text files...')   
# os.chdir(original_directory)
 
# database = input("Please input the path of the sqlite3 database: ")

# initial_state_con = sqlite3.connect(database)
# cur = initial_state_con.cursor()
# concentration_dict = {}

# for row in cur.execute(sql_get_initial_state):
#       concentration_dict[row[0]] = row[1] #associates each species index with its particle count in that trajectory

# # # folder = "Concentration_files"
# # current_directory = os.getcwd()
# # path = os.join(current_directory, folder)
# # os.mkdir(path)
# # os.chdir(path)

# os.chdir(new_dir)

# total_num_particles = float(sum(concentration_dict.values()))
# numbers = "1 1 1 "
# for name in name_dict.keys():
#     filename = name + ".txt"
#     with open(filename, "w") as f:
#         f.write("# column(x) row(y) layer(z) value")
#         f.write('\n')
#         name_index = key_dict[name]
#         species_particle_number = float(concentration_dict[name_index])
#         ratio = species_particle_number / total_num_particles
#         concentration = ratio * 8.076
#         if concentration >= 1: #deals with significant figures
#             rounded_concentration = round(concentration, 3)
#         elif concentration < 0.01:
#             rounded_concentration = round(concentration, 5)
#         else:
#             rounded_concentration = round(concentration, 4)
#         to_write = numbers + str(rounded_concentration)
#         f.write(to_write)
    
# print('Done!') 

# end = time.time()
# total = end - start
# time_min = total/60
# time_min = round(time_min, 2)
# print('named ', len(name_dict), 'molecules and took', time_min, ' minutes')

        # initial_state_array = np.zeros(
        #     self.number_of_species,
        #     dtype=int
        # )

        # for i in range(self.number_of_species):
        #     initial_state_array[i] = initial_state_dict[i]

        # if self.initial_state_dict == {} and self.initial_state_array == {}:
        #     self.initial_state_dict = initial_state_dict
        #     self.initial_state_array = initial_state_array
        # else:
        #     for i in range(self.number_of_species):
        #         if initial_state_array[i] > self.initial_state_array[i]:
        #             self.initial_state_array[i] = initial_state_array[i]
        #             self.initial_state_dict[i] = initial_state_dict[i]
#save to the excel file
# dumpfn(name_dict, 'named_molecules.json')

# print('Done! ', len(mpcule_ids), ' reactions total')

# kinetiscope_reaction_list = []
