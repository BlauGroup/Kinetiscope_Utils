# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 12:55:07 2024

@author: jacob
"""
import pickle
from monty.serialization import loadfn
import os
import sys
import networkx as nx
import pymongo
import copy
from collections import Counter
from pymatgen.core.structure import Molecule
sys.path.append('../analyze_kinetiscope_data')
from build_top_reaction_dict import generate_rxn_data_dict
sys.path.append('../common')
from utilities import correct_path_change_dir
from Rxn_classes import HiPRGen_Reaction, RecombMolEntry
sys.path.append('../Classify_HiPRGen_reactions')
from reaction_classification_utilities import (
    find_mpculeid_charge,
    find_mpculeid_spin
)
import copy

def get_nonzero_freq_reactions(select_freq_file, reaction_name_file, 
                           start_index, end_index):
    def generate_kscope_rxns(select_freq_file, reaction_name_file, 
                             start_index, end_index):
        return generate_rxn_data_dict(
            select_freq_file, reaction_name_file, start_index, end_index
        )[1]

    kscope_rxns = generate_kscope_rxns(
        select_freq_file, reaction_name_file, start_index, end_index
    )

    return {
        index: rxn for index, rxn in kscope_rxns.items()
        if getattr(rxn, 'selection_freq', 0) != 0
    }


def find_times_radical_formed(kscope_rxns, chemical_dict):
    relevant_radical_dict = {}
    
    for kscope_rxn in kscope_rxns.values():
        if kscope_rxn.selection_freq > 0:
            for species in kscope_rxn.products:
                mpculeid = chemical_dict.get(species, None)
                if mpculeid:
                    charge = find_mpculeid_charge(mpculeid)
                    spin = find_mpculeid_spin(mpculeid)
                    if charge == 0 and spin == 2:
                        if mpculeid not in relevant_radical_dict:
                            relevant_radical_dict[mpculeid] = kscope_rxn.selection_freq
                        else:
                            relevant_radical_dict[mpculeid] += kscope_rxn.selection_freq

    return relevant_radical_dict


def find_relevant_radicals(relevant_radical_dict, threshold):
    def find_total_radical_formation_reactions(relevant_radical_dict):
        total = 0

        for formation_freq in relevant_radical_dict.values():
            total += formation_freq
        return total

    relevant_radical_set = set()

    total = find_total_radical_formation_reactions(relevant_radical_dict)

    for mpculeid, formation_freq in relevant_radical_dict.items():
        test = formation_freq/total
        if test >= threshold:
            relevant_radical_set.add(mpculeid)
    return relevant_radical_set


def find_underbonded_sites(relevant_radical_set, new_to_old, mol_entries):
    # def find_mol_entries(relevant_radical_set, mol_entries):
    #     mol_entry_lookup = {}
        
    #     for mol_entry in mol_entries:
    #         if mol_entry.entry_id in relevant_radical_set:
    #             mol_entry_lookup[mol_entry.entry_id] = mol_entry
    #     return mol_entry_lookup

    underbonded_sites = {}
    
    
    for relevant_radical in relevant_radical_set:
        underbonded_atoms = []
        lookup_mpculeid = new_to_old[relevant_radical]
        document = collection.find_one({"molecule_id": lookup_mpculeid})
        spins = document["partial_spins"]["DIELECTRIC=3,00"]["nbo"]
        for atom_index, spin in enumerate(spins):
            species  = document["species"][atom_index]
            if spin > 0.09 and species != "H":
                underbonded_atoms.append(atom_index)
        underbonded_sites[relevant_radical] = underbonded_atoms
    return underbonded_sites

def find_mol_entries(relevant_radical_set, mol_entries):
    mol_entry_lookup = {}
    
    for mol_entry in mol_entries:
        if mol_entry.entry_id in relevant_radical_set:
            mol_entry_lookup[mol_entry.entry_id] = mol_entry
    return mol_entry_lookup

def collect_lists_from_nested_dict(d):
    """
    Collect all list items from a nested dictionary.

    This function traverses a nested dictionary and collects all items
    from lists contained within the dictionary. It recurses through nested
    dictionaries and aggregates items from all lists found.

    Parameters:
    -----------
    d : dict
        The nested dictionary from which lists are to be collected.

    Returns:
    --------
    list
        A list containing all items collected from lists within the nested
        dictionary.
    """
    def recurse_through_dict(d, collected_items):
        """
        Recursively traverse a nested dictionary to collect items from lists.

        This helper function is used by `collect_lists_from_nested_dict` to
        perform a recursive traversal of the dictionary. It extends the list
        of collected items with elements from any lists found in the dictionary.

        Parameters:
        -----------
        d : dict
            The nested dictionary to traverse.
        collected_items : list
            The list to which items from nested lists are added.
        """

        for value in d.values():

            if isinstance(value, dict):

                recurse_through_dict(value, collected_items)

            elif isinstance(value, list):

                collected_items.extend(value)

    collected_items = []
    recurse_through_dict(d, collected_items)

    return collected_items

def build_radicals_sets(mol_entries, new_to_old, relevant_radical_set):
    def find_starting_species_graphs(mol_entries):
        starting_species_mpculeids = {
            "TPS": "bfed458e642b8daa1eab6bc02d5e5682-C18H15S1-1-1",
            "Nf": "9a8a88b8b92c714d7f65b8526ffabc7a-C4F9O3S1-m1-1",
            "4-CNBZ": "17f31f89123edbaa0e3b9c7eb49d26f3-C8H4N1O2-m1-1",
            "PHS": "00a7dcc352b0d613f58e850935bf5609-C10H14O1-0-1",
            "PtBMA": "4864aee73a83d357c31fadd50b81e3cd-C10H20O2-0-1"
        }
        starting_species_graphs = []
        
        for mol_entry in mol_entries:
            lookup_mpculeid = new_to_old[mol_entry.entry_id]
            if lookup_mpculeid in starting_species_mpculeids.values():
                starting_species_graphs.append(nx.Graph(mol_entry.graph))
        return starting_species_graphs
    
    def is_neutral_radical(mol_entry):
        return mol_entry.spin_multiplicity == 2 and mol_entry.charge == 0
    
    def radical_is_relevant(mol_entry, relevant_radical_set):
        return mol_entry.entry_id in relevant_radical_set
    
    def radical_isnt_fragment(mol_entry, starting_species_graphs):
        nm = nx.isomorphism.categorical_node_match("specie", None) # ensures isomorphic graphs must have the same atoms
        for graph in starting_species_graphs:

            # ismorphism.subgraph_is_isomorphic returns True if a subgraph of
            # G1 is isomorphic to G2.

            isomorphism = nx.isomorphism.GraphMatcher(graph, nx.Graph(mol_entry.graph), node_match = nm)
            if isomorphism.subgraph_is_isomorphic():
                return False
        return True

    all_radicals = set()
    # all_pi_radicals = []
    potentially_problematic = set()

    starting_species_graphs = find_starting_species_graphs(mol_entries)

    for mol_entry in mol_entries:

        if is_neutral_radical(mol_entry) and radical_is_relevant(mol_entry, relevant_radical_set):

            all_radicals.add(mol_entry.entry_id)
            if radical_isnt_fragment(mol_entry, starting_species_graphs):
                potentially_problematic.add(mol_entry.entry_id)
    total_number_radicals = len(all_radicals)
    number_problematic_radicals = len(potentially_problematic)
    # nonproblematic = all_radicals - potentially_problematic
    # for mpculeid in nonproblematic:
    #     print(new_to_old[mpculeid])
    #     print()
    # sys.exit()
    # # for mpculeid in potentially_problematic:
    # #     print(new_to_old[mpculeid])
    # #     print()
    # # sys.exit()
    # print(f"total number of radicals: {total_number_radicals}")
    # print(f"total number of problematic radicals: {number_problematic_radicals}")
    # print(f"estimated number of reactions: {total_number_radicals * number_problematic_radicals}")
    # sys.exit()
    return all_radicals, potentially_problematic



def get_radical_node_mapping(radical_nodes, problem_radical_nodes):
    offset = len(problem_radical_nodes)
    
    return {index: index+offset for index in radical_nodes}

def collect_graph_and_nodes(radical_mpculeid, radical_mol_entries):
    def ensure_nodes_are_sequential_ints(int_list):
        """
        Checks if a list contains sequential integers.
        
        Parameters:
            int_list (list): A list of integers to check.
            
        Raises:
            ValueError: If the list does not contain sequential integers.
        """
        for i in range(1, len(int_list)):
            # Check if each element is an integer
            if not isinstance(int_list[i], int):
                raise ValueError(f"Element at index {i} is not an integer: {int_list[i]}")
            # Check if the current element is exactly 1 greater than the previous one
            if int_list[i] != int_list[i - 1] + 1:
                raise ValueError(f"List is not sequential at index {i}: {int_list[i-1]} -> {int_list[i]}")

    mol_entry = radical_mol_entries[radical_mpculeid]
    graph = mol_entry.graph
    graph_node_list = list(graph.nodes)
    ensure_nodes_are_sequential_ints(graph_node_list)
    return graph, graph_node_list

def build_combined_graph(problem_radical, radical, radical_mol_entries, underbonded_sites):
        problem_graph, problem_node_list = collect_graph_and_nodes(
            problem_radical, radical_mol_entries
        )
        radical_graph, radical_node_list = collect_graph_and_nodes(
            radical, radical_mol_entries
        )
        mapping = get_radical_node_mapping(radical_node_list, problem_node_list)
        
        combined_graph = nx.disjoint_union(problem_graph, radical_graph)
        nodes = list(combined_graph.nodes(data='specie'))
        return combined_graph, nodes, mapping


def write_new_formula(node_list):
    atom_counts = Counter(tup[1] for tup in node_list)

    # strings (which we sort here) are sorted alphabetically by default

    sorted_atoms = sorted(atom_counts.items())
    return ''.join(f"{atom}{count}" for atom, count in sorted_atoms)

def generate_recombinant_graph(combined_graph, problem_underbonded_atom, atom_in_combined_graph):
    graph_copy = copy.deepcopy(combined_graph)
    atom_in_combined_graph = mapping[original_atom]
    graph_copy.add_edge(
        problem_underbonded_atom, atom_in_combined_graph
    )
    
    return graph_copy
    
    
def WriteNewMpculeid(graph, formula, charge_int, spin):
    def GenerateGraphHash(graph):
        graph_hash = nx.weisfeiler_lehman_graph_hash(
           nx.Graph(graph), node_attr="specie"
        )
        return graph_hash

    def WriteChargeString(charge_int):
        charge_string = str(charge_int) if charge_int >= 0 else "m" + str(abs(charge_int))
        return charge_string
    
    
    
    graph_hash = GenerateGraphHash(graph)
    charge_string = WriteChargeString(charge_int)
    
    to_combine = [graph_hash, formula, charge_string, str(spin)]
    
    mpculeid = "-".join(to_combine)
    
    return mpculeid


def write_recombination_reaction(reactant_list, recombinant):
    return " + ".join(reactant_list) + " => " + recombinant

print("Connecting to MongoDB...")

mongodb_info = {
    "database": "sb_qchem",
    "collection": "new_tasks",
    "admin_user": "smblau_lbl.gov_readWrite",
    "admin_password": "baffler-underranger-sanguinely-distent-flukeworm",
    "host": "mongodb03.nersc.gov",
    "port": 27017,
    "aliases": {},
    "authSource": "sb_qchem"
}

client = pymongo.MongoClient(
    host=mongodb_info["host"],
    username=mongodb_info["admin_user"],
    password=mongodb_info["admin_password"],
    authSource=mongodb_info["authSource"]
)

database = client.sb_qchem
collection = database.euvl_mar_summary

# load current names and mpculeids

print("Done!")

print("Loading mol entries...")

os.chdir("C:/Users/jacob/Kinetiscope_Utils/classify_HiPRGen_reactions")
current_mpculeids_and_names = loadfn("name_full_mpculeid_092124.json")

# loads HiPRGen mol_entry objs

with open('mol_entries.pickle', 'rb') as f:

    try:

        mol_entries = pickle.load(f)

    except AttributeError:

        print("Need to use the pickle file from phase 2, not phase 1")
        print("The 'electron species' in phase 1 causes issues here")
        sys.exit()

mpculeid_mapping = {}
new_mpculeid_set = set()

for mol_entry in mol_entries:

    formula = mol_entry.formula.replace(" ", "")

    mpculeid = WriteNewMpculeid(
        mol_entry.graph, formula, mol_entry.charge, mol_entry.spin_multiplicity
    )

    mpculeid_mapping[mol_entry.entry_id] = mpculeid
    mol_entry.entry_id = mpculeid
    new_mpculeid_set.add(mpculeid)

print('Done!')

new_to_old = {value: key for key, value in mpculeid_mapping.items()}

print("Finding radicals that formed in kinetiscope...")

kinetiscope_files_dir = (
    r"G:\My Drive\Kinetiscope\production_simulations_092124"
)

correct_path_change_dir(kinetiscope_files_dir)

select_freq_file = "excitation_selection_freq_092524.txt"
reaction_name_file = "excitation_reactions_092524.txt"

relevant_reactions = get_nonzero_freq_reactions(
    select_freq_file,
    reaction_name_file,
    start_index=7,
    end_index=5384
)

chemical_dict = loadfn("name_full_mpculeid_092124.json")
corrected_mpculeids = {
    name: mpculeid_mapping[mpculeid]
    for name, mpculeid in chemical_dict.items()
    }

relevant_radical_dict = find_times_radical_formed(relevant_reactions,
                                                  corrected_mpculeids)

relevant_radical_set = find_relevant_radicals(relevant_radical_dict,
                                              threshold=0.01)

print("Done!")

underbonded_sites = find_underbonded_sites(
    relevant_radical_set, new_to_old, mol_entries
)

radical_mol_entries = find_mol_entries(relevant_radical_set, mol_entries)

all_radicals, radicals_to_recombine = build_radicals_sets(
    mol_entries, new_to_old, relevant_radical_set
)

recombinant_graphs = []
new_mol_entries = []
recombinant_mpculeids = set()
recombination_reactions = set()

for problem_radical in radicals_to_recombine:
    for radical in all_radicals:

        combined_graph, node_list, mapping = build_combined_graph(
            problem_radical, radical, radical_mol_entries, underbonded_sites
        )

        formula = write_new_formula(node_list)
        problem_underbonded_atoms = underbonded_sites[problem_radical]
        radical_underbonded_atoms = underbonded_sites[radical]

        for problem_underbonded_atom in problem_underbonded_atoms:
            for original_atom in radical_underbonded_atoms:

                atom_in_combined_graph = mapping[original_atom]

                recombinant_graph = generate_recombinant_graph(
                    combined_graph,
                    problem_underbonded_atom,
                    atom_in_combined_graph
                )

                # we're combining neutral radicals with neutral radicals, so
                # the resulting species will always be neutral and closed-shell

                mpculeid = WriteNewMpculeid(recombinant_graph, formula, 0, 1)

                if mpculeid not in new_mpculeid_set:

                    recombinant_mpculeids.add(mpculeid)
                    recombinant_graphs.append(recombinant_graph)
                    reactant_list = sorted([radical, problem_radical])
                    reaction = write_recombination_reaction(
                        reactant_list, mpculeid
                    )
                    mol = Molecule()
                    recombination_reactions.add(reaction)
                    # sys.exit()
print(f"Estimated number of recombinants: {len(radicals_to_recombine) * len(all_radicals)}")
# estimated number: 3721, w/o removing duplicates
print(f"Total number of new recombinants: {len(recombinant_mpculeids)}")
print(f"Total number of new reactions: {len(recombination_reactions)}")
        # combined_list = list(combined_graph.nodes(data=True))
        # for index, entry in enumerate(problem_list):
        #     corresponding_entry = combined_list[index]
        #     assert entry == corresponding_entry
        # sys.exit()
        # for site in problem_underbonded_atoms:
        #     for atom in radical_underbonded_atoms:
                
        
# for each recomb radical,
# generate a recomb rxn,
#  save the recomb mpculid and
# name
# if nbo_spins[atom_ind] > 0.09 and str(mol_graph.molecule[int(atom_ind)].specie) != "H": #how did we choose 0.09 as the cutoff?