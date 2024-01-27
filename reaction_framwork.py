# -*- coding: utf-8 -*-
from monty.serialization import loadfn, dumpfn
from monty.json import MSONable
import bson


"""
Created on Sat Nov 11 18:35:47 2023

@author: jake
"""

class ChemicalProcess(MSONable):
    def __init__(self, name, reaction_dict):
        self.name = name
        self.reaction_dict = reaction_dict
        self.reactant_charge = None
        self.product_charge = None
        self.reactant_spin = None
        self.product_spin = None

    def set_charges_and_spins(self):
        # Placeholder method, you should implement this based on your needs
        pass

    def as_dict(self):
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "name": self.name,
            "reaction_dict": self.reaction_dict,
            "reactant_charge": self.reactant_charge,
            "product_charge": self.product_charge,
            "reactant_spin": self.reactant_spin,
            "product_spin": self.product_spin
        }

    @classmethod
    def from_dict(cls, d):
        instance = cls(name=d["name"], reaction_dict=d["reaction_dict"])
        instance.reactant_charge = d["reactant_charge"]
        instance.product_charge = d["product_charge"]
        instance.reactant_spin = d["reactant_spin"]
        instance.product_spin = d["product_spin"]
        return instance
        
class TransitionState(ChemicalProcess):
    def __init__(self, name, rxn_set, forward_tuple, reverse_tuple,
                  TS_molecule, forward_complex, reverse_complex,
                  TS_self_test, TS_other_test):
        self.forward_tuple = forward_tuple
        self.reverse_tuple = reverse_tuple
        self.TS_molecule = TS_molecule
        self.forward_complex = forward_complex
        self.reverse_complex = reverse_complex
        self.TS_self_test = TS_self_test
        self.TS_other_test = TS_other_test

    def __str__(self):
        return f"TransitionState(name={self.name}, rxn_dict={self.rxn_dict}, " \
                f"rxn_set={self.rxn_set}, forward_tuple={self.forward_tuple}, " \
                f"reverse_tuple={self.reverse_tuple}, TS_molecules={self.TS_molecules}, " \
                f"forward_molecules={self.forward_molecules}, reverse_molecules={self.reverse_molecules}, " \
                f"TS_self_test={self.TS_self_test}, TS_other_test={self.TS_other_test})"
                
def process_ts_entry(entry, job_record):
    """
    Process a transition state (TS) entry from Jaguar data and updates the job record.
    
    Parameters
    ----------
    entry : dict
        The entry containing TS information from Jaguar data.
    job_record : dict
        The dictionary to update with TS information.
    
    Returns
    -------
    None, but updates the dictionary with this information.
    """
    
    name = entry["name"]

    # Determine if the optimization is completed based on where or not there is an output geometry present
    geo_opt_completed = entry["output"].get("geopt", False) or entry["input"].get("geopt", False)

    if geo_opt_completed:
        
        source = "output"
        
    else:
        
        source = "input"

    final_geo_opt = entry[source]["molecule"]

    if name in job_record:
        
        job_record[name].append({"ts_molecule": final_geo_opt, "geo_from": source})
    else:
        
        job_record[name] = [{"ts_molecule": final_geo_opt, "geo_from": source}]

def process_opt_entry(entry, job_record):
    """
    Process an optimization (OPT) entry from Jaguar data and update the job record.
    
    Parameters
    ----------
    entry : dict
        The entry containing optimization information from Jaguar data.
    job_record : dict
        The dictionary to update with optimization information.
    
    Returns
    -------
    None, but updates the dictionary with this information.
    """
    
    split_name = entry["name"].split()
    base_name = split_name[0] + " " + split_name[1]
    direction = "forward" if split_name[2] == "forwards" else split_name[2]

    direction_mol_name = direction + "_mol"

    duplicate_geom = any(geo_dict.get(direction_mol_name, False) for geo_dict in job_record.get(base_name, []))

    if entry["output"].get("geopt", False):
        
        direction_mol = entry["output"]["molecule"]
        source = "output"
        
    else:
        direction_mol = entry["input"]["molecule"]
        source = "input"

    if base_name in job_record:
        
        job_record[base_name].append({direction_mol_name: direction_mol, "geo_from": source, "duplicate": duplicate_geom, 'mpculeids': {}})
        
    else:
        job_record[base_name] = [{direction_mol_name: direction_mol, "geo_from": source, "duplicate": duplicate_geom, 'mpculeids': {}}]

def load_id_name_dict(autoTS_queue_bson):
    """
    Load an ID-name dictionary from an AutoTS queue BSON file.

    Parameters
    ----------
    autoTS_queue_bson : str
        The path to the AutoTS queue BSON file.

    Returns
    -------
    dict
        The ID-name dictionary initialized with information from the BSON file.
    """
    
    id_name_dict = {}

    with open(autoTS_queue_bson, "rb") as f:
        id_name_data = bson.decode_file_iter(f)
        for entry in id_name_data:
                id_name_dict[str(entry["rxnid"])] = {"mpcule_name": entry["name"], "rxn_id": "", "TS_molecules": [],
                                                      "forward_molecules": [], "reverse_molecules": [], "TS_self_test": False,
                                                      "TS_other_test": (False, 0)}

    return id_name_dict

def update_id_name_dict(id_name_dict, rxn_id, mol, basis):
    """
    Update an ID-name dictionary with information from an AutoTS output dictionary.

    Parameters
    ----------
    id_name_dict : dict
        The ID-name dictionary to update.
    rxn_id : str
        The reaction ID.
    mol : dict
        The molecule information dictionary.
    basis : str
        The basis set information.

    Returns
    -------
    None
        Updates the ID-name dictionary with the provided information.
    """
    
    jaguar_id_num = rxn_id.split("_")[0]
    current_molecules_key = ""

    if mol.get("ts_molecule", False):
        current_molecules_key = "TS_molecules"
    elif mol.get("reverse_mol", False):
        current_molecules_key = "reverse_molecules"
    elif mol.get("forward_mol", False):
        current_molecules_key = "forward_molecules"
    else:
        raise ValueError("Molecule is somehow not a TS, forward, or reverse!")

    current_geometries = id_name_dict[jaguar_id_num][current_molecules_key]
    updated_mol = mol.copy()
    updated_mol["basis"] = basis
    current_geometries.append(updated_mol)
    id_name_dict[jaguar_id_num][current_molecules_key] = current_geometries
        
input_bson = "G:/.shortcut-targets-by-id/1e9g_4Hj-KL1HwtAf4qE0vww634guxT9n/euvl/jaguar_data.bson"
output_json = "Found_TSs_and_minima.json"

with open(input_bson, "rb") as f:
    jaguar_data = bson.decode_file_iter(f)
    job_record = {}

    for entry in jaguar_data:
        if entry["job_type"] == "ts":
            process_ts_entry(entry, job_record)
        elif entry["job_type"] == "opt":
            process_opt_entry(entry, job_record)

# dumpfn(job_record, output_json)

autoTS_queue_bson = "G:/.shortcut-targets-by-id/1e9g_4Hj-KL1HwtAf4qE0vww634guxT9n/euvl/autots_queue.bson" 
id_name_dict = {}

with open(autoTS_queue_bson,"rb") as f: #generates a dictionary associating a jaguar id with molecules found for that id
    id_name_data = bson.decode_file_iter(f)
    for entry in id_name_data:
        id_name_dict[str(entry["rxnid"])] = {"mpcule_name": entry["name"], "rxn_id": "", "TS_molecules":[], "forward_molecules":[], "reverse_molecules": [], "TS_self_test": False, "TS_other_test": (False, 0)}

autoTS_queue_bson = "G:/.shortcut-targets-by-id/1e9g_4Hj-KL1HwtAf4qE0vww634guxT9n/euvl/autots_queue.bson"
id_name_dict = load_id_name_dict(autoTS_queue_bson)

TS_json = "G:/My Drive/CRNs/AutoTS_output/Found_TSs_and_minima.json"
TS_dict = loadfn(TS_json)
rxn_dict_path = "G:/My Drive/CRNs/AutoTS_output/rxn_dict.json"

for rxn_id, mol_list in TS_dict.items():
    for mol in mol_list:
        underscore_split = rxn_id.split("_")
        space_split = underscore_split[1].split()  # don't need the space
        basis = space_split[1]
        jaguar_id_num = underscore_split[0]

        update_id_name_dict(id_name_dict, rxn_id, mol, basis)

# dumpfn(id_name_dict, rxn_dict_path)