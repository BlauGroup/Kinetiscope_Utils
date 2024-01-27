from monty.serialization import loadfn, dumpfn
from monty.json import MSONable
import copy
import pickle
import time
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash
import sqlite3

class Reaction(MSONable):
    def __init__(self, name, rxn_dict):
        self.name = name
        self.rxn_dict = rxn_dict
        self.reactant_charge = None
        self.product_charge = None
        self.reactant_spin = None
        self.product_spin = None
        
        self.set_charge_spin_values()
        
    def set_charge_spin_values(self):
        charge_spin_info = Reaction.get_charge_spin_info(self.name)
        self.reactant_charge = charge_spin_info['reactants_charge']
        self.product_charge = charge_spin_info['products_charge']
        self.reactant_spin = charge_spin_info['reactants_spin']
        self.product_spin = charge_spin_info['products_spin']

    def as_dict(self):
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "name": self.name,
            "rxn_dict": self.rxn_dict,
            "reactant_charge": self.reactant_charge,
            "product_charge": self.product_charge,
            "reactant_spin": self.reactant_spin,
            "product_spin": self.product_spin
        }

    @classmethod
    def from_dict(cls, d):
        instance = cls(name=d["name"], rxn_dict=d["rxn_dict"])
        instance.set_charge_spin_values()
        return instance
    
    @staticmethod
    def get_charge_spin_info(reaction_str):
        def parse_molecule(molecule_str):
            graph_hash, charge, spin = molecule_str.split('-')
            if 'm' in charge:
                charge = -int(charge[1:])
            else:
                charge = int(charge)
            spin = int(spin)
            return charge, spin
    
        reactants, products = reaction_str.split(' -> ')
        reactants = [parse_molecule(mol) for mol in reactants.split(' + ')]
        products = [parse_molecule(mol) for mol in products.split(' + ')]
    
        reactants_charge = sum(charge for charge, _ in reactants)
        products_charge = sum(charge for charge, _ in products)
    
        reactants_spin = 1 if all(spin == 1 for _, spin in reactants) else 2
        products_spin = 1 if all(spin == 1 for _, spin in products) else 2
    
        return {
            'name': reaction_str,
            'reactants_charge': reactants_charge,
            'products_charge': products_charge,
            'reactants_spin': reactants_spin,
            'products_spin': products_spin
        }

def write_reaction(reaction_dict): #writes reactions by mpculeids; currently order matters when comparing a reaction to itself
    """
    Takes a dictionary of associating reactants and products with their
    mpculeids and returns a reaction written as a string

    Parameters
    ----------
    reaction_dict : dictionary
        dictionary whose keys are "reactants" and "products" and whose values
        are a list of mpculeids associated with a given reaction

    Returns
    -------
    rxn_eqn: string
        the desired reaction written in the form: mpculeid1 + mpculeid 2 -> 
        mpculeid 3 + mpculeid 4
    """
    rxn_eqn = ""
    
    if len(reaction_dict["reactants"]) == 1:
        
        name = reaction_dict["reactants"][0]
        rxn_eqn = name + " "
        
    else:
        
        first_name = reaction_dict["reactants"][0]
        rxn_eqn = first_name + " + "
        second_name = reaction_dict["reactants"][1]
        rxn_eqn = rxn_eqn + second_name + " "
        # else:
        #     rxn_eqn = rxn_eqn + " "
            
    rxn_eqn = rxn_eqn + "-> "
    
    if len(reaction_dict["products"]) == 1:
        
        name = reaction_dict["products"][0]
        rxn_eqn = rxn_eqn + name
        
    else:
        
          first_name = reaction_dict["products"][0]
          second_name = reaction_dict["products"][1]
          rxn_eqn = rxn_eqn + first_name + " + " + second_name
          # else:
          #     rxn_eqn = rxn_eqn + first_name
          
    return rxn_eqn         
                                                                            
def gen_hashes_standardize_mpculeids(entry_graph, mpculeid, mpcule_dict, hash_dict, hash_set):
    '''
    Between versions, graph hashes from NetworkX can change. Mpculeids are generated
    from those hashes, and thus can also change. This function standardizes both
    the mpculeids and graph hashes associated with species in the network.
    
    Parameters
    ----------
    entry_graph : undirected NetworkX Graph
        graph associated with a species in the HiPRGen network
    mpculeid : string
        a string identifying a particular species in the network
    mpcule_dict : dict
        a dictionary associating old mpculeids with new ones
    hash_dict : dict
       a dictionary associating an mpculeid with its graph hash
    hash_set: set
        a set containing the hashes of all species in the network

    Returns
    -------
    None, modifies mpcule_dict, hash_dict, and hash_set
    '''
    entry_hash = weisfeiler_lehman_graph_hash(entry_graph, node_attr="specie")
    list_for_new_mpculeid = mpculeid.split('-') #first entry of this list is the old mpculeid
    list_for_new_mpculeid[0] = entry_hash 
    new_mpculeid = ''.join([list_for_new_mpculeid[0], list_for_new_mpculeid[1], '-', list_for_new_mpculeid[2], '-', list_for_new_mpculeid[3]]) #an mpculeid is a string of the form: hashcomposition-charge-spin so we need to rewrite the mpculeid that way
    mpcule_dict[mpculeid] = new_mpculeid
    hash_set.add(entry_hash) #TODO: whereever this is called, just test membership of reverse_hash_dict instead
    hash_dict[new_mpculeid] = entry_hash

def generate_sorted_reaction_tuple(rxn_copy, hash_dict):
    """
    Generates a 2-tuple of sorted tuples containing hash values and their frequencies associated with the reaction.

    Parameters
    ----------
    rxn_copy : dict
        Dictionary whose keys are "reactants" and "products" and whose values
        are lists of molecule ids associated with a given reaction.

    hash_dict : dict
        Dictionary that maps molecule ids to their corresponding hash values.

    Returns
    -------
    sorted_reaction_tuple : tuple
        A 2-tuple containing two sorted tuples of 2-tuples: the first for reactants and the second for products.
    """
    reactant_hashes = {}
    product_hashes = {}

    for identity, side in rxn_copy.items():
        for member in side:
            if hash_dict.get(member, False):
                h = hash_dict[member]
                current_hashes = reactant_hashes if identity == "reactants" else product_hashes
                current_hashes[h] = current_hashes.get(h, 0) + 1

    reactant_tuples = tuple((hash_value, freq) for hash_value, freq in reactant_hashes.items())
    product_tuples = tuple((hash_value, freq) for hash_value, freq in product_hashes.items())

    # Sort reactant and product tuples separately
    sorted_reactant_tuples = tuple(sorted(reactant_tuples, key=lambda x: (x[0], x[1])))
    sorted_product_tuples = tuple(sorted(product_tuples, key=lambda x: (x[0], x[1])))

    return sorted_reactant_tuples, sorted_product_tuples

# def update_rxn_copy(rxn_dict, mpcule_dict):
#     """
#     Updates a reaction dictionary copy with the most recent molecule ids from an mpcule dictionary.

#     Parameters
#     ----------
#     rxn_dict : dict
#         Original reaction dictionary whose keys are "reactants" and "products" and whose values
#         are lists of molecule ids associated with a given reaction.

#     mpcule_dict : dict
#         Dictionary that maps old molecule ids to their most recent ones.

#     Returns
#     -------
#     rxn_copy : dict
#         Updated reaction dictionary copy with the most recent molecule ids.
#     """
#     rxn_copy = copy.deepcopy(rxn_dict)

#     for identity, side in rxn_dict.items():
#         for index, member in enumerate(side):
#             if mpcule_dict.get(member, False):
#                 rxn_copy[identity][index] = mpcule_dict.get(member)
#             else:
#                 if not member: #catches the "electron" species
#                     rxn_copy[identity].pop(index)
#                 else:
#                     raise KeyError('mpculeid not found in the mpcule_dict')

#     return rxn_copy

def update_rxn_copy(rxn_dict, mpcule_dict, reverse_mpcule_dict):
    """
    Updates a reaction dictionary copy with the most recent molecule ids from an mpcule dictionary.

    Parameters
    ----------
    rxn_dict : dict
        Original reaction dictionary whose keys are "reactants" and "products" and whose values
        are lists of molecule ids associated with a given reaction.

    mpcule_dict : dict
        Dictionary that maps old molecule ids to their most recent ones.

    Returns
    -------
    rxn_copy : dict
        Updated reaction dictionary copy with the most recent molecule ids.
    """
    # Assuming mpcule_dict and reverse_mpcule_dict are dictionaries
    mpcule_dict_first_key_type = type(next(iter(mpcule_dict), None))
    reverse_mpcule_dict_first_key_type = type(next(iter(reverse_mpcule_dict), None))
    
    print(f"Type of first key in mpcule_dict: {mpcule_dict_first_key_type}")
    print(f"Type of first key in reverse_mpcule_dict: {reverse_mpcule_dict_first_key_type}")
    
    rxn_copy = copy.deepcopy(rxn_dict)
    
    for identity, side in rxn_dict.items():
        for index, member in enumerate(side):
            print(f"Processing: identity={identity}, side={side}, index={index}, member={member}")
    
            # Add the following print statements
            print(f"Type of first key in mpcule_dict: {mpcule_dict_first_key_type}")
            print(f"Type of first key in reverse_mpcule_dict: {reverse_mpcule_dict_first_key_type}")
            print(f"Type of member: {type(member)}")
    
            if mpcule_dict.get(member, False):
                updated_member = mpcule_dict.get(member)
                rxn_copy[identity][index] = updated_member
                print(f"Updated {member} to {updated_member} at rxn_copy[{identity}][{index}]")
            elif reverse_mpcule_dict.get(member, False):
                rxn_copy[identity][index] = member
                print(f"Kept {member} unchanged at rxn_copy[{identity}][{index}]")
            else:
                if not member:  # catches the "electron" species
                    rxn_copy[identity].pop(index)
                    print(f"Removed {member} at rxn_copy[{identity}][{index}]")
                else:
                    raise KeyError(f'mpculeid {member} not found in mpcule dicts')

    return rxn_copy

def process_sqlite_reaction_tuple(row, ind_dict):
    reactants = [ind_dict.get(row[0], False)]
    reactants += [ind_dict.get(row[1], False)] if row[1] != -1 or ind_dict.get(row[1], 0) else []
    products = [ind_dict.get(row[2], False)]
    products += [ind_dict.get(row[3], False)] if row[3] != -1 or ind_dict.get(row[3], 0) else []

    return {"reactants": reactants, "products": products}
 
if __name__ == "__main__":
       
    start = time.time()
    
    print("Generating dictionaries...")
    
    mol_pickle = "mol_entries.pickle"
    
    with open(mol_pickle, "rb") as f:
        mol_entries = pickle.load(f)
    
    hash_set = set()
    hash_dict = {}
    mpcule_dict = {}
    ind_dict = {}
    
    for entry in mol_entries:
        
        graph = entry.graph
        mpculeid = entry.entry_id
        if not mpculeid: #skips the "electron" species, which has no id
            continue
        ind = entry.ind
        ind_dict[ind] = mpculeid
        gen_hashes_standardize_mpculeids(graph, mpculeid, mpcule_dict, hash_dict, hash_set)
        
    reverse_hash_dict = {value: key for key, value in hash_dict.items()}
    reverse_mpcule_dict = {value: key for key, value in mpcule_dict.items()}
    
    requested_json = "/global/home/groups/lr_mp/smblau/jrmilton/HiPRGen/euvl_TSreactions_041823.json" #currently contains 134 duplicates
    requested_list = loadfn(requested_json)
    
    p1_json = "/global/home/groups/lr_mp/smblau/jrmilton/HiPRGen/reaction_tally_p1.json"
    p1_dict = loadfn(p1_json)
    p2_json = "/global/home/groups/lr_mp/smblau/jrmilton/HiPRGen/reaction_tally_p2.json"
    p2_dict = loadfn(p2_json)
    combined_dict = {**p1_dict, **p2_dict}
    
    unique_requested_list = []
    
    for item in requested_list:
        if item not in unique_requested_list:
            unique_requested_list.append(item)
    
    requested_rxn_dict = {}
    
    for rxn_dict in unique_requested_list:
        
        rxn_copy = update_rxn_copy(rxn_dict, mpcule_dict, reverse_mpcule_dict)
        rxn_tuple = generate_sorted_reaction_tuple(rxn_copy, hash_dict)
        
        # Check if the reaction tuple is already present
        
        if rxn_tuple in requested_rxn_dict:
            
            existing_entry = requested_rxn_dict[rxn_tuple]
            new_reaction_result = write_reaction(rxn_copy)
            existing_entry.append(Reaction(new_reaction_result, rxn_copy))
                
        else:
            
            requested_rxn_dict[rxn_tuple] = [Reaction(write_reaction(rxn_copy), rxn_copy)]
    
    tally_rxn_dict = {}
    
    for rxn_dict in combined_dict['reactions'].values():
        
        rxn_copy = update_rxn_copy(rxn_dict, mpcule_dict, reverse_mpcule_dict)
        rxn_tuple = generate_sorted_reaction_tuple(rxn_copy, hash_dict)
        
        # Check if the reaction tuple is already present
        
        unrequested_reaction = rxn_tuple not in requested_rxn_dict
        tuple_present_not_reaction = rxn_tuple in requested_rxn_dict and Reaction(write_reaction(rxn_copy), rxn_copy) not in requested_rxn_dict[rxn_tuple]
        
        if unrequested_reaction or tuple_present_not_reaction:
            # Want this dict to contain only unrequested reactions from the tally files
        
            if rxn_tuple in tally_rxn_dict:
                
                existing_entry = tally_rxn_dict[rxn_tuple]
                new_reaction_result = write_reaction(rxn_copy)
                existing_entry.append(Reaction(new_reaction_result, rxn_copy))
                
            else:
                
                tally_rxn_dict[rxn_tuple] = [Reaction(write_reaction(rxn_copy), rxn_copy)]
        
    total_entries = sum(len(entry) for entry in requested_rxn_dict.values())
    assert total_entries == len(unique_requested_list), f"Mismatch between total reactions ({total_entries}) and member_dict length ({len(unique_requested_list)})"
    
    database_p1 = "/global/home/groups/lr_mp/smblau/HiPRGen/euvl_feb_phase1/phase_1/rn.sqlite"
    database_p2 = "/global/home/groups/lr_mp/smblau/HiPRGen/euvl_feb_phase2/phase_2/rn.sqlite"
    
    full_rxn_dict = {}
    
    p1_con = sqlite3.connect(database_p1)
    cur = p1_con.cursor()
    
    sql_get_reactions = """
        SELECT reactant_1, reactant_2, product_1, product_2 FROM reactions;
    """
    
    p1_full_rxns = cur.execute(sql_get_reactions)
    
    for row in cur.execute(sql_get_reactions):
        
        rxn_dict = process_sqlite_reaction_tuple(row, ind_dict)        
        rxn_copy = update_rxn_copy(rxn_dict, mpcule_dict, reverse_mpcule_dict)
        rxn_tuple = generate_sorted_reaction_tuple(rxn_copy, hash_dict)
            
        # Check if the reaction tuple is already present
        
        rxn_tuple_not_in_either = rxn_tuple not in requested_rxn_dict or rxn_tuple not in tally_rxn_dict
        tuple_present_not_reaction_requested = rxn_tuple in requested_rxn_dict and Reaction(write_reaction(rxn_copy), rxn_copy) not in requested_rxn_dict[rxn_tuple]
        tuple_present_not_reaction_tallies = rxn_tuple in tally_rxn_dict and Reaction(write_reaction(rxn_copy), rxn_copy) not in tally_rxn_dict[rxn_tuple]
        
        if rxn_tuple_not_in_either or (tuple_present_not_reaction_requested or tuple_present_not_reaction_tallies):
            # Want this dict to contain only unrequested reactions from the tally files
        
            if rxn_tuple in full_rxn_dict:
                
                existing_entry = full_rxn_dict[rxn_tuple]
                new_reaction_result = write_reaction(rxn_copy)
                existing_entry.append(Reaction(new_reaction_result, rxn_copy))
            else:
                
                full_rxn_dict[rxn_tuple] = [Reaction(write_reaction(rxn_copy), rxn_copy)]
        
    p2_con = sqlite3.connect(database_p2)
    cur = p2_con.cursor()
       
    for row in cur.execute(sql_get_reactions):
        
        rxn_dict = process_sqlite_reaction_tuple(row, ind_dict)        
        rxn_copy = update_rxn_copy(rxn_dict, mpcule_dict, reverse_mpcule_dict)
        rxn_tuple = generate_sorted_reaction_tuple(rxn_copy, hash_dict)
            
        # Check if the reaction tuple is already present
        
        rxn_tuple_not_in_either = rxn_tuple not in requested_rxn_dict or rxn_tuple not in tally_rxn_dict
        tuple_present_not_reaction_requested = rxn_tuple in requested_rxn_dict and Reaction(write_reaction(rxn_copy), rxn_copy) not in requested_rxn_dict[rxn_tuple]
        tuple_present_not_reaction_tallies = rxn_tuple in tally_rxn_dict and Reaction(write_reaction(rxn_copy), rxn_copy) not in tally_rxn_dict[rxn_tuple]
        
        if rxn_tuple_not_in_either or (tuple_present_not_reaction_requested or tuple_present_not_reaction_tallies):
            # Want this dict to contain only unrequested reactions from the tally files
        
            if rxn_tuple in full_rxn_dict:
                
                existing_entry = full_rxn_dict[rxn_tuple]
                new_reaction_result = write_reaction(rxn_copy)
                existing_entry.append(Reaction(new_reaction_result, rxn_copy))
                
            else:
                full_rxn_dict[rxn_tuple] = [Reaction(write_reaction(rxn_copy), rxn_copy)]
    
    requested_json = "/global/home/groups/lr_mp/smblau/jrmilton/HiPRGen/requested_rxns.json"
    tally_json = "/global/home/groups/lr_mp/smblau/jrmilton/HiPRGen/tally_rxns.json"
    total_json = "/global/home/groups/lr_mp/smblau/jrmilton/HiPRGen/total_rxns.json"
    
    json_files = [requested_json, tally_json, total_json]
    
    for json_file in json_files:
        if json_file == requested_json:
            current_dict = requested_rxn_dict
        elif json_file == tally_json:
            current_dict = tally_rxn_dict
        else:
            current_dict = full_rxn_dict
    
        str_keys_dict = {str(key): value for key, value in current_dict.items()}
        dumpfn(str_keys_dict, str(json_file))
    
    end = time.time()
    time_min = (end - start)/60
    print(round(time_min, 2))