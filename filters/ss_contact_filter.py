#!/usr/bin/env python3

import json
import sys

"""
Select ssRNA in interaction with protein.
"""

def ss_contact_filter(pdb_info, chain_id):
    """
    Input: pdb_info: information of the pdb in json format
           chain_id: the name of the RNA chain
    Output: nuclfrag, a set of nucleotide of interest
    """

    ss_set = set(["L", "T", "S", "J", "B", "I"])

    chain_id = "chain_" + chain_id
    # print(pdb_info['ss'][chain_id])
    ss = pdb_info['ss'][chain_id]

    interf_prot = pdb_info["interface_protein"]["model_1"][chain_id]
    # print(interf_prot)
    for element in interf_prot:
        element = element.split("_")[1]

    # List of interface nucleotide resid
    nucl_interf = list(interf_prot.keys())
    nucl_interf = [element.split("_")[1] for element in nucl_interf]
    # List of single-stranded nucleotides resid
    nucl_ss = [n.split("_")[1] for n in ss if ss[n][0] in ss_set]

    # List of single-stranded interface nucleotide resid
    nucl_interf_ss = [int(n) for n in nucl_interf if n in nucl_ss]
    nucl_interf_ss.sort()

    return nucl_interf_ss


def main():
    js = json.load(open(sys.argv[1])) #structures.json

    for pdb_id in js:
        pdb_info = js[pdb_id]
        for chain_id in pdb_info['nachains']:
            nuclfrag = ss_contact_filter(pdb_info, chain_id)
            print("{} {} {}".format(pdb_id, chain_id, nuclfrag))

if __name__ == '__main__':
    main()
