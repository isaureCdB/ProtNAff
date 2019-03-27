#!/usr/bin/env python3
import json, sys, os
import argparse
from collections import OrderedDict
'''
search sequences ssRNA and dsRNA and write a new json file.
'''

parser = argparse.ArgumentParser(description=__doc__,
                          formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("file", type=str, help="the json file on which you want to find ssRNA and dsRNA.")
args = parser.parse_args()

# Put input in variable names
json_file=args.file

js = json.load(open(json_file)) #x3dna.json

def query(pdb_id, chain_id):
    interf_prot = pdb_id["interface_protein"]["1"][chain_id]
    # List of interface nucleotide resid
    nucl_interf = [n for n in interf_prot if len(interf_prot[n]) > 0]
    if len(nucl_interf) > 0:
        return True
    else:
        return False

js_dict = OrderedDict()

for k in js:
    pdb_id = js[k]
    for chain_id in pdb_id['nachains']:
        keep = query(pdb_id, chain_id)
        if keep:
            js_dict[k] = pdb_id
            break

with open("{}_test.json".format(json_file.split(".")[0]), 'w') as outfile:
    json.dump(js_dict, outfile, indent=2)
