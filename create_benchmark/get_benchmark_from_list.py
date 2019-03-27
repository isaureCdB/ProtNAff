#!/usr/bin/env python3
import json, sys, os
import argparse
from queries.query_ss import query_ss
import shutil
from benchmark_scripts.extract_chains import extract_chains
from benchmark_scripts.create_pdb import create_pdb
'''
search sequences ssRNA and dsRNA and write a new json file.
'''

parser = argparse.ArgumentParser(description=__doc__,
                          formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("file", type=str, help="The json file on which you want to find ssRNA and dsRNA.")
parser.add_argument("list", type=str, help="The list on which you want to apply your query.")
parser.add_argument("output", type=str, help="The directory where you want your output.")
args = parser.parse_args()

# Put input in variable names
json_file=args.file

outpdir = args.output # ~/projets/ssRNA/benchmark/
listPDB = args.list

js = json.load(open(json_file)) #x3dna.json



# For each PDB file in the json
for pdb_id in open(listPDB, 'r'):

    try:

        pdb_info = js[pdb_id]

        dirPDB = "{}/{}".format(outpdir, pdb_id)

        dirList = [ dir for dir in os.listdir(dirPDB) if not os.path.isfile(os.path.join(dirPDB,dir)) ]

        for dir in dirList:
            


        with open("{}/query_successful.txt".format(outpdir), 'a') as f:
            f.write("{}\n".format(pdb_id))


#shutil.rmtree(tmpPDB)
