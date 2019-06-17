#!/usr/bin/env python3

import os, sys

# def in_list(L, element):
#     for couple in L:
#         if element in couple:
#             return True
#     return False
ATTRACT=os.environ.copy()["ATTRACTDIR"]
def extract_chains(pdb_id, outputdir, verify_hetatm=False):

    os.system('$SCRIPTS/pdb_download_biological_assembly {} {}'.format(pdb_id, outputdir))

    i = 1

    if verify_hetatm:
        list_hetatm = []
        with open(ATTRACT+"/../allatom/rna-mutate.list", 'r') as f:
            lines = f.readlines()
            for line in lines:
                list_hetatm.append(line.split()[0])

    interaction_chains = {}
    while True:
        try:
            with open('{}/{}.pdb{}'.format(outputdir, pdb_id, i), 'r') as file:
                interaction_chains[i] = {"protein":[], "rna":[]}
                lines = file.readlines()
                for line in lines:
                    if line.startswith('ATOM'):
                        if line[13:15] == "CA":
                            if line[21] not in interaction_chains[i]["protein"]:
                                interaction_chains[i]["protein"].append(line[21])
                        elif line[13:16] == "O2\'":
                            if line[21] not in interaction_chains[i]["rna"]:
                                interaction_chains[i]["rna"].append(line[21])
                    if line.startswith('HETATM'):
                        if line[13:16] == "O2\'" and line[18:20] in list_hetatm:
                            if line[21] not in interaction_chains[i]["rna"]:
                                interaction_chains[i]["rna"].append(line[21])
            i += 1
        except Exception as e:
            print(e, file=sys.stderr   )
            break
    print(interaction_chains)
    return interaction_chains
