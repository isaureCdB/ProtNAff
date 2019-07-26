"""Adds PFAM and UniProt information to structures.json
Requires the existence of /tmp/pdb_chain_pfam.tsv.gz
(downloadable from ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv )
"""

import json, os
import pandas as pd

filename = "/tmp/pdb_chain_pfam.tsv.gz"
if not os.path.exists(filename):
    msg = "Requires the existence of /tmp/pdb_chain_pfam.tsv.gz, downloadable from ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv"
    raise Exception(msg)
pdb_chain_pfam = pd.read_csv(filename, dtype=str, skiprows=1, delimiter="\t")

mapping = {}                                                                                    

for pdb, chain, uniprot, pfam in zip(
  pdb_chain_pfam["PDB"], pdb_chain_pfam["CHAIN"], 
  pdb_chain_pfam["SP_PRIMARY"], pdb_chain_pfam["PFAM_ID"]
): 
    mapping[pdb.upper(), chain] = uniprot, pfam 

with open("structures.json") as f:
    structures = json.load(f)

for pdb, struc in structures.items():
    pfam = {}
    uniprot = {}
    for chain in struc["protchains"]:
        key = mapping.get((pdb, chain))
        if key is None:
            continue
        uniprot[chain] = key[0]
        pfam[chain] = key[1]
        
    struc["pfam"] = pfam
    struc["uniprot"] = uniprot
    
with open("structures.json", "w") as f:
    json.dump(structures, f, indent=2, sort_keys=True)
