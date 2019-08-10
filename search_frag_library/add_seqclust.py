"""Adds sequence cluster information
Requires the existence of /tmp/bc-X.out, where X is 30, 40, 50, 70, 90, 95, 100.
(downloadable from ftp://resources.rcsb.org/sequence/clusters/)
"""

import json, os

mapping = {}
seqids = 30, 40, 50, 70, 90, 95, 100
for seqid in seqids:
    filename = "/tmp/bc-%d.out" % seqid
    if not os.path.exists(filename):
        msg = "Requires the existence of %s, downloadable from ftp://resources.rcsb.org/sequence/clusters/"
        raise Exception(msg % filename)
for seqid in seqids:
    curr_mapping = {}
    mapping[seqid] = curr_mapping
    filename = "/tmp/bc-%d.out" % seqid
    with open(filename) as f:
        clustnr = 0
        for lnr, line in enumerate(f):
            fields = line.split()
            if not len(fields):
                continue
            clustnr += 1
            for pdb_chain in fields:
                pdb, chain = pdb_chain.split("_")
                curr_mapping[pdb.upper(), chain] = clustnr

with open("structures.json") as f:
    structures = json.load(f)

for pdb, struc in structures.items():
    seqclust = {}
    for chain in struc["protchains"]:
        curr_seqclust = {}
        seqclust[chain] = curr_seqclust
        for seqid in seqids:
            curr_mapping = mapping[seqid]
            clustnr = curr_mapping.get((pdb[:4], chain))
            if clustnr is None:
                continue
            curr_seqclust[seqid] = clustnr
    struc["seqclust"] = seqclust

with open("structures.json", "w") as f:
    json.dump(structures, f, indent=2, sort_keys=True)
