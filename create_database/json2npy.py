#!/usr/bin/env python3
import numpy as np, sys

dtype = [
    ("motif", "S3"),
    ("frag", np.uint32),
    ("structure", "S4"),
    ("model", np.uint16),
    ("chain", "S1"), #nachains?
    ("resid", "S5", 3),
    ("indices", np.uint32, 3), #residue index
    ("missing_atoms", np.uint8, 3),
    ("seq", "S3"),
    ("clust0.2Aaa", np.uint16),
    ("clust0.2Aaa_center", bool),
    ("clust1Aaa", np.uint16),
    ("clust1Aaa_center", bool),
    ("clust2Aaa", np.uint16),
    ("clust2Aaa_center", bool)
]

dtype = np.dtype(dtype)
import json
jfrags=json.load(open(sys.argv[1])) #fragments_clust-aa_missing.json
size = 0
for motif in jfrags:
    size += len(jfrags[motif])
print(size)
frags = np.zeros(size,dtype)

count = 0
for motif in jfrags:
    for jfragnr, jfrag in jfrags[motif].items():
        frag = frags[count]
        frag["motif"] = motif
        frag["frag"] = int(jfragnr)
        frag["chain"] = jfrag["chain"]
        frag["clust0.2Aaa"] = jfrag["clust0.2Aaa"]
        frag["clust0.2Aaa_center"] = jfrag["clust0.2Aaa\\_center"]
        frag["clust1Aaa"] = jfrag["clust1Aaa"]
        frag["clust1Aaa_center"] = jfrag["clust1Aaa\\_center"]
        frag["clust2Aaa"] = jfrag["clust2Aaa"]
        frag["clust2Aaa_center"] = jfrag["clust2Aaa\\_center"]
        frag["indices"] = jfrag["indices"]
        frag["missing_atoms"] = jfrag["missing_atoms"]
        frag["model"] = jfrag["model"]
        frag["resid"] = jfrag["resid"]
        frag["seq"] = jfrag["seq"]
        frag["structure"] = jfrag["structure"]
        count += 1
print(count)
outp = sys.argv[1].split(".json")[0] + ".npy"
np.save(outp, frags)
