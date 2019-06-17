#!/usr/bin/env python3
import numpy as np, sys
dtype = [
    ("motif", "S3"),
    ("frag", np.uint32),
    ("structure", "S5"),
    ("model", np.uint16),
    ("chain", "S1"), #nachains?
    ("resid", "S5", 3),
    ("indices", np.uint32, 3), #residue index
    ("missing_atoms", np.uint8, 3),
    ("seq", "S3"),
    ("clust0.2", np.uint16),
    ("clust0.2_center", bool),
    ("clust1.0", np.uint16),
    ("clust1.0_center", bool),
    ("clust2.0", np.uint16),
    ("clust2.0_center", bool)
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
        frag["clust0.2"] = jfrag["clust0.2"]
        frag["clust0.2_center"] = jfrag["clust0.2_center"]
        frag["clust1.0"] = jfrag["clust1.0"]
        frag["clust1.0_center"] = jfrag["clust1.0_center"]
        frag["clust2.0"] = jfrag["clust2.0"]
        frag["clust2.0_center"] = jfrag["clust2.0_center"]
        frag["indices"] = jfrag["indices"]
        frag["missing_atoms"] = jfrag["missing_atoms"]
        frag["model"] = jfrag["model"]
        frag["resid"] = jfrag["resid"]
        frag["seq"] = jfrag["seq"]
        frag["structure"] = jfrag["structure"]
        count += 1
print(count)
np.save(sys.argv[2], frags) # fragments_clust-aa_missing.npy
