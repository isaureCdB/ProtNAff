import sys

infile = sys.argv[1]
outfile = sys.argv[2]

import numpy as np

dtype = [
    ("motif", "S3"),
    ("frag", np.uint32),
    ("structure", "S4"),
    ("model", np.uint16),
    ("chain", "S1"),  # nachains?
    ("resid", "S5", 3),
    ("indices", np.uint32, 3),  # residue index
    ("missing_atoms", np.uint8, 3),
    ("seq", "S3"),
    ("clust0.2A", np.uint16),
    ("clust0.2A_center", bool),
    ("clust0.5A", np.uint16),
    ("clust0.5A_center", bool),
    ("clust1A", np.uint16),
    ("clust1A_center", bool),
]
dtype = np.dtype(dtype)

import json

jfrags = json.load(open(infile))
size = 0
for motif in jfrags:
    size += len(jfrags[motif])
print(size)
frags = np.zeros(size, dtype)

count = 0
for motif in jfrags:
    for jfragnr, jfrag in jfrags[motif].items():
        frag = frags[count]
        frag["motif"] = motif
        frag["frag"] = int(jfragnr)
        for field in dtype.fields:
            if field in ("frag", "motif"):
                continue
            frag[field] = jfrag[field]
        count += 1
print(count)
np.save(outfile, frags)
