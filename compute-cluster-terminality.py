"""
Generates a terminality array for each fragment cluster.

The terminality array is 5 bools [A,B,C,D,E]
A: the fragment is present at a 5'-terminus (beginning)
B: the fragment is present next to 5'-terminus (beginning + 1)
C: the fragment is present in the middle (at a position not next to a terminus)
D: the fragment is present next to 3'-terminus (end - 1)
E: the fragment is present at a 3'-terminus (end)

Input: 

- fragments_clust.npy 
(fragments.json => assign_clusters.py => query/json2npy.py => fragments_clust.npy)
- structures.json (to get the sequence length.)
- name of cluster field in each fragment (e.g clust0.5A)

Generates files {X}-{Y}-terminality.npy, where:
    - X is the library sequence ([A/C][A/C][A/C], e.g. ACA)
    - Y is the cluster field name
"""

import sys, json, numpy as np

fragments_clust_npy, structures_json, cluster_field = sys.argv[1:4]

fragments = np.load(fragments_clust_npy)
structures = json.load(open(structures_json))

seqlen = {}
for pdb_code, struc in structures.items():
    for chain, seq in struc["sequence"].items():
        assert chain.startswith("chain_"), chain
        chain = chain[len("chain_") :]
        seqlen[pdb_code, chain] = len(seq)

for motif in ("AAA", "AAC", "ACA", "ACC", "CAA", "CAC", "CCA", "CCC"):
    frags = fragments[fragments["motif"] == motif.encode()]
    nclust = frags[cluster_field].max()
    terminality = np.zeros((nclust, 5), bool)
    for frag in frags:
        code = frag["structure"].decode()
        if code not in structures:
            continue
        seql = seqlen[code, frag["chain"].decode()]
        index = frag["indices"][0]
        clust = frag[cluster_field]  # count from one
        if clust == 0:  # for some reason, unclustered
            continue
        term = terminality[clust - 1]
        middle = True
        assert index > 0
        if index == 1:
            term[0] = True
            middle = False
        if index == 2:
            term[1] = True
            middle = False
        if index == seql - 2:
            term[4] = True
            middle = False
        if index == seql - 3:
            term[3] = True
            middle = False
        if middle:
            term[2] = True
    print(f"{motif}, {len(terminality)} clusters, sums: {terminality.sum(axis=0)}")
    filename = f"{motif}-{cluster_field}-terminality.npy"
    np.save(filename, terminality)
