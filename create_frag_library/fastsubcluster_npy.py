#!/usr/bin/env python3

import sys
import numpy as np
from scipy.spatial.distance import pdist, squareform
from math import sqrt
from npy import npy3to2

npy_file = sys.argv[1]
clustfile = sys.argv[2]
cutoff = float(sys.argv[3])
output_subclust = sys.argv[4]
output_superclust = sys.argv[5]

outp_npy = None
if len(sys.argv) > 6:
    outp_npy = sys.argv[6]

def read_clustfile(clustfile):
    clust = []
    for l in open(clustfile):
        ll = l.split()[3:]
        clust.append([int(v)-1 for v in ll])
    return clust

def write_clustfile(clust, clustfile):
    cf = open(clustfile, "w")
    for cnr, c in enumerate(clust):
        print("Cluster %d ->"%(cnr+1), end=' ', file=cf)
        # !!! changed to start at 1 !!! 18/12/17
        for cc in c: print(cc + 1, end=' ', file=cf)
        print("", file=cf)

rootclusters = read_clustfile(clustfile)

superclust = []
subclust = []
maxstruc = 100000

coors = npy3to2(np.load(npy_file))
nstruc, natom = coors.shape[:2]
lim = cutoff * cutoff * natom

for r in rootclusters:
    if len(r) == 1:
        subclust.append(r)
        superclust.append([len(subclust)])
        continue
    leafclust = []
    clust_struc = coors[r]
    d = squareform(pdist(clust_struc, 'sqeuclidean'))
    d2 = d<lim
    clustered = 0
    while clustered < len(r):
        neigh = d2.sum(axis=0)
        heart = neigh.argmax()
        leaf = np.where(d2[heart])[0]
        for cs in leaf:
            d2[cs,:] = False
            d2[:, cs] = False
        leaf = [heart+1] + [v+1 for v in leaf if v != heart]
        leafclust.append(leaf)
        clustered += len(leaf)
    mapped_root = []
    for leaf in leafclust:
        mapped_leaf = [r[n-1] for n in leaf]
        subclust.append(mapped_leaf)
        mapped_root.append(len(subclust))
    superclust.append(mapped_root)


s = [s[0] for s in subclust]
subclustnpy = coors[s]
write_clustfile(superclust, output_superclust)
write_clustfile(subclust, output_subclust)
if outp_npy is not None:
    np.save(outp_npy,subclustnpy)
