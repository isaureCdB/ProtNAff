#!/usr/bin/env python3
import sys

ll = open(sys.argv[1]).readlines()
clusters = [ [i for i in l.split()[3:]] for l in ll]
size = [len(c) for c in clusters]

print(ll[0][:-1], end="\n")
for nc in range(len(clusters)):
    if nc==0: continue
    if sum(size[nc:]) > size[0]:
        print(ll[nc][:-1], end="\n")
    else:
        print("cluster %i => "%(nc+1), end="")
        for clust in clusters[nc:]:
            for c in clust:
                print(c, end=" ")
        break
