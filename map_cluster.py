#!/usr/bin/env python3

import sys

clust = sys.argv[1] #$m-dr0.2r-clust1.0
indices = [ int(l) for l in open(sys.argv[2]).readlines()] #$m-dr0.2.list
mapping = { str(ni+1) : i for ni, i in enumerate(indices)}

for nl, l in enumerate(open(clust).readlines()):
    old = l.split()[3:]
    print('cluster %s =>'%str(nl+1), end=' ')
    for i in old:
        ii = mapping[str(i)]
        print(ii, end=' ')
    print('')
