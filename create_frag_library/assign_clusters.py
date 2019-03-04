#!/usr/bin/env python3
# Author: Isaure Chauvot de Beauchene (CNRS)

import numpy as np
import sys, argparse, json

########################
parser =argparse.ArgumentParser(description=__doc__,
formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('inp', help="fragments.json")
parser.add_argument('outp', help="fragments_clust.json")
parser.add_argument('na', help="rna or dna")
#parser.add_argument('outp', help="fragments_clust.json")
parser.add_argument('--clustfiles', nargs="+")
parser.add_argument('--clustnames', nargs="+")
args = parser.parse_args()
########################

inp = json.load(open(args.inp))
outp = {}
na = args.na

assert len(args.clustnames) == len(args.clustfiles)

def get_clust(filename):
    print(filename, file=sys.stderr)
    ll = [l.split()[3:] for l in open(filename)]
    clusters = [ l if len(l)>1 else [l[0]] for l in ll]
    return clusters

s = ["C", "A"]
if na == 'dna':
    s = ["G", "T"]

mutations = {'G': 'A', 'U': 'C', 'T': 'C'}
mutpattern = [[a, b, c] for a in [0, 1] for b in [0, 1] for c in [0, 1] ]

def mutate(seq, pattern):
    motif = []
    for i, p in enumerate(pattern):
        m = seq[i]
        if p == 1:
            m = mutations[seq[i]]
        motif.append(m)
    return "".join(motif)

count = 0
#for (a, b, c) in [(a, b, c) for a in s for b in s for c in s]:
for (a, b, c) in [("A", "A", "A")]:
    motif = a+b+c
    print(motif, file=sys.stderr)
    for frag in inp[motif]:
        outp[motif] = inp[motif]
        for name in args.clustnames:
            outp[motif][frag]['%s_center'%name] = False
            outp[motif][frag][name] = 0
    dr = get_clust(motif + "-" + args.clustfiles[0])
    clust1 =  get_clust(motif + "-" + args.clustfiles[1])
    clust2 =  get_clust(motif + "-" + args.clustfiles[2])
    drname, name1, name2 = args.clustnames
    for nd, d in enumerate(dr):
        print(d)
        center = d[0]
        outp[motif][str(center)]['%s_center'%drname] = True
        for frag in d:
            outp[motif][frag][drname] = nd+1
    for nc, cl in enumerate(clust1):
        for c in cl:
            d = dr[int(c)-1]
            center = dr[int(cl[0])-1][0]
            outp[motif][str(center)]['%s_center'%name1] = True
            for frag in d:
                outp[motif][frag][name1] = nc+1
    for nc2, cl2 in enumerate(clust2):
        center_c1 = clust1[ int(cl2[0])-1 ][0]
        center = dr[int(center_c1)-1][0]
        outp[motif][str(center)]['%s_center'%name2] = True
        for c2 in cl2:
            cl1 = clust1[int(c2)-1]
            for c in cl1:
                d = dr[int(c)-1]
                for frag in d:
                    outp[motif][frag][name2] = nc2+1

json.dump(outp, open(args.outp,'w'), indent = 2, sort_keys = True)
