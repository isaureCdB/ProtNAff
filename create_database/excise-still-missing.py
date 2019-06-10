#!/usr/bin/env python3
# Author: Isaure Chauvot de Beauchene (CNRS)

"""
Removes XXXX and anything that is not RA,RC,RG,RU
Add to excise.json the missing atoms (usually 5PH) that could not be completed
by nalib because all conformers in nalib clash:
    7 <=> deleted (too many missing atoms)
"""

import sys, os, json, copy
from collections import defaultdict

d={}
def process(missing, mapped_resi, L, lines):
    global d, c
    d['missing_atoms'][c][mapped_resi] = missing
    if missing == 0:
        for ll in lines:
            L.append(ll)
    else:
        d['breaks'][c].append(mapped_resi)
    if not missing:
        for l in lines:
            L.append(l)
    return L

#print("excise-pdb-missings", file=sys.stderr)
inpdir = sys.argv[1]    #cleanPDB
inplist = sys.argv[2]   #still_missing.list
outplist = sys.argv[3]  #clean-iniparse-aa.list (cleanPDB/xxxx(|[A-F])[A-Z]-(|[12])([0-9])-iniparse-aa.pdb)
jsonfile = sys.argv[4]  #excise.json
na = sys.argv[5]        # dna/rna

js = json.load(open(jsonfile))

for filename in [l.strip() for l in open(inplist)]:
    name = filename.split('/')[-1].split('-')
    struc = name[0][:-1]    #structure name
    c = name[0][-1]         #chain
    m = name[1]             #model
    d = js[struc]
    mapping = js['mapping'][c]
    if not os.path.exists(filename):
        print('ERROR: %s does not exist'%inp, file=sys.stderr)
        raise
    target = "%s/%s%s-%i-iniparse-aa.pdb"%(inpdir, struc, c, m)
    if os.path.exists(target):
        print('%s already exists'%target, file=sys.stderr)
        continue
    print("processing %s %s"%(struc, c), file=sys.stderr)
    L, lines, seq = [], [], []
    resi = -999
    for l in open(inp, 'r'):
        resname = l[17:20].strip()
        resid = l[22:26].strip()
        atomname = l[13:16].strip()
        if resid != resi:
            d['missing_atoms'][c][resi] = missing
            L = process(missing, mapping[resi], L, lines)
            lines = []
            resi = resid
            missing = 0
        if l[31].strip() == "X": #missing atom from --manual mode
            missing = 7
            continue
        lines.append(l)
    #parse last residue
    L = process(missing, mapping[resi], L, lines)

    outf = open(target, "w")
    print(target)
    for l in L:
        print(l, end='', file = outf)
    outf.close()

json.dump(js, open(jsonfile, "w"), indent = 2, sort_keys = True)
