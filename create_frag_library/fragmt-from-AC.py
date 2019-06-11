#!/usr/bin/python3
# Copyright 2017 - 2018 Isaure Chauvot de Beauchene (CNRS)

import sys
import os, json, argparse
import numpy as np

def get_coor(l):
    if len(l)<54:
        raise
    return [float(l[30:38]), float(l[38:46]),float(l[46:54])]

def check_breaks(index, d, cc):  #19/09/18
    breaks = d['breaks'][cc]
    if breaks is None:
        return breaks, 0
    print(breaks)
    if (index+3) in breaks or (index+2) in breaks:
        return breaks, 1
    else:
        return breaks, 0

#################################################################################
a = argparse.ArgumentParser(prog="fragmt-from-AC.py")
a.add_argument("inp")            # structures.json
a.add_argument("frags")          # fragments.json
a.add_argument("na")             # dna / rna
a.add_argument("listofseq")      # motifs.list
a.add_argument("directory")      # cleanPDB
args = a.parse_args()
#################################################################################
directory = args.directory
bases = ['C', 'A']
mutations = {'A': 'A',
            'G': 'A',
            'U': 'C',
            'C': 'C',
            'T': 'C'
            }

listofseq = open(args.listofseq, 'w')
sequences = []
for a in bases:
    for b in bases:
        for c in bases:
            sequences.append(a+b+c)
            print(a+b+c, file = listofseq)
listofseq.close()

coor_bases = {}
count = {}
all_frag = {}
for s in sequences:
    coor_bases[s] = []
    count[s] = 1
    all_frag[s] = {}
prog=0
inp = json.load(open(args.inp))
for struct in sorted(inp.keys()):
    print(struct, file=sys.stderr)
    d = inp[struct]
    print(d['missing_atoms'], file=sys.stderr)
    for m in range(1, d['Nmodels']+1):
        mm = "model_%i"%m
        for c in d['nachains']:
            cc = "chain_"+c
            pdb = '%s/%s%s-%i-iniparse-aa-AC.pdb'%(directory, struct,c,m)
            if not os.path.exists(pdb):
                print('%s does not exist'%pdb)
                continue
            seq = d['sequence'][cc]
            mu = [mutations[i] for i in seq]
            museq = "".join(mu)
            mapp = d['mapping'][cc]
            missings = d['missing_atoms'][cc]
            miss = []
            #print(missings, file=sys.stderr)
            resid, res, coors, seqpdb = 0, "", [], []
            for l in open(pdb,'r').readlines():
                coor = get_coor(l)
                if l[13]=='P':
                    if res in missings:
                        miss.append(missings[res])
                    seqpdb.append(l[19])
                    resid += 1
                    res = 'res_%i'%resid
                    coors.append([])
                coors[-1].append(coor)
            for i in range(len(coors)-2):
                try:
                    breaks, val = check_breaks(i, d, cc)
                except:
                    print(struct, c, len(coors), coors[-1][-1])
                    raise
                if val: continue
                s = museq[i:i+3]
                ss = "".join(seqpdb[i:i+3])
                assert s == ss, (struct, m, c, i, s, ss)
                oriseq = seq[i:i+3]
                coordinates = [c for cc in coors[i:i+3] for c in cc]
                coor_bases[s].append(coordinates)
                ind = str(count[s])
                all_frag[s][ind] = {}
                all_frag[s][ind]['structure'] = struct
                all_frag[s][ind]['chain'] = c
                all_frag[s][ind]['model'] = m
                all_frag[s][ind]['indices'] = (i+1, i+2, i+3)
                all_frag[s][ind]['resid'] = (mapp[str(i+1)], mapp[str(i+2)], mapp[str(i+3)])
                all_frag[s][ind]['seq'] = oriseq
                all_frag[s][ind]['missing_atoms'] = miss[i:i+3]
                count[s]+=1

json.dump(all_frag, open(args.frags, 'w'), indent = 2)

for s in sequences:
    if not len(coor_bases[s]): continue
    nat = len(coor_bases[s][0])
    nstruc = len(coor_bases[s])
    nn = np.array(coor_bases[s])
    print("%s: %i atoms, %i structures"%(s, nat, nstruc), file=sys.stderr)
    np.save('%s-all-aa.npy'%s, np.array(coor_bases[s]))
