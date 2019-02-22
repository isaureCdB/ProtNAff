#!/usr/bin/python3
# Copyright 2017 - 2018 Isaure Chauvot de Beauchene (CNRS)

import sys
import os, json, argparse
import numpy as np

def run(command):
    sys.stdout.flush()
    os.system(command)

def pp(*x):
    for i in x[:-1]:
        print(i, file=sys.stderr, end=' ')
    print(x[-1])

'''
ATOM     63  P    DG     4       8.352 -38.065-100.532   90   0.780 0 1.00
'''
def get_coor(l):
    if len(l)<54:
        raise
    return [float(l[30:38]), float(l[38:46]),float(l[46:54])]

def check_breaks(index, d, cc):  #19/09/18
    breaks = d['breaks'][cc]
    if breaks is None:
        return breaks, 0
    #miss = d['missing_atoms'][cc]
    print(breaks)
    if (index+3) in breaks or (index+2) in breaks:
        return breaks, 1
    else:
        return breaks, 0

def OLD_check_breaks(index, d, cc): #19/09/18
    breaks = d['breaks'][cc]
    try:
        mapp = d['mapping'][cc]
    except:
        print((d.keys(), cc))
        raise
    miss = d['missing_atoms'][cc]
    if mapp[str(index+3)] in breaks or mapp[str(index+2)] in breaks:
        return 1
    else:
        return 0

def check_breaks_from_coor(pdb, struc, cc ):
    ll = open(pdb,'r').readlines()
    P = np.array([ get_coor[l] for l in ll if l[13]=="P" ][1:])
    O = [ get_coor[l] for l in ll if l[13]=="O3'" ][1:]
    X = zip(O,P)
    ss = [ sum([(b[i]-a[i])**2 for i in range(3)]) for a, b in X]
    broken = [n for n, s in enumerate(ss) if s > 3 ]
    breaks_all[struc][cc] = broken

#$RNAFRAG3/fragmt-from-AC.py x3dna.json fragments_ori.json $na listofseq
#################################################################################
a = argparse.ArgumentParser(prog="fragmt-from-AC.py")
a.add_argument("x3dna")          # x3dna.json
a.add_argument("frags")          # fragments_ori.json
a.add_argument("outp")          # structures.json
a.add_argument("na")             #
a.add_argument("listofseq")      #
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
templates = {}
for s in sequences:
    coor_bases[s] = []
    count[s] = 1
    all_frag[s] = {}
    templates[s] = 0
all_templates = 0
prog=0
x3dna = json.load(open(args.x3dna))
for struct in sorted(x3dna.keys()):
    #for struct in ['5HCH']:
    ################################
    #if struct == '1FLJA': prog=1
    #if not prog: continue
    ################################
    pp(struct)
    d = x3dna[struct]
    d['fragments'] = {}
    for m in range(1, d['Nmodels']+1):
        mm = "model_%i"%m
        d['fragments'][mm] = {}
        for c in d['nachains']:
            cc = "chain_"+c
            pdb = '%s/%s%s-%i-iniparse-aa-AC.pdb'%(directory, struct,c,m)
            if not os.path.exists(pdb):
                print('%s does not exist'%pdb)
                continue
            d['fragments'][mm][cc] = []
            fragments = d['fragments'][mm][cc]
            seq = d['sequence'][cc]
            mu = [mutations[i] for i in seq]
            museq = "".join(mu)
            mapp = d['mapping'][cc]
            res, ll, coors, seqpdb = 0, [], [], []
            for l in open(pdb,'r').readlines():
                coor = get_coor(l)
                if l[13]=='P':
                    seqpdb.append(l[19])
                    res+=1
                    coors.append([])
                    if not all_templates: ll.append([])
                coors[-1].append(coor)
                if not all_templates: ll[-1].append(l)
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
                if not templates[s]:
                    out = open('templates/%s.pdb'%s, 'w')
                    for res in ll[i:i+3]:
                        for l in res:
                            print(l[:-1], file=out),
                    out.close()
                    templates[s] = 1
                    all_templates = min([templates[s] for s in sequences])
                coordinates = [c for cc in coors[i:i+3] for c in cc]
                coor_bases[s].append(coordinates)
                fragments.append((s, count[s]))
                ind = str(count[s])
                all_frag[s][ind] = {}
                all_frag[s][ind]['structure'] = struct
                all_frag[s][ind]['chain'] = c
                all_frag[s][ind]['model'] = m
                all_frag[s][ind]['indices'] = (i+1, i+2, i+3)
                all_frag[s][ind]['resid'] = (mapp[str(i+1)], mapp[str(i+2)], mapp[str(i+3)])
                all_frag[s][ind]['seq'] = oriseq
                count[s]+=1


json.dump(all_frag, open(args.frags, 'w'), indent = 2)
json.dump(x3dna, open(args.outp, 'w'), indent = 2)

for s in sequences:
    if not len(coor_bases[s]): continue
    nat = len(coor_bases[s][0])
    nstruc = len(coor_bases[s])
    nn = np.array(coor_bases[s])
    print("%s: %i atoms, %i structures"%(s, nat, nstruc), file=sys.stderr)
    np.save('%s-all-aa.npy'%s, np.array(coor_bases[s]))
