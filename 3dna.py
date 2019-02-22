#!/usr/bin/env python3
import sys, os, json

'''
TODO: Fix bug: nucleotides with weird angles are discarded by 3dna,
and replaced by a "&" in the sequence ("bseq").
'''

chainsmodels = json.load(open(sys.argv[1])) #chainsmodels.json
script = sys.argv[2]    #x3dna.sh?

def x3dna(struct, chains, nmodels):
    for m in range(nmodels):
        inp = "/tmp/%s-%i-3dna.pdb"%(struct,m+1)
        os.system("grep ATOM interface/%s-%i.pdb > %s "%(struct,m+1,inp))
        for c in chains:
            os.system("echo TER >> %s"%inp)
            os.system("cat cleanPDB/%s%s-%i-iniparse-aa.pdb >> %s"%(struct,c,m+1,inp))
        os.system("%s %s %i"%(script, struct, m+1))

for struct in sorted(chainsmodels.keys()):
    d = chainsmodels[struct]
    x ="3dna/%s-1-dssr.json"%struct
    if os.path.exists(x) and os.path.getsize(x) > 0:
        continue
    print(struct, file=sys.stderr)
    Nmodels = d['Nmodels']
    chains = d['nachains']
    x3dna(struct, chains, Nmodels)
