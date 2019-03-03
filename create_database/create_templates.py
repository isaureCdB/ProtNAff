#!/usr/bin/python3
import sys, os


outdir = sys.argv[1]
na = sys.argv[2]

dirname = os.path.dirname(os.path.abspath(__file__))
indir = "%s/%slib"%(dirname,na)

def create_template(seq):
    out = open('%s/%s.pdb'%(outdir, seq), 'w')
    for i, s in enumerate([s for s in seq]):
        for ll in open('%s/%s.pdb'%(indir, s), 'r'):
            print("%s%i%s"%(ll[:25], i+1, ll[26:-2]), file=out)
    out.close()

x = "U"
if na == 'dna':
    x = "T"

s = ['A','C','G',x]
sequences = [a+b+c for a in s for b in s for c in s ]

for seq in sequences:
    create_template(seq)
