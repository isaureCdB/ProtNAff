#!/usr/bin/env python3
# Author: Isaure Chauvot de Beauchene (CNRS)

import numpy as np
import sys, argparse

'''
Remove trinucleotides with inter-nucleotide clashes
inputs : NNN.npy nat1 nat2 nat3 outp1 outp2
    outp1 = npy of non-clashing coordinates
    outp2 = list of non-clashing indices

$SCRIPTS/discard_clashing_fragments.py $m-aa-fit.npy 2 $m-aa-fit-noclash.npy $m-aa-fit.noclash

to get the list of clashing indices, use:

    a=`$SCRIPTS/shape.py $m-all-aa.npy |awk -F ',' '{print substr($1, 2)}'`
    seq $a > /tmp/seq
    python3 $SCRIPTS/select-lines.py /tmp/seq $m-aa.noclash --indices --reverse|
    sort -nk1 > $m-aa.clashes
'''

def mindist(at1, at2):
    #print(at1.shape)
    #print(at2.shape)
    dif = at1 - at2
    dif2 = np.sum(dif*dif, axis=3)**0.5
    mdif2 = dif2.min(axis=(1,2))
    return mdif2

nnn = np.load(sys.argv[1]) #$m-aa-fit.npy
nat = [int(i) for i in sys.argv[2:5]]
cutoff = float(sys.argv[5])
outp1 = sys.argv[6] #$m-aa-noclash.npy
outp2 = sys.argv[7] #$m-aa.noclash

out2 = open(outp2, 'w')

assert len(nnn[0]) == sum(nat)

n1 = nnn[:, :nat[0]-1, :] # do not include O3'
n23 = nnn[:, nat[0]: , :]
d1 = mindist(n1[:, :, None], n23[:, None, :])

n1 = nnn[:, :nat[0], :]
n2 = nnn[:, nat[0]+1 : nat[0]+nat[1]-1 , :] # do not include O3' not P
n3 = nnn[:, nat[0]+nat[1]: , :]
d2 = mindist(n1[:, :, None], n2[:, None, :])
d3 = mindist(n3[:, :, None], n2[:, None, :])

ddd = np.minimum(d1,d2,d3)

noclashes = np.where(ddd > cutoff)
print(len(noclashes[0]))

for i in noclashes[0]:
    print(i+1, file=out2)
out2.close()

clashes = np.where( ddd <= cutoff )
print(len(clashes[0]))

with open('clashes', 'w') as outp3:
    for i in clashes[0]:
        print('%i %s'%(i+1, ddd[i]), file=outp3)
outp3.close()

keep = nnn[noclashes]
np.save(outp1, keep)
