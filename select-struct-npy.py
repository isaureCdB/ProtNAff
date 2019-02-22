#!/usr/bin/env python3
#Copyright 2018 Isaure Chauvot de Beauchene

import numpy as np
import sys, argparse

########################
parser =argparse.ArgumentParser(description=__doc__,
formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('npyfile')
parser.add_argument('outpfile')
parser.add_argument('--atom', '--atoms', type=int, nargs='+', help="indices of atoms to select")
parser.add_argument('--structure', '--structures', '--struct', dest='structure',
                    type=int, nargs='+', help="indices of structures to select")
parser.add_argument('--atname', '--atnames','--atomname', dest='atname', type=str,
                    nargs='+', help="name of atoms to select")
parser.add_argument('--top', help="number of top structures to select", type=int)
parser.add_argument('--template', help="pdb template to select atom indices")

args = parser.parse_args()
########################

npy = np.load(args.npyfile)
reshape = False
if len(npy.shape) == 3:
    assert npy.shape[2] == 3
else:
    assert len(npy.shape) == 2 and npy.shape[1]%3 == 0
    npy = npy.reshape(npy.shape[0], int(npy.shape[1]/3),3)
    reshape = True

if args.structure:
    sel = [ int(i)-1 for i in args.structure]
    npy = npy[sel, : , :]

sel = range(npy.shape[1])
if args.atom:
    sel = [ int(i)-1 for i in args.atom]
elif args.atname:
    if args.template is None:
        print("atname option requires a pdb template (--template)", file=sys.stderr)
        raise
    pdb = [l for l in open(args.template).readlines() if l.startswith('ATOM') ]
    names = set(args.atname)
    sel = [ lnr for lnr, l in enumerate(pdb) if l[13:16].strip() in names ]

npy = npy[:, sel , :]

if reshape:
    npy = npy.reshape(npy.shape[0],npy.shape[1]*3)

np.save(args.outpfile, npy)
