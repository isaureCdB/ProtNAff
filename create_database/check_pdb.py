#!/usr/bin/env python3
# Author: Isaure Chauvot de Beauchene (CNRS)

from Bio.PDB import *
import numpy as np
import sys
import os, json
from collections import defaultdict
import argparse

'''
Check resolution, presence of DNA/RNA, if all models have the same atoms...
If alternative atomic positions => subdivide in several chains.
'''

def select_alternate_conf(pdb):
    struc = pdb.split('.pdb')[0]
    alternames = ["A", "B", "C", "D", "E", "F"]
    ll = [l for l in open(pdb).readlines()]
    alternates = set([])
    for l in ll:
        if len(l)<54:
            continue
        if l[16] in alternames:
            alternates.add(l[16])

    outp = {}
    for a in alternates:
        outp[a] = open("%s%s.pdb"%(struc, a),"w")

    for l in ll:
        if not l.startswith("ATOM"):
            for o in outp:
                print(l, file=outp[o], end="")
        elif l[16]==' ':
            for o in outp:
                print(l, file=outp[o], end="")
        else:
            newl = l[:16] + " " + l[17:]
            print(newl, file=outp[l[16]], end="")

    for o in outp:
        outp[o].close()

class MySelect(Select):
    def __init__(self, model_toselect, res_toselect):
        self.res_toselect = set([id(r) for r in res_toselect])
        self.model_toselect = model_toselect
        Select.__init__(self)
    def accept_residue(self, residue):
        if id(residue) in self.res_toselect:
            return 1
        else:
            return 0
    def accept_model(self, model):
        if model is self.model_toselect:
            return 1
        else:
            return 0
    def accept_atom(self, atom):
        if atom.get_id()[0]=="H":
            return 0
        else:
            return 1

parser = argparse.ArgumentParser(description=__doc__,
                          formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("inpdir",help="brutPDBs")
parser.add_argument("inplist",help="pdbcodes.txt")
parser.add_argument("corrupted",help="corrupted_pdb_files.list")
parser.add_argument("tofix",help="tofix.list")
parser.add_argument("checked",help="checked.list")
parser.add_argument("splitted",help="splitted.list")
parser.add_argument("mutatelist",help="[d/r]nalib/mutate.list")
parser.add_argument("na",help="dna or rna")
args = parser.parse_args()

parser = PDBParser()

modified = set([l.split()[0] for l in open(args.mutatelist)])
canonical = set(["DT", "DC", "DG", "DA", "C", "T", "A", "G", "PSU", "7AT",
                 "RU", "RC", "RG", "RA", "U"])
included = canonical | modified

def testpath(filename):
    if os.path.exists(filename):
        ll = [ l.strip() for l in open(filename).readlines()]
        return open(filename, 'a'), ll
    else:
        return open(filename, 'w'), []

corrupted, lcorrupted = testpath(args.corrupted)
tofix, ltofix = testpath(args.tofix)
checked, lchecked = testpath(args.checked)
splitted, lsplitted = testpath(args.splitted)

done = lcorrupted + ltofix + lchecked + lsplitted

def check_models_chains(struc):
    chains = []
    for model in structure:
        chains.append([])
        for chain in model:
            chains[-1].append(chain.id)
        if set(chains[-1]) != set(chains[0]):
            return 1
    return 0

def check_alternate(filename, nachains):
    global inpdir, included
    alternates = set()
    for l in open(filename):
        if not l.startswith("ATOM"):
            continue
        if l[21] in nachains:
            if l[17:20].strip() in included and l[16] in ['A','B','C','D','E','F']:
                alternates.add(l[16])
    return alternates

for l in open(args.inplist):
    try:
        struc = l.strip()
        filename = "%s/%s.pdb"%(args.inpdir, struc)
        if struc in done:
            print('%s done'%struc)
            continue
        if not os.path.exists(filename) or os.stat(filename).st_size==0:
            print("ERROR: %s/%s.pdb does not exist"%(args.inpdir, struc))
            continue
        print('-------------------------------' + struc, file=sys.stderr)
        structure = parser.get_structure(struc, filename)
        if structure.header['structure_method'] == "x-ray diffraction":
            if float(structure.header['resolution']) > 3.5:
                print("resolution > 3.5A", file=sys.stderr)
                print(struc, file = corrupted)
                continue
        if check_models_chains(structure):
            print("%s has different MODELS => to fix"%struc)
            print(struc, file = tofix)
            continue
        nachains = []
        for chain in structure[0]:
            c = chain.get_id()[0]
            countrna = 0
            for r in chain:
                if r.has_id("O2'") and args.na == 'rna' :
                    countrna += 1
                if r.has_id("O3'") and args.na == 'dna' and not r.has_id("O2'") :
                    countrna += 1
                if countrna == 3:
                    nachains.append(c)
                    break
        if not len(nachains):
            print("%s contain no %s chains"%(struc, args.na))
            print(struc, file = corrupted)
            continue
        alternates = check_alternate(filename, nachains)
        if alternates:
            print("%s has alternate atoms => splited"%struc)
            select_alternate_conf("%s/%s.pdb"%(args.inpdir, struc))
            print(struc, file = splitted)
            for a in alternates:
                print('%s%s'%(struc, a), file = checked)
        else:
            print(struc, file = checked)
    except:
        print(struc, file = tofix)
        print("##################################", file=sys.stderr)
        print("unknown error in %s"%struc, file=sys.stderr)
        print("##################################", file=sys.stderr)

checked.close()
tofix.close()
corrupted.close()
