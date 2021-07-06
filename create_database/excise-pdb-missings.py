#!/usr/bin/env python3
# Author: Isaure Chauvot de Beauchene (CNRS)

"""
Removes XXXX and anything that is not RA,RC,RG,RU
Create excise.json describing missing atoms in residues:
    1 = missing atom(s) in phostphate
    2 = missing atom(s) in sugar
    4 = missing atom(s) in base
    ----------
    3 = missing atom(s) in sugar & phosph
    5 = missing atom(s) in base & ph
    6 = missing atom(s) in base & sugar
    3, 6, 7 <=> deleted (too many missing atoms)

!!! ADDED 01-03-2018: If P but no PO1/PO2: remove P
"""

import sys, os, json, copy, argparse
from collections import defaultdict

def testpathjson(jsonfile):
    if os.path.exists(jsonfile):
            #print >> sys.stderr, dictname +" = json.load(open('"+jsonfile+"'))"
        dictname = json.load(open(jsonfile))
    else:
        dictname = {}
        #exec( dictname + "= {}" )
    return dictname

def del_chain(struc, c):
    print('del %s %s'%(struc, c), file=sys.stderr)
    global out
    d = out[struc]
    cc="chain_%s"%c
    for k in d.keys():
        if isinstance(d[k], dict):
            if cc in d[k].keys():
                del d[k][cc]
    if cc in d["nachains"]:
        d["nachains"].remove(cc)
    if len(d["nachains"]) == 0:
        del out[struc]
        print('deleted struc %s'%struc, file=sys.stderr)
        return 1
    return 0

d={}
def process(Nmissings, missing, resi, lines, seq):
    newL = []
    if missing in [2, 3, 6] and Nmissings > 5:
        missing = 7
    if missing != 7:
        for ll in lines:
            newL.append(ll)
    if missing == 7:
        if len(seq):
            seq = seq[:-1]
        else:
            seq = []
    #if missing != 0:
    #print("%i missing atoms in res %s, code %i"%(Nmissings, resi, missing), file=sys.stderr)
    return seq, newL, missing

parser = argparse.ArgumentParser(description=__doc__,
                          formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("inpdir",help="cleanPDB/")
parser.add_argument("inpjson",help="chainsmodels.json")
parser.add_argument("outjson",help="excise.json")
parser.add_argument("outlist",help="excised.list")
parser.add_argument("na",help="dna or rna")
parser.add_argument("--delete", help="remove from outp entries not in inp", action="store_true")
parser.add_argument("--replace", help="list of entries to replace in outp")
parser.add_argument("--list", help="force processing all and only entries in list")
parser.add_argument("--checkinp", help="remove entries for which input pdb doesn't exist", action="store_true")

args = parser.parse_args()

chainsmodels = json.load(open(args.inpjson))

out_init = testpathjson(args.outjson)
if args.delete:
    out = {}
else:
    out = out_init

replace = []
if args.replace:
    replaces = [l.split()[0] for l in open(args.replace)]

resn_list = ["RA","RC","RG","RU"]
if args.na=='dna':
    resn_list = ["DA","DC","DG","DT"]

ph = ["P","O1P","O2P","OP1", "OP2","O5'"]
sug = ["C5'", "C4'", "C3'", "C2'", "C1'", "O3'", "O2'", "O4'"]
base = ["N1","C2","N2","O2","N3","N4","O4","C4","C5","N6","C6","O6","N7","C7","C8","N9"]
code = { 1:ph, 2:sug, 4:base}

missingcode = defaultdict(lambda: 0)
for a in [1, 2, 4]:
    for b in code[a]:
        missingcode[b] = a

if args.list is not None:
    list = [l.split()[0] for l in open(args.list).readlines()]
else:
    list = sorted(chainsmodels.keys())

for struc in list:
    if struc in out_init and not struc in replace and args.list is None \
        and not args.checkinp:
        out[struc] = out_init[struc]
        continue
    print(struc, file=sys.stderr)
    try:
        d = chainsmodels[struc]
        chains = copy.deepcopy(d["nachains"])
        if len(chains) == 0: continue
        assert d["Nmodels"]
        out[struc] = d
        d["missing_atoms"] = {"chain_%s"%c:{} for c in d["nachains"]}
        d["breaks"] = {"chain_%s"%c:[] for c in d["nachains"]}
        delstruc = 0
        for c in chains:
            cc="chain_%s"%c
            for m in range(1, d["Nmodels"]+1):
                inp = "%s/%s%s-%i-iniparse.pdb"%(args.inpdir, struc, c, m)
                if not os.path.exists(inp):
                    print('%s does not exist'%inp, file=sys.stderr)
                    if args.checkinp:
                        delstruc = del_chain(struc, c)
                    continue
                outfile = "%s/%s%s-%i-iniparse-excise.pdb"%(args.inpdir, struc, c, m)
                if os.path.exists(outfile):
                    print('%s already exists'%outfile, file=sys.stderr)
                    if m > 1 or ("sequence" in d and c in d["sequence"]):
                        continue
                #print("processing %s %s"%(struc, c), file=sys.stderr)
                L, lines, seq = [], [], []
                resi = -999
                Nmissings, missing = 0, 0
                for l in open(inp, 'r'):
                    if not l.startswith("ATOM") and not l.startswith("HETATM"):
                        L.append(l)
                        continue
                    resname = l[17:20].strip()
                    resid = l[22:26].strip()
                    atomname = l[13:16].strip()
                    if atomname == "O5T": continue
                    if resname not in resn_list: continue
                    if resid != resi:
                        seq, newL, missings = process(Nmissings, missing, resi, lines, seq)
                        L += newL
                        if missings != 0:
                            d["missing_atoms"][cc]["res_%s"%resi] = missing
                            if missings == 7:
                                d["breaks"][cc].append(resi)
                        seq.append(resname[1])
                        lines = []
                        resi = resid
                        Nmissings = 0
                        missing = 0
                    if l[31].strip() == "X": #missing atom from --manual mode
                        #if atomname == "O5T": continue
                        Nmissings += 1
                        missing = missingcode[atomname] | missing
                        continue
                    lines.append(l)
                #### parse last residue
                seq, newL, missings = process(Nmissings, missing, resi, lines, seq)
                L += newL
                if missings != 0:
                    d["missing_atoms"][cc]["res_%s"%resi] = missing
                    if missings == 7:
                        d["breaks"][cc].append(resi)
                ########################
                outf = open(outfile, "w")
                for l in L:
                    print(l, end='', file = outf)
                outf.close()
                ########################
                if m == 1:
                    if 'sequence' not in d.keys():
                        d["sequence"] = {}
                    d["sequence"]["chain_%s"%c] = "".join(seq)
            if delstruc:
                continue
    except:
        continue

json.dump(out, open(args.outjson, "w"), indent = 2, sort_keys = True)

outlist = open(args.outlist, 'w')
for struct in out:
    d = out[struct]
    for c in d["nachains"]:
        for m in range(1, d["Nmodels"]+1):
            outfile = "%s/%s%s-%i-iniparse-excise.pdb"%(args.inpdir, struct, c, m)
            print(outfile, file = outlist)
outlist.close()
