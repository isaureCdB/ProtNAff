#!/usr/bin/env python3
# Author: Isaure Chauvot de Beauchene (CNRS)

"""
Modify rna/dna pdb writen by interface_pdb.py
All chain identifiers are removed from the chain
Nucleotide-like cofactors and incompatible bases had been excised
Within nucleic acid chains, PSU/7AT residues are also canonized (special since it involves atom re-mapping)
"""

import sys, os, json, copy, argparse

def testpathjson(jsonfile):
    if os.path.exists(jsonfile):
        #print >> sys.stderr, dictname +" = json.load(open('"+jsonfile+"'))"
        dictname = json.load(open(jsonfile))
    else:
        dictname = {}
        print('%s does not exist'%jsonfile, file=sys.stderr)
        #exec( dictname + "= {}" )
    return dictname

def del_chain(struc, c):
    print('del %s %s'%(struc, c), file=sys.stderr)
    global js
    d = js[struc]
    cc = "chain_"+c
    for k in d.keys():
        if isinstance(d[k], dict):
            if cc in d[k].keys():
                del d[k][cc]
    d['nachains'].remove(cc)
    if len(d['nachains']) == 0:
        del js[struc]

def process_struc(struc, out, args):
    print("process %s"%struc, file=sys.stderr)
    d = js[struc]
    out[struc] = d
    d['canonized'] = {}
    mut = []
    chains = copy.deepcopy(d['nachains'])
    for c in chains:
        cc = "chain_"+c
        for m in range(1, d['Nmodels']+1):
            outfile = "%s/%s%s-%i.pdb"%(args.outdir, struc, c, m)
            if os.path.exists(outfile):
                continue
            #print("processing %s %s"%(struc, c), file=sys.stderr)
            inpfile = "%s/%s%s-%i.pdb"%(args.inpdir, struc, c, m)
            if not os.path.exists(inpfile):
                print("%s does not exist"%inpfile, file=sys.stderr)
                #del_chain(struc, c)
                continue
            inf = open(inpfile,'r')
            prev = 0
            L = []
            for l in inf:
                if l.split()[0] not in ("ATOM", "HETATM"): continue
                res = int(l[23:26])
                resn = l[17:20].strip()
                atomname = l[12:16].strip()
                if resn == "23": #A23, C23, etc. remove oxygen atoms
                    if atomname in ("O2'", "O3'"): continue
                elif resn in modif.keys():
                    base = modif[resn][0][args.na]
                    dict_at = modif[resn][1]
                    if atomname.strip() in dict_at:
                        atomname = dict_at[atomname.strip()]
                    l = "%s%4.4s   %s%s" % (l[:12], atomname, base, l[20:])
                elif resn == "H2U" and atomname in ["C2", "C4", "C5", "C6", "N3"]: continue
                elif resn == "CMU" and atomname in ["C5", "C6", "C5M"]: continue
                if resn in modified or resn in modif.keys():
                    if cc not in d['canonized']:
                        d['canonized'][cc] = set()
                    d['canonized'][cc].add(res)
                prev = res
                L.append(l)
            outf = open(outfile, 'w')
            for l in L:
                print(l[:16] + " " + l[17:-1], file=outf)
            outf.close()
        if struc in js.keys():
            if cc in d['canonized']:
                d['canonized'][cc] = list(d['canonized'][cc])
    return out

parser = argparse.ArgumentParser(description=__doc__,
                          formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("inpdir",help="chainsmodels")
parser.add_argument("outdir",help="cleanPDB")
parser.add_argument("outlist",help="outp: cleanPDB.list")
parser.add_argument("mutatelist",help="list of modified residue and their equivalent")
parser.add_argument("inpjson",help="chainsmodels.json")
parser.add_argument("outjson",help="clean_rna.json")
parser.add_argument("na",help="dna or rna")
parser.add_argument("--delete", help="remove from outp entries not in inp", action="store_true")
parser.add_argument("--replace", help="replace entries already in outp", action="store_true")
parser.add_argument("--files", help="process entries for which output processed \
                            PDB file does not exist (even if the entry already \
                            exists in the output json file)", action="store_true")
parser.add_argument("--subset", help="list of entries to process and replace")


args = parser.parse_args()

js = json.load(open(args.inpjson))
modified = [l.split()[0] for l in open(args.mutatelist).readlines() if len(l.strip())]

out_init = testpathjson(args.outjson)
if args.delete:
    out = {}
else:
    out = out_init

base = ["N1","C2","N2","O2","N3","N4","C4","O4","C5","N6","C6","O6","C7","N7","C8","N9"]
### ADD OTHER MUTATIONS
modif = {
    'BZG': [{ 'rna':'G', 'dna':'G'},
            {"OS'": "O4'",
            "CP'": "C2'",
            "CT'": "C1'",
            "NE": "N9",
            "CO": "C8",
            "NN": "N7",
            "CM": "C5",
            "CF": "C4",
            "NG": "N3",
            "CH": "C2",
            "NI": "N2",
            "NJ": "N1",
            "CK": "C6",
            "OL": "O6",
            }
            ],
    'A3P': [    {'rna':'A', 'dna':'A'},
                {   "P2" : "P",
                    "O4P" : "O1P",
                    "O5P" : "O2P",
                },
            ],
    'AD2': [    {'rna':'A', 'dna':'A'},
                {   "P1" : "P",
                },
            ],
    'AT7': [    {'rna':'A', 'dna':'A'},
                {   "C7": "N7",
                    "N8": "C8",
                }
            ],
    'PSU': [    {'rna':'U', 'dna':'T'},
                {   "N1": "C5",
                    "C5": "N1",
                    "C2": "C4",
                    "C4": "C2",
                    "O4": "O2",
                    "O2": "O4",
                }
            ],
    }

if args.subset:
    todo = [l.split()[0] for l in open(args.subset).readlines()]
else:
    todo = sorted(js.keys())

for struc in todo:
    if struc in out_init and not args.replace:
        chainx = out_init[struc]["nachains"][0]
        outx = "%s/%s-1.pdb"%(args.outdir,struc)
        if not os.path.exists(outx) and args.files:
            out = process_struc(struc, out, args)
            continue
        if args.delete:  # do not replace already existing entries
            out[struc] = out_init[struc]
        print("skipping %s"%struc, file=sys.stderr)
        continue
    out = process_struc(struc, out, args)

outlist = open(args.outlist, 'w')
for struct in out:
    d = out[struct]
    for c in d['nachains']:
        for m in range(1, d['Nmodels']+1):
            outfile = "%s/%s%s-%i.pdb"%(args.outdir, struct, c, m)
            print(outfile, file=outlist)
outlist.close()

json.dump(out, open(args.outjson, "w"), indent = 4)
print('%s dumped'%args.outjson, file=sys.stderr)
