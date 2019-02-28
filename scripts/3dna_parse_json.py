#!/usr/bin/env python3
# Author: Isaure Chauvot de Beauchene (CNRS)

import sys, json, os
from collections import defaultdict
import copy
'''
"is_broken" in "nts" of 3dna json indicates a break in the chain
=> could be used instead of distance-based script for break detection
(But use distance-based if the user skips 3dna-analysis)

!!! Compatible with x3dna-dssr version before 2018oct09
'''
#TODO: add bulges and intraloop

print("save the kakapo", file=sys.stderr)

def pp(*x):
    for i in x[:-1]:
        print(i, file=sys.stderr, end=' ')
    print(x[-1])

def parse_nts(j):
    nts = [ j[x].split(".") for x in ["nt1","nt2"] ]
    resint = [ int(x[3]) for x in nts ]
    c = [ x[1] for x in nts ]
    resid = [ x[3] for x in nts ]
    return nts, resid, resint, c

def dump_indiv(f, struct):
    json.dump(f, open("3dna/3dna_%s.json"%struct, "w"), indent = 2, sort_keys=True)

def testpathjson(jsonfile):
    pp(jsonfile)
    if os.path.exists(jsonfile):
            #print >> sys.stderr, dictname +" = json.load(open(""+jsonfile+""))"
        dictname = json.load(open(jsonfile))
    else:
        dictname = {}
        #exec( dictname + "= {}" )
    return dictname

def map_indices(struct, c):
    # convert residues indexing:
    # n1 = original PDB
    # n2 = in "interface_protein" and in xxxX-1.mapping
    # n3 = in 3dna and in xxxX-1iniparse-aa.mapping
    dict_n3to1, dict_n2to3, dict_n2to1 = {}, {}, {}
    for l in open("cleanPDB/%s%s-1.mapping"%(struct, c)):   #outp from aareduce.py in parse_*.sh
        ll = l.split()
        dict_n2to1[ll[1]] = ll[0]
    for l in open("cleanPDB/%s%s-1-iniparse-aa.mapping"%(struct, c)): #outp from aareduce.py in parse_*.sh
        ll = l.split()
        dict_n2to3[ll[0]] = ll[1]
        dict_n3to1[ll[1]] = dict_n2to1[ll[0]]
    #
    nind = [int(i) for i in dict_n3to1]
    firstres, lastres = min(nind), max(nind)
    #
    return firstres, lastres, dict_n3to1, dict_n2to3

def map_missing(d, cc, dict_n2to3):
    if "missing_atoms" in d:
        if cc in d["missing_atoms"]:
            missings = d["missing_atoms"][cc]
            new = { dict_n2to3[i]: missings[i] for i in missings.keys() if i in dict_n2to3.keys() }
            d["missing_atoms"][cc] = new
    return d

def map_interf(d, cc, dict_n3to1):
    if cc not in d["interface_protein"]["model_1"]:
        d["interface_protein"]["model_1"][cc] = []
    interf = d["interface_protein"]["model_1"][cc]  # n1 {B: {"20A": {interface} }}
    dict_n1to3 = { "res_"+v:"res_"+k for k,v in dict_n3to1.items()} # {B: {"20A": "1"}}
    new_interf = {}
    if len(interf) == 0:
        return new_interf
    for i in sorted(interf.keys()):
        if i not in dict_n1to3.keys():
            continue
        new_interf[ dict_n1to3[i] ] =  interf[i]
    return new_interf

def initialise_all(d, js):
    d["mapping"] = {"chain_"+c:{} for c in chains}
    d["breaks"] = {}
    d["ss"] = {"chain_"+c:{} for c in chains}
    #
    keys = ["hbonds","nonPairs","pairs","stems","hairpins","junctions","ssSegments"]
    for k in keys:
        if k not in js:
            js[k] = []
    #
    outputs = ["bptype", "intraRNA_hb", "stacking", "ss", "RNAprot_hb",\
    "DNApairing"] # +"breaks"
    for k in outputs:
        if k not in d:
            d[k] = { "chain_"+c:{} for c in chains}
    return d, js

def initialise_chain(d, c, m, dict_n3to1):
    global codenames
    for l, r in enumerate(dict_n3to1):
        cc, rr = "chain_"+c, "res_"+r
        d["ss"][cc][rr] = ["S", l+1, len(dict_n3to1)]
    return d

def update_intraRNA(d_cc, rr, code, n, value):
    # distinguish hbond of nucl n1 with nucl n2 in [n-2, n-1, other, n+1, n+2]
    if rr not in d_cc :
        d_cc[rr] = {}
    if code not in d_cc[rr]:
        d_cc[rr][code] = {}
    if n not in d_cc[rr][code]:
        d_cc[rr][code][n] = 0
    d_cc[rr][code][n] += value
    return d_cc

def update_stacking(d_cc, rr, n):
    # distinguish hbond of nucl n1 with nucl n2 in [n-2, n-1, other, n+1, n+2]
    if rr not in d_cc :
        d_cc[rr] = {}
    d_cc[rr][n] = 1
    return d_cc

def update_RNAprot(d_cc, rr, code, value):
    if rr not in d_cc :
        d_cc[rr] = {}
    if code not in d_cc[rr]:
        d_cc[rr][code] = 0
    d_cc[rr][code] += value
    return d_cc

def check_breaks(js, c):
    breaks = []
    cc = "chain_"+c
    if "dbn" not in js.keys(): return breaks
    if cc not in js["dbn"].keys(): return breaks
    bseq = js["dbn"][cc]["bseq"].strip()
    nts = [ jj for jj in js["nts"] if jj["chain_name"] == c]
    nt_id = [ int(x["nt_id"].split(".")[3]) for x in nts ]
    for nr, n in enumerate(nt_id):
        if bseq[nr] == "&":
            breaks.append(n)
            if len(bseq) > nr:
                bseq = bseq[:nr] + bseq[nr+1:]
                if len(breaks) > 0 : pp((c, "breaks", breaks))
                return breaks

def EX_check_breaks(js, c): #changed on 18/09/2018
                    # those are breaks in the backbone. Does not account for missing residues
                    # such as those excised by aareduce
                    nts = [ jj for jj in js["nts"] if jj["chain_name"] == c]
                    breaks = [n["nt_resnum"] for n in nts if n["is_broken"]]
                    if len(breaks) > 0 : pp((c, "breaks", breaks))
                    return breaks

def weight_hbond(js):
    # weight the Hbond by its atom-atom distance
    # shorter = more probable
    dist = float(js["distance"])
    donnor = js["donAcc_type"]
    value = 0.5
    if dist < 3.5:
        value = 0.75
        if dist < 3.0:
            value = 1.0
    if donnor == "questionable":
        value *= 0.5
    return value

def get_hbonds(js_hbond, d_RNAprot_hb, d_intraRNA_hb):
    global chains
    for j in js_hbond:
        if j["residue_pair"] == "aa:aa":
            # ignore prot-prot hbonds
            continue
        value = weight_hbond(j)
        # !!! some nt are not in the main chains (e.g mono-nucleotides) !!!
        # Then they are not in residues, as j[x].split("@")[1][2] != "R"
        # Ex of atom description by x3dna: "atom2_id": "N6@.M.RA.8.",
        atoms = [ j[x].split("@") for x in ["atom1_id", "atom2_id"]]
        residues = [ at[1].split(".") for at in atoms]
        nt_indices = [ xr for (xr, x) in enumerate(residues) if x[2][0] == prefix and x[1] in chains]
        res = [ residues[i] for i in nt_indices]
        if len(res) == 0: continue
        if j["residue_pair"] == "nt:nt":
            if len(res) != 2: #hbond with another nucleic acid type
                continue
            # 2 nucl in correct nucleic acid type
            at = [ a[0].split(".")[0] for a in atoms ]
            c = [x[1] for x in res ]
            cc = [ "chain_"+x[1] for x in res ]
            resint = [ int(x[3]) for x in res ]
            rr  = [ "res_%i"%resint[i] for i in [0, 1]]
            resname = [ x[2] for x in res ]
            for (a, b) in [(0, 1), (1, 0)]:
                # get the nucleotide part involved in the hbond
                part = dictcode[at[a]]
                code = codenames[part]
                # default = the two nucl are no neighbors
                n = "other"
                if c[0] == c[1] and c[0] in chains:
                    # neighbors nucl in same chain
                    dist_nucl = resint[b]-resint[a]
                    if abs(dist_nucl) <= 2:
                        n = neighbors[dist_nucl + 2]
                d_intraRNA_hb[cc[a]] = update_intraRNA(d_intraRNA_hb[cc[a]], rr[a], code, n, value)
        else:
            #"nt:aa" or "aa:nt":
            at = [ atoms[i] for i in nt_indices][0]
            cc, rr = "chain_"+res[0][1], "res_"+res[0][3]
            atname = at[0].split(".")[0]
            part = dictcode[atname]
            code = codenames[part]
            d_RNAprot_hb[cc] = update_RNAprot(d_RNAprot_hb[cc], rr, code, value)
    return d_RNAprot_hb, d_intraRNA_hb

def get_nonPairs_hbonds(j_nonPairs, d_intraRNA_hb):
    # the interaction is an H-bond
    # Ex: "hbonds_desc": "OP1*OP2[2.77],OP2*N7[3.33]",
    hbonds = j_nonPairs["hbonds_desc"]
    for hbond in hbonds.split(","):
        if not "-" in hbond:
            # weird "hbonds" have * as separator instead of -
            continue
        at = [ hbond.split("(")[0].split("[")[0].split("-")[0], hbond.split("-")[1].split("(")[0].split("[")[0] ]
        for x in ind_in_chains:
            part = dictcode[at[x]]
            code = codenames[part]
            cc = "chain_"+c[x]
            rr = "res_"+resid[x]
            #intraRNA_hb[resid[a]-1].keys() = ["sug", "ph", "base"]
            #codenames = {1:"ph", 2:"sug", 3:"base"}
            d_intraRNA_hb[cc] = update_intraRNA(d_cc, rr, code, neighbors[pos[x]], 1)
    return d_intraRNA_hb

def get_nonPairs(js_nonPairs, d_intraRNA_hb, d_stacking):
    # stacking and nonpairing H-bonds:
    # TODO: check the chain of the stacking residue. Cf 1CVJ chainsP_res5 - chainM_res7
    for j in js_nonPairs:
        # non-paired interacting pairs of nt
        #if "hbonds_desc" in list(j.keys()):
        #   d_intraRNA_hb = get_nonPairs_hbonds(j, d_intraRNA_hb)
        if "stacking" in list(j.keys()):
            nts, resid, resint, c = parse_nts(j)
            ind_in_chains = [ x for x in range(2) if c[x] in chains]
            cc = ["chain_"+c[x] for x in ind_in_chains]
            rr = ["res_" + resid[x] for x in ind_in_chains]
            # position in the 5-nb vector of the nucleotide
            # pos = 2 : per default, it is not a neighbors
            pos = [2, 2]
            if c[0] in chains and c[1] == c[0]:
                if abs(resint[0]-resint[1]) <= 2:
                    # the nucl are neighbors in the same chain
                    pos = [resint[1] - resint[0] + 2, resint[0] - resint[1] + 2 ]
                    update_stacking(d["stacking"][cc[0]], rr[0], neighbors[pos[0]])
                    update_stacking(d["stacking"][cc[1]], rr[1], neighbors[pos[1]])
            else:
                for x in ind_in_chains:
                    update_stacking(d["stacking"][cc[x]], rr[x], "other")
    return d_intraRNA_hb, d_stacking

def get_pairs(js_pairs, d_ss, d_bptype):
    for j in js_pairs:
        nts, resid, resint, c = parse_nts(j)
        for i in [0, 1]:
            cc = "chain_"+c[i]
            rr = "res_"+resid[i]
            if c[i] in chains:
                if cc not in d_bptype:
                    d_bptype[cc] = {}
                if rr not in d_bptype[cc]:
                    d_bptype[cc][rr] = []
                d_bptype[cc][rr].append(j["name"])
                # per default, if a base makes a basepair,
                # the structure is "D" for undetermined double-stranded.
                # It can be detailed further in the code (stem, helix, junction ...)
                d_ss[cc][rr][0] = "D"
    return d_ss, d_bptype

def get_hairpins(js_hairpins, d_ss):
    for j in js_hairpins:
        nts = [ x.split(".") for x in j["nts_long"].split(",") ]
        for nr, n in enumerate(nts):
            c = nts[nr][1]
            resid = int(nts[nr][3])
            if nr==0 or nr==len(nts)-1:
                d_ss["chain_%s"%c]["res_%i"%resid] = ["D", 0, 0]
            else:
                d_ss["chain_%s"%c]["res_%i"%resid] = ["L", nr+1, j["num_nts"]-2]
    return d_ss

def get_junctions(js_junctions, d_ss):
    for j in js_junctions:
        for jj in j["bridges"]:
            nts = [ x.split(".") for x in jj["nts_long"].split(",") ]
            l = jj["num_nts"]
            for nr, n in enumerate(nts):
                if len(n) < 2 : continue
                resid, resint, c = n[3], int(n[3]), n[1]
                d_ss["chain_"+c]["res_%i"%resint] = ["J", nr , l]
    return d_ss

def get_ss(js_ssSegments, d_ss):
    #TODO : distinguish T from full-ss RNA
    for j in js_ssSegments:
        nts = [ x.split(".") for x in j["nts_long"].split(",") ]
        c = nts[0][1]
        if c not in chains: continue
        l = j["num_nts"]
        for nr, n in enumerate(nts):
            resid, resint, c = n[3], n[3], n[1]
            d_ss["chain_"+c]["res_%i"%resint] = ["S", nr + 1 , l]
            if int(resint)+1 in d["breaks"]["chain_"+c]:
                d_ss["chain_"+c]["res_%i"%resint] = ["T", l - nr , l]
        if int(nts[-1][3]) == lastres:
            #3'-term ssSegments
            resint = nts[-1][3]
            nr = len(nts)
            d_ss["chain_"+c]["res_%i"%resint] = ["T", l - nr , l]
        if int(nts[0][3]) == firstres:
            #5'-term ssSegments
            resint = nts[0][3]
            nr = 0
            d_ss["chain_"+c]["res_%i"%resint] = ["T", l - nr , l]

#def check_ss(d_ss):
#    for k in d_ss:

#dssrlist = open(sys.argv[1]).readlines()  #list of "3dna/xxxX-y.json"
chainsmodels = json.load(open(sys.argv[1])) # chainsmodels.json
outfile = sys.argv[2]                       # x3dna.json
out = testpathjson(outfile)
na = sys.argv[3]                            # "rna" or "dna"

if True: # for editing purpose
    ph = ["P","O1P", "O2P","OP1", "OP2", "O5'", "C5'"]
    sug = ["C4'", "C3'", "C2'", "C1'", "O2'", "O4'", "O3'"]
    base = ["N1","C2","N2","O2","N3","N4","O4","C4","C5","N6","C6","O6","N7","C7","C8","N9"]
    resnames = ["RG", "RC", "RA", "RU"]
    prefix = "R"
    if na == "dna":
        prefix = "D"
        resnames = ["DG", "DC", "DA", "DT"]
        sug = ["C4'", "C3'", "C2'", "C1'", "O4'", "O3'"]
    code = {1:ph, 2:sug, 3:base}
    codenames = {1:"ph", 2:"sug", 3:"base"}
    dbn = {
        "(" : "1",
        "." : "0",
        ")" : "1"
        }
    neighbors = ["n-2", "n-1", "other", "n+1", "n+2"]


dictcode = defaultdict(lambda: 0)
for a in range(1,4):
    for b in code[a]:
        dictcode[b] = a
#
    #intraRNA_hb  = Nb of RNA-NA H-bond for each nmap part (sug, ph, base)
    #RNAprot_hb  = Nb of RNA-NA/prot Hb for each nmap part
    #intraRNA_hb["res_"+resid[a]].keys() = ["sug", "ph", "base"]
    #loops: {L : hairpinloop, J : jonction, B : Bulge, I : intraloop}
    #
    # position in the 5-nb vector of the nucleotide:
    # {1: {"base": [0, 0, 1, 0, 0], "ph": [0, 0, 0, 0, 0], "sug": [0, 0, 0, 0, 0]},

# Nb of interactions with nt at [i-2, i-1, other, i+1, i+2]
vect = list(range(5))
for struct in sorted(chainsmodels.keys()):
    if struct in out.keys(): # structure already in the output
        continue
    pp(struct)
    d = chainsmodels[struct]
    out[struct] = d
    # get data from x3dna output
    inp = "3dna/%s-1-dssr.json"%struct #output from 3dna-dssr
    if not os.path.exists(inp):
        print("%s does not exist"%inp)
        continue
    js = json.load(open(inp))
    # initialise dictionary for the current structure
    chains = d["nachains"]  #list of IDs for rna/dna chains
    d, js = initialise_all(d, js)
    for c in chains:
        cc="chain_"+c
        firstres, lastres, dict_n3to1, dict_n2to3 = map_indices(struct, c)
        d["mapping"][cc] = dict_n3to1 # {"1": "20A"}
        for m in range(1, d["Nmodels"]+1): # for each model
            d = map_missing(d, cc, dict_n2to3)
            d["interface_protein"]["model_1"][cc] = map_interf(d, cc, dict_n3to1)
            d["breaks"][cc] = check_breaks(js, c)
            d = initialise_chain(d, c, m, dict_n3to1)
    d["RNAprot_hb"], d["intraRNA_hb"] = get_hbonds(js["hbonds"], d["RNAprot_hb"], d["intraRNA_hb"])
    d["intraRNA_hb"], d["stacking"] = get_nonPairs(js["nonPairs"], d["intraRNA_hb"], d["stacking"])
    d["ss"], d["bptype"] = get_pairs(js["pairs"], d["ss"], d["bptype"])
    d["ss"] = get_hairpins(js["hairpins"], d["ss"])
    d["ss"] = get_junctions(js["junctions"], d["ss"])

json.dump(out, open(outfile, "w"), indent = 2)
print("done", file=sys.stderr)

#((((((((((((((((((((((((((..((((....))))....))))))..(((..).)).......((((....)))).((((((...)))))).
