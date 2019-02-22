#!/usr/bin/env python3
# Copyright 2017 - 2018 Isaure Chauvot de Beauchene (CNRS)

from Bio.PDB import *
import numpy as np
import sys
from scipy.spatial import cKDTree as KDTree
from scipy.sparse import dok_matrix
import re, os, json
from collections import defaultdict
import argparse

'''
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TODO:
    Distinguish between NA and protein chains to compute interface !!!
    Now it does not distinguish contacts with protein or with another NA chain
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'''

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

def pp(s):
    print(s, file=sys.stderr)

def testpathjson(jsonfile):
    if os.path.exists(jsonfile):
            #print >> sys.stderr, dictname +" = json.load(open('"+jsonfile+"'))"
        dictname = json.load(open(jsonfile))
    else:
        dictname = {}
        #exec( dictname + "= {}" )
    return dictname

def get_interf_res(atomlist, reftree, model):
    if len(atomlist) == 0 :
        return([], [], [])
    array = np.array([ a.get_coord() for a in atomlist])
    tree = KDTree(array)
    #list of list of interfacing nucl-atom per aa-atom
    interf = tree.query_ball_tree(reftree, args.cutoff)
    if max([len(i) for i in interf]) == 0:
        return([], [], [])
    # mask on aa-atoms: interface / non-interface
    mask = [len(j) > 0 for j in interf ]
    interf = array[np.array(mask)]
    # list of interfacing aa-atoms
    interf_at = set([id(a) for anr,a in enumerate(atomlist) if mask[anr]])
    interf_res = []
    for residue in model.get_residues():
        for atom in residue:
            if id(atom) in interf_at:
                interf_res.append(residue)
                break
    interftree = KDTree(interf)
    nainterf = natree.sparse_distance_matrix(interftree, args.cutoff)
    A = nainterf.toarray()
    A[A==0] = 99
    mininterf = np.min(A, axis=1)
    return(mininterf, interf, interf_res)

def interf_dist(mask, distances, indices):
    global threshold
    mindist = []
    for i in range(0, len(indices)-1, 2):
        selection = list(range(indices[i],indices[i+1]))
        try:
            d = min([ distances[i] for i in selection if mask[i] ])
        except:
            d = 999
        if d > threshold and d != 999: d = 99
        if d <= threshold : d = round(d,2)
        mindist.append(d)
    return(mindist)

def get_indexres_resid(model, nachains):
    resid, indexres = {}, {}
    i = 0
    for chain in model:
        c = chain.get_id()[0]
        indexres[c], resid[c] = [], []
        if c in nachains:
            for res in chain:
                #discarde non-nucleotide
                if len(res) <= 4: continue
                rr = str(res.get_resname())
                resname = rr.strip()
                if resname not in included:
                    continue
                indexres[c].append(i)
                for a in res:
                    if a.get_id() not in allatoms: continue
                    i += 1
                indexres[c].append(i)
                resid[c].append(str(res.get_id()[1]))
    return indexres, resid

def process_nachain(chain, struc, nuclpart, naatoms):
    global included, allatoms
    global part_code
    new_nuclpart = []
    new_naatoms = []
    c = chain.get_id()[0]
    NA_c = []
    for res in chain:
        #discarde non-nucleotide
        if len(res) <= 4: continue
        rr = str(res.get_resname())
        resname = rr.strip()
        if resname not in included:
            print("exclude: '%s in %s'"%(resname, struc))
            continue
        NA_c.append(res)
        for a in res:
            nat = a.get_id()
            if nat not in allatoms: continue
            naatoms.append(a)
            nuclpart.append(part_code[a.get_name()])
    return nuclpart, naatoms, NA_c

def process_protchain(chain, hetatoms, protatoms):
    global exclude
    for res in chain:
        #discarde ions, water and hetatoms
        if len(res) == 1 or res.get_resname() in exclude:
            continue
        for a in res:
            if a.get_id()[0] == "H": continue
            if res.get_resname()[0] == "MSE": continue
            if res.get_full_id()[3][0][:2] == "H_":
                hetatoms.append(a)
            else:
                protatoms.append(a)
    return hetatoms, protatoms

def get_chains(model, na):
    nachains, protchain = [], []
    for chain in model:
        c = chain.get_id()[0]
        countna = 0
        for r in chain:
            if r.has_id("O2'") and na == 'rna' :
                countna += 1
                if countna == 3:
                    nachains.append(c)
                    break
            if r.has_id("O3'") and na == 'dna' and not r.has_id("O2'") :
                countna += 1
                if countna == 3:
                    nachains.append(c)
                    break
            if r.has_id("CA"):
                protchain.append(c)
                break
    return nachains, protchain

def process_prot_interf(nuclpart, protmin, resid_c, indexres_c):
    assert 2 * len(resid_c) == len(indexres_c), (len(resid_c), len(indexres_c))
    d_interf = {}
    for k in ['ph', 'sug', 'base']:
        mask = [ j == k for j in nuclpart]
        mindist = interf_dist(mask, protmin, indexres_c)
        assert len(resid_c) == len(mindist), (len(resid_c), len(mindist))
        for ri, i in enumerate(mindist):
            if i != 99 :
                res = "res_%s"%resid_c[ri]
                if res not in d_interf:
                    d_interf[res] = {}
                d_interf[res][k] = i
    return d_interf

def process_het_interf(naatoms, hetmin, resid_c, indexres_c):
    assert 2 * len(resid_c) == len(indexres_c), (len(resid_c), len(indexres_c))
    d_hetatm = {}
    maskrna = [ True for j in naatoms ]
    mindist_het = interf_dist(maskrna, hetmin, indexres_c)
    for ri, i in enumerate(mindist_het):
        if i != 99 :
            res = resid_c[ri]
            d_hetatm["res_%s"%res] = i
    return d_hetatm

parser = argparse.ArgumentParser(description=__doc__,
                          formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("cutoff",help="5A cutoff for interface contact", type=float)
parser.add_argument("inpdir",help="brutPDBs")
parser.add_argument("outdir",help="chainsmodels")
parser.add_argument("inlist",help="pdbcodes.txt")
parser.add_argument("mutatelist",help="$ATTRACTTOOLS/..//allatoms/$na-mutate.list")
parser.add_argument("outjson",help="chainsmodels.json")
parser.add_argument("na",help="dna or na")
parser.add_argument("--delete", help="remove from outpjson entries not in inplist", action="store_true")
parser.add_argument("--replace", help="replace entries already in outp", action="store_true")

args = parser.parse_args()

parser = PDBParser()
pdbl = PDBList()
io = PDBIO()
threshold = 8

inlist = [l.strip() for l in open(args.inlist,'r').readlines()]
modified = set([l.split()[0] for l in open(args.mutatelist).readlines()])

####################################
# initialise lists and dictionaries
####################################
if True:
    exclude = set([ "FHU", "ATP", "ADP", "AMP", "GTP", "GDP", "GMP", "CTP", \
                "CDP", "CMP", "UTP", "UDP", "UMP", "SAM","O2I", "2BA", \
                "574", "6MZ", "A2P", "APZ", "C2E", "CG1", "DU", "GH3", \
                "MGT", "URU", "SUC", "02I", "HOH", "SO4","CYS", "12A"])
    if args.na == 'dna':
        canonical = set(["DT", "DC", "DG", "DA", "C", "T", "A", "G", "PSU", "7AT"])
        sug = ["C4'", "C3'", "C2'", "C1'", "O4'"]
    elif args.na == 'rna':
        canonical = set(["RU", "RC", "RG", "RA", "C", "U", "A", "G", "PSU", "7AT"])
        sug = ["C4'", "C3'", "C2'", "C1'", "O2'", "O4'"]
    included = canonical | modified
    ph = ["P","O1P", "O2P","OP1", "OP2", "O5'", "O3'", "C5'"]
    base = ["N3", "C2", "N2", "N1", "C4", "O6", "C5", "N7", "C8", "N6", "C6","N9","N4", "O2", "O4"]
    code = { 'ph':ph, 'sug':sug, 'base':base}
    allatoms = ph + sug + base

    part_code = defaultdict(str)
    for a in ['ph', 'sug', 'base']:
        for b in code[a]:
            part_code[b] += a

#################################################################

out_init = testpathjson(args.outjson)  # chainsmodels.json
if args.delete:
    pp("delete entries not in inplist")
    out = {}
else:
    out = out_init

pdbcodes = set(out_init.keys())        # already processed PDBs

for struc in inlist:
    if struc in out_init and not args.replace and args.delete:
        out[struc] = out_init[struc]
        continue
    print(struc, file=sys.stderr)
    if not os.path.exists(args.inpdir+"/"+struc+'.pdb'):
        print("ERROR: "+args.inpdir+"/"+struc+".pdb does not exist")
        continue
    try:
        structure = parser.get_structure(struc, args.inpdir+"/"+struc+'.pdb')
        out[struc] = {}
        d = out[struc]
        d['method'] = structure.header['structure_method']
        d['resolution'] = structure.header['resolution']
        d['interface_protein'] = {}
        d['interface_hetatoms'] = {}
        d['hetnames'] = {}
        nachains, protchain = get_chains(structure[0], args.na)
        Nmodels = len(structure)
        indexres, resid = get_indexres_resid(structure[0], nachains)
        m = 0
        for model in structure:
            m+=1
            mm = "model_%i"%m
            NA = {}
            naatoms, protatoms, hetatoms, nuclpart = [], [], [], []
            for chain in model:
                c = chain.get_id()[0]
                if c in nachains:
                    nuclpart, naatoms, NA_c = process_nachain(chain, struc, nuclpart, naatoms)
                    if NA_c is None: continue
                    NA[c] = NA_c
                else:
                    hetatoms, protatoms = process_protchain(chain, hetatoms, protatoms)
            naarray = np.array([a.get_coord() for a in naatoms])
            natree = KDTree(naarray)
            protmin, _ , prot_interf = get_interf_res(protatoms, natree, model)
            hetmin, _ , het_interf = get_interf_res(hetatoms, natree, model)
            io.set_structure(model)
            selection = MySelect(model, prot_interf + het_interf)
            io.save("interface/%s-%i.pdb"%(struc, m), selection)
            io.set_structure(model)
            d['interface_protein'][mm] = {}
            for c in nachains:
                cc = "chain_%s"%c
                # Interface with protein
                d_prot_interf = process_prot_interf(nuclpart, protmin, resid[c], indexres[c])
                if len(d_prot_interf.keys()):
                    d['interface_protein'][mm][cc] = d_prot_interf
                io.save("%s/%s%s-%i.pdb"%(args.outdir, struc, c, m), MySelect(model,NA[c]))
                # Interface with HETATM
                if len(het_interf) > 0:
                    d_het_interf = process_het_interf(naatoms, hetmin, resid[c], indexres[c])
                    if len(d_het_interf.keys()):
                        if str(m) not in d['interface_hetatoms']:
                            d['interface_hetatoms'][mm] = {}
                        d['interface_hetatoms'][mm][cc] = d_het_interf
                    if model is structure[0]:
                        hetnames = [ r.get_resname() for r in het_interf ]
                        d['hetnames'][mm] = hetnames
        d['nachains'] = nachains
        d['protchain'] = protchain
        d['Nmodels'] = m
    except:
        print('ERROR')
        pass
        raise

json.dump(out, open(args.outjson, "w"), indent = 2, sort_keys = True)
