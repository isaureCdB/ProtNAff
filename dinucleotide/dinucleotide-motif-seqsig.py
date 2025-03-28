"""
Find the sequence signature of a dinucleotide stacking motif
"""

import os
import sys

import numpy as np

_3_to_1 = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "MSE": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
}


def get_seqsig(motif, resids, refe_ca):
    resid1, resid2 = int(motif[8]), int(motif[10])
    s = []
    for resid in (resid1, resid2):
        ss = ["-"] * 15
        atoms = refe_ca[(resids >= resid - 7) & (resids <= resid + 7)]
        lastr = None
        for atom in atoms:
            if atom["resid"] != lastr:
                lastr = atom["resid"]
                pos = lastr - resid + 7
                # print(pos, lastr, resid)
                resname = atom["resname"].decode()
                ss[pos] = _3_to_1.get(resname, "X")
        s.append("".join(ss))
    return s[0] + "..." + s[1]


find_double_stacking_out = sys.argv[1]
find_double_stacking = [l.split() for l in open(find_double_stacking_out).readlines()]
find_double_stacking.sort()

curr_code = None
curr_chain = None
for motif in find_double_stacking:

    if motif[0] != curr_code:
        curr_code = motif[0]
        curr_chain = None
        pdbf = f"../database/pdbnpy/{motif[0]}.npy"
        if not os.path.exists(pdbf):
            print("MISS", curr_code)
            miss = True
            continue
        miss = False
        refe_pdb = np.load(pdbf)
        refe_pdb = refe_pdb[refe_pdb["model"] == 1]
    elif miss:
        continue
    chain = motif[6]
    if chain != curr_chain:
        curr_chain = chain
        refe_prot = refe_pdb[refe_pdb["chain"] == chain.encode()]
        refe_ca = refe_prot[refe_prot["name"] == b"CA"]
        resids = refe_ca["resid"]
    print(get_seqsig(motif, resids, refe_ca))
