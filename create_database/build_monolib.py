#!/usr/bin/env python3
import sys, os
import numpy as np
sys.path.insert(0, os.environ["ATTRACTTOOLS"])
import rmsdlib

refedir = sys.argv[1]   # data/${na}lib
filelist = sys.argv[2]  # clean-iniparse.list
outpdir = sys.argv[3]   # data/${na}lib

assert os.path.exists(refedir)
assert os.path.exists(filelist), filelist
assert os.path.exists(outpdir)

def pp(*x):
    for i in x[:-1]:
        print(i, file=sys.stderr, end=' ')
    print(x[-1])

def _apply_matrix_multi(atoms_array, pivots, rotmats, offsets):
    cen = (atoms_array - pivots[:, None, :])
    newcoor = np.einsum("ijk,ikl->ijl", cen, rotmats) #diagonally broadcasted form of cen.dot(rotmats)
    newcoor += (pivots + offsets)[:, None, :]
    return newcoor

def fit_residue(coor, refnpy, fitatoms, refatoms):
    fitmask = np.array([(a in fitatoms) for a in refatoms])
    #fitcoor = coor[fitmask]np.compress(fitmask, coor, axis=1)
    #refe =  np.compress(fitmask, refnpy, axis=0)
    fitcoor = coor[:, fitmask]
    refe =  refnpy[0, fitmask]
    rotmats, offsets, _ = rmsdlib.multifit(fitcoor, refe)
    rotmats = rotmats.swapaxes(1,2)
    offsets = -offsets
    pivots = np.sum(fitcoor,axis=1) / float(fitcoor.shape[1])
    fitted_coor = _apply_matrix_multi(coor, pivots, rotmats, offsets)
    return fitted_coor

def get_coor(filelist):
    count = {}
    res_coor = {}
    for l in open(filelist):
        pdb = l.strip()
        print(pdb)
        res = []
        resi, Xresid = None, None
        for l in open(pdb):
            if l[30:38] == " XXXXXXX":
                #missing atom from --manual mode
                resi = None
                Xresid = l[21:27]
                continue
            resid = l[21:27]
            if resid == Xresid:
                #skip residue with missing atoms
                continue
            if resid != resi:
                # new residue
                if resi != None:
                    # process previous residue
                    nat = len(res)
                    res = np.array(res)
                    res_coor[resn][count[resn]] = res
                    count[resn] += 1
                resi = resid
                res = []
            resi = resid
            resn = l[19]
            if resn not in res_coor:
                refe = np.load("%s/%s-refe.npy"%(refedir, resn))
                res_coor[resn] =  np.zeros((500000, refe.shape[1], 3))
                count[resn] = 0
            x = float(l[30:38])
            y = float(l[38:46])
            z = float(l[46:54])
            res.append([x,y,z])
    for resn in res_coor:
        res_coor[resn] = res_coor[resn][:count[resn]]
    return res_coor

sugatoms = ["C1'", "C2'", "C3'", "C4'","C5'", "O2'", "O3'", "O4'"]
fitatoms = ["C1'", "C2'", "C3'", "C4'", "O4'"]

res_coor = get_coor(filelist)
ncoor = sum([len(res_coor[i]) for i in res_coor])

sugar = np.zeros((ncoor, len(sugatoms), 3))
nsug = 0

fitted = {}
for resn in res_coor:
    # parse reference
    refnpy = np.load("%s/%s-refe.npy"%(refedir, resn))
    refpdb = "%s/%s.pdb"%(refedir, resn)
    refatoms = [l[13:16].strip() for l in open(refpdb)]
    #fit nucleotide to reference
    fitted[resn] = fit_residue(res_coor[resn], refnpy, fitatoms, refatoms)
    np.save("%s/%s-fit.npy"%(outpdir, resn), fitted[resn])
    #extract sugar
    sugmask = np.array([(a in sugatoms) for a in refatoms])
    n = len(fitted[resn])
    sugar[nsug: nsug+n] = np.compress(sugmask, fitted[resn], axis=1)
    nsug += n

np.save("%s/sugar-fit.npy"%outpdir, sugar)
