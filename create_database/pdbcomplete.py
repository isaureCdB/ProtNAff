#!/usr/bin/env python3
# Copyright Sjoerd J. De Vries (TUM), Isaure Chauvot de Beauchene (CNRS)

from __future__ import print_function
import sys, os, tempfile, json
import scipy.spatial
import numpy as np
from copy import deepcopy
import topology, traceback

def pp(i):
    print(i, file=sys.stderr)

currdir = os.path.abspath(os.path.split(__file__)[0])

class nalib(object):pass

class FixError(Exception):pass

# pdb2pqr returns another atom name for that RNA atom than the opls ff.
# check_pdb will chock on it.
map_atnames = {"HO2'": "H2''"}
na_resnames = {"RA","RC","RG","RU","DA","DC","DG","DT","A","T","C","G","U"}
x5pho = {"C5'","H5'", "H5''","P","O1P","O2P","O5'","O5T","H5T"} #############

def run_pdb2pqr(pdblines):
    """
    Runs PDB2PQR on the input
    """
    oldsyspath = sys.path
    pqrhandle, pqrfile = tempfile.mkstemp()
    tmphandle, tmpfile = tempfile.mkstemp()
    tmpf = open(tmpfile, "w")
    err = "/tmp/pdb2pqr.err"
    for l in pdblines:
        if l.find("XXX") == -1:
            print(l, file=tmpf)
    tmpf.close()

    try:
        import pdb2pqr
        args = [pdb2pqr.__file__, "--ff=AMBER", tmpfile, pqrfile]
        #if pdb2pqr.PACKAGE_PATH != "":
        #  sys.path.extend(pdb2pqr.PACKAGE_PATH.split(":"))
        oldstdout = sys.stdout
        sys.stdout = sys.stderr
        oldargv = list(sys.argv)
        sys.argv[:] = args
        try:
            pdb2pqr.mainCommand(args)
        except ValueError as exc:
            msg = exc.args[0]
            if msg.startswith("This PDB file is missing too many"):
                raise FixError(msg)
            else:
                raise
        finally:
            sys.argv[:] = oldargv
            sys.stdout = oldstdout
        pqr = os.fdopen(pqrhandle)
        pqrlines = pqr.readlines()
        pqr.close()
    except:
        os.system("pdb2pqr --ff=AMBER %s %s 2> %s"%(tmpfile, pqrfile, err))
        os.system("cat %s"%err)
        pqrlines = list(open(pqrfile))
    finally:
        os.remove(tmpfile)
        os.remove(pqrfile)

    result = []
    repl = (
      (" H  ", " HN "),
      (" H  ", " HT1"),
      (" H2 ", " HT2"),
      (" H3 ", " HT3"),
    )
    for l in pqrlines:
        result.append(l)
        atom2 = l[12:16]
        for pin, pout in repl:
            if atom2 == pin:
                p = l[:12] + pout + l[16:]
                result.append(p)
    return result

def update_patches(pdb, top_patches):
    """
    Updates PDB His/Cys patches based on the presence of hydrogens
    """
    for res in pdb:
        if res.resname == "HIS":
            protons = []
            for proton in ("HE1", "HE2", "HD1", "HD2"):
                if proton in res.coords:
                    protons.append(proton)

            if protons == ["HE2", "HD1"]: #HIS+ from aareduce
                protons = ["HE1", "HE2", "HD1", "HD2"]
            elif protons == ["HE2"]: #HISE from aareduce
                protons = ["HE1", "HE2", "HD2"]
            elif protons == ["HD1"]: #HISD from aareduce
                protons = ["HE1", "HD1", "HD2"]

            if len(protons) < 3:
                continue
            if "HE2" not in protons or "HE1" not in protons:
                res.topology.patch(top_patches["HISD"])
            if "HD1" not in protons or "HD2" not in protons:
                res.topology.patch(top_patches["HISE"])
        if res.resname == "CYS":
            if "HG" not in res.coords:
                res.topology.patch(top_patches["DISU"])

def pdbcomplete(pdb, other):
    """
    Completes a PDB by taking missing coordinates from the other PDB
    """
    assert len(pdb) == len(other), (len(pdb), len(other))
    for r1, r2 in zip(pdb, other):
        for atom in r2.coords:
            if atom not in r1.coords:
                r1.coords[atom] = r2.coords[atom]

def pdbfix(pdb, refe):
    """
    Fixes His and Cys hydrogens in a PDB that are present in the reference
     Hydrogens are filled in by computing the N-H or S-H vector from the reference
     and applying it to the N or S in the reference
    Also copies the OXT position from the reference if it is missing
    """
    assert len(pdb) == len(refe)
    for res, ref in zip(pdb, refe):
        if "oxt" in res.topology.atomorder:
            if "OXT" not in res.coords and "OXT" in ref.coords:
                res.coords["OXT"] = ref.coords["OXT"]
        if res.resname == "HIS":
            for proton in ("HD1", "HE2"):
                if proton in res.coords: continue
                if proton not in ref.coords: continue
                nitro = "N" + proton[1:]
                if nitro not in res.coords or nitro not in ref.coords: #can't fix...
                    continue
                nref = ref.coords[nitro]
                href = ref.coords[proton]
                vref = href[0]-nref[0],href[1]-nref[1],href[2]-nref[2]
                nres = res.coords[nitro]
                hres = nres[0] + vref[0], nres[1] + vref[1], nres[2] + vref[2]
                res.coords[proton] = hres
        if res.resname == "CYS":
            proton = "HG"
            sulfur = "SG"
            if proton in res.coords: continue
            if proton not in ref.coords: continue
            if sulfur not in res.coords or sulfur not in ref.coords: #can't fix...
                continue
            sref = ref.coords[sulfur]
            href = ref.coords[proton]
            vref = href[0]-sref[0],href[1]-sref[1],href[2]-sref[2]
            sres = res.coords[sulfur]
            hres = sres[0] + vref[0], sres[1] + vref[1], sres[2] + vref[2]
            res.coords[proton] = hres

def load_nalib(libname):
    lib = nalib()
    lib.dir = currdir + "/" + libname
    lib.sugar = {}
    lib.base = {}
    lib.nucl = {}
    lib.ph = {}
    lib.pho5 = {}
    ph = {
        # !!! keep the order of atoms in the library!
        "all_atoms":     ["P", "O1P", "O2P", "O5'", "C5'", "C4'", "N3", "C5", "C3'"], # for masking (corresponds to "coor")
        "atoms":     ["P", "O1P", "O2P", "O5'"], # for completion
        "fit_atoms": ["P", "O1P", "O2P", "O5'", "C5'", "C4'", "C3'"], #atoms to fit on
        "rmsd_atoms": ["P", "O1P", "O2P", "O5'", "C5'", "C4'", "N3", "C5", "C3'"], # to compute best RMSD (avoid base-phosphate clashes),
        "max_missing": 4,
        }
    pho5 = {
    "coor": np.load(lib.dir + "/5PHO.npy"),
    "atoms": [l[12:16].strip() for l in open(lib.dir + "/5PHO.pdb") if l.startswith("ATOM")],
    "max_missing": 5,
    }
    pho5["fit_atoms"] = pho5["atoms"]
    pho5["rmsd_atoms"] = pho5["atoms"]
    sugar = {
        "coor": np.load(lib.dir + "/sugar.npy"),
        "atoms": [l[12:16].strip() for l in open(lib.dir + "/sugar.pdb") if l.startswith("ATOM")],
        "max_missing": 2,
        }
    sugar["all_atoms"] = sugar["atoms"]
    sugar["fit_atoms"] = sugar["atoms"]
    sugar["rmsd_atoms"] = sugar["atoms"]
    bases = ["A","C","G","U"]
    if libname == "dnalib":
        bases = ["A","C","G","T"]
    for nuc in bases:
        lib.sugar[nuc] = sugar
        lib.nucl[nuc] = {
            "coor": np.load(lib.dir + "/%s.npy" % nuc),
            "missing": 23,
        }
        nuclatoms = [l[12:16].strip() for l in open(lib.dir + "/%s.pdb" % nuc) if l.startswith("ATOM")]
        lib.nucl[nuc]["all_atoms"] = nuclatoms
        lib.nucl[nuc]["atoms"] = nuclatoms
        lib.nucl[nuc]["fit_atoms"] = lib.sugar[nuc]["atoms"]
        lib.nucl[nuc]["rmsd_atoms"] = nuclatoms

        lib.base[nuc] = {}
        #lib.base[nuc]["coor"] = np.array([(float(l[30:38]),float(l[38:46]),float(l[46:54])) for l in open(base) if l.startswith("ATOM")])
        lib.base[nuc]["coor"] = np.load(lib.dir + "/base%s.npy" % nuc)[None]
        baseatoms = [l[12:16].strip() for l in open(lib.dir + "/base%s.pdb" % nuc) if l.startswith("ATOM")]
        lib.base[nuc]["atoms"] = baseatoms
        lib.base[nuc]["fit_atoms"] = baseatoms
        lib.base[nuc]["rmsd_atoms"] = baseatoms
        lib.base[nuc]["max_missing"] = len(baseatoms) - 3

        lib.ph[nuc] = ph.copy()
        lib.ph[nuc]["coor"] = np.load(lib.dir + "/%s.npy" % nuc)
        phlist = set(lib.ph[nuc]["all_atoms"])
        ph_indices = [ anr for anr,a in enumerate(lib.nucl[nuc]["all_atoms"]) if a in phlist]
        lib.ph[nuc]["coor"] = lib.nucl[nuc]["coor"][:, ph_indices]

        lib.pho5[nuc] = pho5
    return lib


def check_missing(pdb):
    for nr, res in enumerate(pdb):
        top = res.topology
        for a in top.atomorder:
            aa = a.upper()
            if aa[0].startswith("H"): continue
            if aa not in res.coords:
                return 1
    return 0

def process_pdb(pdb):
    other_coor, atom_to_residue = [], []
    for nr, res in enumerate(pdb):
        top = res.topology
        for a in top.atomorder:
            aa = a.upper()
            if aa[0].startswith("H"): continue
            if aa not in res.coords:
                continue
            xyz = res.coords[aa]
            other_coor.append(xyz)
            atom_to_residue.append(nr)
    return np.array(other_coor), np.array(atom_to_residue)

def _apply_matrix(atoms_array, pivots, rotmats, offsets):
    cen = (atoms_array - pivots[:, None, :])
    newcoor = np.einsum("ijk,ikl->ijl", cen, rotmats) #diagonally broadcasted form of cen.dot(rotmats)
    newcoor += (pivots + offsets)[:, None, :]
    return newcoor

def check_clashes(nl, libconf, lib_complete_indices, new_at, tree, atom_to_residue):
    #for nl, libconf in enumerate(libcoor_fitted_sorted):
    lib_complete = libconf[lib_complete_indices]
    # get clashes with any original residue
    neighbors = tree.query_ball_point(lib_complete, r=2)
    neighbors = np.concatenate(neighbors).astype(int)
    # ! the completed residue clashes w. himself
    clash_res = np.unique(atom_to_residue[neighbors])
    # get clashes with other added atom
    new_clash = 0
    if len(new_at):
        new_atoms = np.array(new_at)
        dd = new_atoms[None,:,:] - lib_complete[:,None,:]
        dist_new = np.sum(dd*dd, axis=2)
        new_clash = np.sum(dist_new < 4)
        #m = np.min(dist_new)
        #print('nb new clash: %i, min %f'%(new_clash, m))
    if len(clash_res) < 2 and not new_clash:
        return 0
    return 1

def rank_library(res, libcoor, fit_atoms, rmsd_atoms, all_atoms):
    '''
    rank nucleotides in library by best-fitting to our target residue
    '''
    import rmsdlib
    from scipy.spatial.distance import cdist
    #select atoms in the order of the library atoms
    coor_atoms = [a for a in all_atoms if a in res.coords and a in fit_atoms]
    coor = np.array([res.coords[a] for a in all_atoms if a in res.coords and a in fit_atoms])
    fit_mask = np.array([(a in fit_atoms and a in res.coords) for a in all_atoms])
    #fit the mononucl lib on the nucl to repair:
    ## get the rotation-translation to apply to the rank_library
    libcoor_fit = libcoor[:,fit_mask]
    rotmats, offsets, rmsds = rmsdlib.multifit(libcoor_fit,coor)
    rotmats = rotmats.swapaxes(1,2)
    offsets = -offsets
    x = rmsds.argsort()
    pivots = libcoor_fit.sum(axis=1)/libcoor_fit.shape[1]
    ## fit on fit_atoms and evaluate the rmsd on the rmsd_atoms
    common_atoms = set([a for a in rmsd_atoms if a in res.coords])
    rmsd_mask = np.array([(a in common_atoms) for a in rmsd_atoms])
    libcoor_rmsd_unfitted = libcoor[:,rmsd_mask]
    libcoor_rmsd_fitted = _apply_matrix(libcoor_rmsd_unfitted, pivots, rotmats, offsets)
    coor_rmsd = np.array([res.coords[a] for a in rmsd_atoms if a in common_atoms])
    dist = libcoor_rmsd_fitted - coor_rmsd
    d = np.einsum("ijk,ijk->i", dist, dist)
    ### d = rmsds ### use the fit rmsds as an approx. for debugging
    # sort library indices by RMSD
    lib_indices = np.argsort(d)
    libcoor_fitted = _apply_matrix(libcoor, pivots, rotmats, offsets)
    libcoor_fitted_sorted = libcoor_fitted[lib_indices]
    return libcoor_fitted_sorted


def apply_nalib(pdb, lib, manual, heavy=True):
    """
    Adds missing atoms using a nucleotides lib
    """
    print("apply_nalib", file=sys.stderr)
    #assert heavy == True #TODO: remove this when the rna library has hydrogens
    syspath = list(sys.path)
    if "ATTRACTTOOLS" in os.environ:
        sys.path.insert(0, os.environ["ATTRACTTOOLS"])
    import rmsdlib
    from scipy.spatial import cKDTree
    sys.path[:] = syspath
    any_missing = False
    #check if any missing atoms in the whole PDB
    #If any atom missing, then process the PDB to get all other atoms
    any_missing = check_missing(pdb)
    if any_missing:
        #coordinates toward which one should test clashes
        other_coor, atom_to_residue = process_pdb(pdb)
        tree  = cKDTree(other_coor)
    new_at_all = []
    for nr, res in enumerate(pdb):
        if res.resname not in na_resnames: continue
        new_at = []
        try:
            while 1: #keep fixing as long as we can
                missing = set()
                top = res.topology
                for a in top.atomorder:
                    aa = a.upper()
                    #if aa[0].startswith("H"): continue ###
                    if aa.startswith("H") or aa == '5pho': continue ###
                    if aa not in res.coords:
                        missing.add(aa)
                if not missing: break
                nuc = res.resname[1]
                fixmode = None
                for fixmode in ("ph", "sugar", "base", "nucl"): #from high to low priority
                    #we can fix if there are any missing atoms, and there are at least three non-lacking atoms
                    #if fixmode == "sugar": continue
                    sublib = getattr(lib, fixmode) # lib.ph or lib.sugar or ...
                    atoms = sublib[nuc]["atoms"]
                    all_atoms = sublib[nuc]["all_atoms"]
                    fit_atoms = sublib[nuc]["fit_atoms"]
                    rmsd_atoms = sublib[nuc]["rmsd_atoms"]
                    libcoor = sublib[nuc]["coor"]
                    assert libcoor.shape[-2] == len(all_atoms)
                    if any([(m in atoms) for m in missing]) and \
                     len([a for a in fit_atoms if a in res.coords]) >= 3:
                        break
                else:
                    msg = 'residue %s could not be fixed'%res.resid
                    print("missing atoms: ", file=sys.stderr)
                    print(missing, file=sys.stderr)
                    break
                    #raise FixError(msg)
                print("fixing %s"%fixmode, file=sys.stderr)
                libcoor_fitted_sorted = rank_library(res, libcoor, fit_atoms, rmsd_atoms, all_atoms)
                lib_complete_indices = []
                for anr, a in enumerate(all_atoms):
                    if a in missing or fixmode == "base":
                        lib_complete_indices.append(anr)
                #TODO: change clashing threshold when not --heavy
                #optimize: if clashes, take next nucleotide in mononucl_library
                print('optimize resid %s %s'%(res.resname, res.resid), file=sys.stderr)
                for nl, libconf in enumerate(libcoor_fitted_sorted):
                    clashing = check_clashes(nl, libconf, lib_complete_indices, new_at_all, tree, atom_to_residue)
                    if not clashing:
                        break
                if nl == len(libconf):
                    raise FixError("all lib conformer clash on resid %s %s"%(res.resname, res.resid))
                #
                for anr, a in enumerate(atoms):
                    if anr in lib_complete_indices: # add "or 1" to replace all atoms, for debugging
                        x,y,z = libconf[anr]
                        res.coords[a] = x,y,z
                        new_at.append([x,y,z])
                if fixmode == "nucleotide" : break
        except FixError as err:
            if manual:
                e = "\n" + "!"*60 + "\n"
                print(e + "WARNING: " + err.args[0] + e, file=sys.stderr)
            else:
                raise
        for a in new_at:
            new_at_all.append(a)


def pdb_lastresort(pdb, heavy=True):
    """
    Last-resort fixes to prevent errors
    """
    for res in pdb:
        if "oxt" in res.topology.atomorder:
            #Put an OXT at the same place as the O
            if "OXT" not in res.coords and "O" in res.coords:
                res.coords["OXT"] = res.coords["O"]
        if not heavy:
            if "H5T" not in res.coords and "O5'" in res.coords:
                res.coords["H5T"] = res.coords["O5'"]
            if "H3T" not in res.coords and "O3'" in res.coords:
                res.coords["H3T"] = res.coords["O3'"]

def get_atomcode_resname(l, mutations, mapnuc, is_dna, is_rna):
    atomcode = l[12:16].strip()
    resname = l[17:20].strip()
    if resname=="HYP" and atomcode=='OD1':atomcode='OG1'
    if resname=="HYP" and atomcode=='HD1':atomcode='HG1'
    if resname in ["HIE", "HIP", "HSD", "HSE", "HYP"]: resname="HIS"
    if resname=="MSE" and atomcode=='SE':atomcode='SD'
    if resname=="MSE":  resname="MET"
    if resname in mutations:
        resname = mutations[resname]
    if resname in mapnuc:
        if is_dna:
            resname = mapnuc[resname][0]
        elif is_rna:
            resname = mapnuc[resname][1]
        elif not heavy:
            raise ValueError("PDB contains a nucleic acid named \"%s\", but it could be either RNA or DNA. Please specify the --dna or --rna option" % resname)
        if resname is None:
            if is_dna: na = "DNA"
            if is_rna: na = "RNA"
            raise ValueError("'%s' can't be %s" % (l[17:20].strip(), na))
    return atomcode, resname

#Mapping of nucleic-acid codes to DNA/RNA
mapnuc = {
  "A": ["DA", "RA"],
  "A3": ["DA", "RA"],
  "A5": ["DA", "RA"],
  "DA3": ["DA", None],
  "DA5": ["DA", None],
  "ADE": ["DA", "RA"],
  "C": ["DC", "RC"],
  "C3": ["DC", "RC"],
  "C5": ["DC", "RC"],
  "DC3": ["DC", None],
  "DC5": ["DC", None],
  "CYT": ["DC", "RC"],
  "G": ["DG", "RG"],
  "G3": ["DG", "RG"],
  "G5": ["DG", "RG"],
  "DG3": ["DG", None],
  "DG5": ["DG", None],
  "GUA": ["DG", "RG"],
  "T": ["DT", None],
  "T3": ["DT", None],
  "T5": ["DT", None],
  "DT3": ["DT", None],
  "DT5": ["DT", None],
  "THY": ["DT", None],
  "U": [None, "RU"],
  "U3": [None, "RU"],
  "U5": [None, "RU"],
  "URA": [None, "RU"],
  "URI": [None, "RU"],
  }
mapnucrev = {
  "DA":"A",
  "RA":"A",
  "DC":"C",
  "RC":"C",
  "DG":"G",
  "RG":"G",
  "DT":"T",
  "RU":"U",
  }
#nucnames = set(mapnucrev.keys()) and set(mapnuc.keys())

class PDBres:
    def __init__(self, chain, resid, resname, topology):
        self.chain = chain
        self.resid = resid
        self.resname = resname
        self.coords = {}
        self.chainfirst = False
        self.chainlast = False
        self.nter = False
        self.cter = False
        self.topology = topology

#TODO: remove non-fixable residues
#TODO: add 5'-phosphate
def read_pdb(pdblines, pdbname, top_residues, is_dna, is_rna, heavy,
            mutations={}, add_termini=False, modbase=False, modres=False):
    repl = (
      ("H","HN"),
      ("HT1","HN"),
      ("OP1","O1P"),
      ("OP2","O2P"),
      ("H1","HN"),
      ("OP3","O5T"),
      ("HO5'","H5T"),
      ("HO3'","H3T"),
      )
    atomlines = []
    ter_indices = []
    pdblines = list(pdblines)
    if (modbase or modres):
        # check for each residue containing some HETATM if it is a base/aa
        # !!! HETATM with any missing C/CA/N or C[1-5]' is discarded !!!
        # => don't do it for residues that have only "ATOM" lines
        res0 = {}
        for l in pdblines:
            if l.startswith("HETATM"):
                resid = l[22:27]
                if resid not in res0:
                    res0[resid] = set()
                atomcode = l[12:16].strip()
                res0[resid].add(atomcode)
        res_ok = set()
        res_atoms = set(["CA", "C", "N"])
        base_atom = set(["C%i'"%i for i in range(1, 6)])
        for r in res0:
            ratoms = res0[r]
            if modres and res_atoms.issubset(ratoms):
                res_ok.add(r)
            elif modbase and base_atom.issubset(ratoms):
                res_ok.add(r)
        for l in pdblines:
            if l.startswith("ATOM") or l.startswith("TER"):
                atomlines.append(l)
            elif l.startswith("HETATM"):
                resid = l[22:27]
                if resid in res_ok:
                    atomlines.append(l)
        pdblines = atomlines
    else:
        pdblines = [l for l in pdblines if l.startswith("TER") or l.startswith("ATOM")]

    atlines = [l for l in pdblines if l.split()[0]!="TER"]
    if len(atlines) == 0:
        raise ValueError("PDB '%s' contains no atoms" % pdbname )
    #
    curr_res = None
    ter_indices = []
    pdbres = []
    for l in pdblines:
        if l.startswith("TER"):
            if curr_res is not None:
                ter_indices.append(curr_res)
            continue
        if l[16] not in (" ", "A"):
            #only keep the first of alternative conformations
            traceback.print_exc("PDB contains alternative positions, those will be deleted in output")
            continue
        if l[30:38] == " XXXXXXX": continue #missing atom from --manual mode
        atomcode, resname = get_atomcode_resname(l, mutations, mapnuc, is_dna, is_rna)
        chain = l[21]
        resid = l[22:27]
        x, y, z = [float(x) for x in [l[30:38], l[38:46], l[46:54]]]
        newres = False
        nter = False
        chainfirst = False
        if curr_res is None:
            newres = True
            chainfirst = True
            if add_termini: nter = True
        elif chain != curr_res.chain:
            newres = True
            chainfirst = True
            curr_res.chainlast = True
            if add_termini:
                nter = True
                curr_res.cter = True
        elif resid != curr_res.resid or resname != curr_res.resname:
            newres = True
        if newres:
            try:
                if resname is None: raise KeyError
                topr = deepcopy(top_residues[resname])
            except KeyError:
                raise KeyError("Residue type %s not known by the topology file" % resname)
            curr_res = PDBres(chain, resid, resname, topr)
            if chainfirst: curr_res.chainfirst = True
            if nter: curr_res.nter = True
            pdbres.append(curr_res)
        curr_res.chain = chain
        curr_res.resid = resid
        curr_res.coords[atomcode] = (x,y,z)
        for pin, pout in repl:
            if atomcode != pin: continue
            curr_res.coords[pout] = (x,y,z)
    if curr_res is not None:
        curr_res.chainlast = True
        if add_termini:
            curr_res.cter = True
    return pdbres, ter_indices

def eval_moltype(pdbres):
    is_na = False
    is_prot = False
    for res in pdbres:
        if res.resname in na_resnames:
            is_na = True
        if "CA" in res.coords.keys():
            is_prot = True
        if is_prot and is_na:
            break
    if is_prot:
        print('input contains protein', file=sys.stderr)
    if is_na:
        print('input contains nucleic acids', file=sys.stderr)
    return is_prot, is_na

def write_pdb(pdbres, args, patches={}, heavy = False,
            one_letter_na = False, write_missing=True, pdb2pqr=False):
    pdblines = []
    mapping = []
    mutations_5ter = {}
    firstchain = True
    atomcounter = 1
    rescounter = 1
    map_atnames = {"HO2'": "H2''",
                    "H5" : "H5''"}
    for res in pdbres:
        top = res.topology
        if res.chainfirst and not pdb2pqr and not firstchain:
            pdblines.append("TER")
        firstchain = False
        for a in top.atomorder:
            atom = top.atoms[a]
            if a.lower().startswith("h"):
                if heavy: continue
                if atom["charge"] == 0: continue
            aa = a.upper()
            at = aa.strip()
            if at in map_atnames:
                at = map_atnames[at]
            x = " XXXXXXX"
            y = x; z = x
            if at in res.coords:
                x,y,z = ("%8.3f" % v for v in res.coords[at])
            elif not write_missing:
                continue
            xyz = x + y + z
            a0 = aa
            if len(a0) < 4:
                a0 = " " + a0 + "   "[len(a0):]
            resname = res.resname
            atomname = a0
            if one_letter_na and resname in mapnucrev:
                resname = mapnucrev[resname]
            chain = res.chain
            resid = res.resid
            if args.renum_res is not None:
                resid = rescounter
            if pdb2pqr:
                if res.chainfirst and res.resid in patches:
                    print(patches[res.resid]) #####################
                    if None in patches[res.resid] or "5pho" in patches[res.resid]:
                        if aa in x5pho:
                            mutations_5ter['PHO'] = resname
                            resname = "PHO"
            pdblines.append("ATOM%7d %4s%4s %s%4s   %s%7.3f 1.00" % \
              (atomcounter, atomname, resname, res.chain, resid, xyz, atom["charge"]))
            atomcounter += 1
        mapping.append((res.resid, rescounter))
        rescounter += 1
    return pdblines, mapping, mutations_5ter

def termini_pdb(pdbres, nter, cter):
    xter = nter, cter
    for n in range(2):
        ter = xter[n]
        for resnr in ter:
            r = [res for res in pdbres if res.resid == resnr]
            if len(r) == 0:
                raise ValueError("Cannot find residue %d" % resnr)
            elif len(r) > 1:
                raise ValueError("Multiple residues %d" % resnr)
            res = r[0]
            if n == 0: res.nter = True
            else: res.cter = True

def patch_pdb(pdbres, patches, top_patches):
    pp(patches)
    for res in pdbres:
        if res.resid in patches:
            for p in patches[res.resid]:
                if p is None: continue
                pp(res.resid)
                pp(p)
                res.topology.patch(top_patches[p.upper()])
        elif len(pdbres) > 1 and "ca" in res.topology.atomorder: #protein
            if res.nter:
                if res.resname == "PRO":
                    res.topology.patch(top_patches["PROP"])
                else:
                    res.topology.patch(top_patches["NTER"])
            if res.cter:
                res.topology.patch(top_patches["CTER2"])
        elif len(pdbres) > 1 and "p" in res.topology.atomorder: #DNA/RNA
            if res.chainfirst:
                if res.nter:
                    res.topology.patch(top_patches["5PHO"])
                else:
                    res.topology.patch(top_patches["5TER"])
            if res.chainlast:
                res.topology.patch(top_patches["3TER"])

def check_pdb(pdbres, heavy=False):
    map_atnames = {"HO2'": "H2''"}
    for res in pdbres:
        top = res.topology
        for a in top.atomorder:
            atom = top.atoms[a]
            if a.lower().startswith("h"):
                if heavy: continue
                if atom["charge"] == 0: continue
            aa = a.upper()
            at = aa.strip()
            if at in map_atnames:
                at = map_atnames[at]
            if at not in res.coords:
                raise ValueError('Missing coordinates for atom "%s" in residue %s %s%s' % (aa.strip(), res.resname, res.chain, res.resid))

def set_reference(pdbres, pdbreferes):
    if len(pdbres) != len(pdbreferes):
        raise ValueError("PDB and reference do not have the same number of residues, %d vs %s" % (len(pdbres), len(pdbreferes)))
    for n in range(len(pdbres)):
        pdbr, refr = pdbres[n], pdbreferes[n]
        if pdbr.resname != refr.resname:
            rsid = pdbr.resid
            if refr.resid != pdbr.resid: rsid = "%s(%s)" % (pdbr.resid, refr.resid)
            raise ValueError("PDB and reference are different at resid %s: %s vs %s" % (rsid, pdbr.resname, refr.resname))
        pdbr.nter = refr.nter
        pdbr.cter = refr.cter
        pdbr.topology = refr.topology

def get_topologies(toplist):
    topfiles = [ currdir + "/oplsx-top.json",
                  currdir + "/dna-rna-top.json",
                ]

    for f in toplist:
        print(f, file=sys.stderr)
        assert os.path.exists(f), f
        topfiles.append(f)

    topologies = []
    for f in topfiles:
        try:
            topologies.append(topology.load(json.load(open(f))))
        except:
            print(f, file=sys.stderr)
            raise
    top_residues, top_patches = topology.merge(topologies)
    return top_residues, top_patches

def run_pdbcomplete(pdbfile, args):
    pp("run pdbcomplete")
    mutations = {}
    for f in args.mutatefiles:
        assert os.path.exists(f), f
        for l in open(f):
            ll = l.split()
            if len(ll) == 0: continue
            assert len(ll) == 2, l
            mutations[ll[0]] = ll[1]

    top_residues, top_patches = get_topologies(args.topfile)
    pdb, ter_indices = read_pdb(open(pdbfile), pdbfile, top_residues,
                    args.dna, args.rna, args.heavy, mutations,
                    add_termini = args.termini, modbase=args.modbase, modres=args.modres)
    is_prot, is_na = eval_moltype(pdb)
    pdblines, mapping, _ = write_pdb(pdb, args)
    termini_pdb(pdb, args.nter, args.cter)
    patches = {}
    if args.readpatch:
        indata = open(os.path.splitext(pdbfile)[0]+'.patch').readlines()
        indata = [line.split() for line in indata]
        args.patches = indata

    for p in args.patches:
        chain, resid = p[0], p[1]
        print(chain, file=sys.stderr)
        if chain != "None":
            resindices = [ri for ri,r in enumerate(pdb) if r.resid.strip() == resid and r.chain.strip() == chain]
        else:
            resindices = [ri for ri,r in enumerate(pdb) if r.resid.strip() == resid]
        if len(resindices) == 0:
            raise ValueError("No residues have chain %s resid %s" %(chain, resid))
        elif len(resindices) > 1:
            raise ValueError("Multiple residues have resid %s" % resid)
        resid2 = pdb[resindices[0]].resid
        if resid2 not in patches: patches[resid2] = []
        pname = p[-1].lower()
        if pname == "none": pname = None
        patches[resid2].append(pname)
    patch_pdb(pdb, patches, top_patches)

    if args.refe:
        refe, _ = read_pdb(open(args.refe), args.refe, top_residues,
                        args.dna, args.rna, args.heavy,
                        mutations, add_termini=args.termini)
        patch_pdb(refe, patches, top_patches)
        if not args.heavy:
            update_patches(refe, top_patches)
        set_reference(pdb, refe)

    if is_na:
        libname = "rnalib"
        if args.dna:
            libname = "dnalib"
        nalib = load_nalib(libname)
        print('apply_nalib', file=sys.stderr)
        apply_nalib(pdb, nalib, args.manual)

    # If hydrogen must be added to DNA/RNA, pdb2pqr must be used after nalib
    # If protein in input, apply pdb2pqr
    if is_prot or not args.heavy:
        print('apply_pdb2pqr', file=sys.stderr)
        pdblines, mapping, mutations_5ter = write_pdb(pdb, args, patches, write_missing=False, pdb2pqr=True)
        with open("/tmp/kakapo", "w") as kakapo:
            for l in pdblines:
                print(l,file=kakapo)
        pqrlines = run_pdb2pqr(pdblines)
        with open("/tmp/orca", "w") as orca:
            for l in pqrlines:
                print(l,end='',file=orca)
        pqr, ter_indices = read_pdb(pqrlines, "<PDB2PQR output from %s>" % pdbfile, top_residues,
                            args.dna, args.rna, args.heavy, mutations=mutations_5ter)
        pdbcomplete(pdb, pqr)
        if not args.heavy and not args.refe:
            update_patches(pdb, top_patches)

    if args.refe:
        pdbfix(pdb, refe)

    if not args.manual:
        pdb_lastresort(pdb, args.heavy)
        check_pdb(pdb, heavy=args.heavy)

    pdblines, _ , _ = write_pdb(pdb, args, patches, heavy=args.heavy, write_missing=False)
    return pdblines, mapping, pdb
