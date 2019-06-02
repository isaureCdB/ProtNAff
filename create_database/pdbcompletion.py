#!/usr/bin/env python3
# Copyright Sjoerd J. De Vries (TUM), Isaure Chauvot de Beauchene (CNRS)

"""
Complete missing atmos from a PDB. Adapted from ATTRACT's aareduce.py script.
Calls pdbcomplete.py

For nucleic acids, the patch "5ter" is applied by default on the 1st nucleotide,
which removes a 5PHO if present in the input. To prevent this, use
"--patch [chain] [X] None" with X the residue id of the 1st nucl; this will keep the
phosphate, but not the O5T. To keep the O5T as well, add --termini or --nter.
"""

#TODO: change clash cutoff in apply_nalib when hydrogen added to lib
#TODO: change the name of 1st phosphate into "5PH0" for hydrogen completion
#        of dna/rna by pdb2pqr, if --patch <first res> None

from __future__ import print_function
import sys, os, json
import pdbcomplete
import topology
from pdbcomplete import run_pdb2pqr, pdbfix, update_patches, pdb_lastresort, FixError, run_pdbcomplete
from copy import deepcopy

has_argparse = False
try:
    import argparse
    has_argparse = True
except ImportError:
    import optparse  #Python 2.6

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

def pp(i):
    print(i, file=sys.stderr)

def read_filelist(filelist):
    ret = []
    for l in open(filelist):
        l = l.strip()
        if not len(l): continue
        assert len(l.split()) == 1, (filelist, l)
        ret.append(l)
    return ret

def check_models(pdbs):
    atoms0 = [l[11:20] for l in pdbs[0]]
    for pdb in pdbs:
        atoms = [l[11:20] for l in pdb]
        assert atoms == atoms0, "All MODELS must have the same atoms in same order"

def read_models(pdb):
    outp_name = None
    pdbs = []
    pdbs_list = []
    for l in open(pdb):
        ll = l.split()
        if not pdbs and ll[0] != "MODEL":
            # there are no MODEL statements
            return [pdb]
        if ll[0] == "MODEL":
            if outp_name is not None:
                outp.close()
                pdbs_list.append(outp_name)
            outp_name = "/tmp/model%s.pdb"%ll[1]
            outp = open(outp_name, "w")
            pdbs.append([])
        if ll[0] not in ["ATOM","TER","HETATM"]: continue
        if l0 == "HETATM" and not args.modres and not args.modbase:
            print("discarde %s"%l, file=sys.stderr) ###
            continue
        print(l, file=sys.stderr) ###
        outp.write(l)
        pdbs[-1].append(l)
    outp.close()
    pdbs_list.append(outp_name)
    #if len(pdbs):
    #    check_models(pdbs)
    return pdbs_list

currdir = os.path.abspath(os.path.split(__file__)[0])
topfiles = [ currdir + "/oplsx-top.json",
              currdir + "/dna-rna-top.json",
            ]

map_atnames = {
    "HO2'": "H2''",
    "H5'" : "H5''"
    }
parser = argparse.ArgumentParser(description=__doc__,
                          formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("pdb",help="PDB file to reduce")
parser.add_argument("output",help="all-atom reduced output PDB file", nargs="?")
parser.add_argument("--heavy",help="Ignore all hydrogens", action="store_true")
parser.add_argument("--refe", "--reference",help="Analyze the hydrogens of a reference file to determine histidine/cysteine states")
parser.add_argument("--autorefe",help="Analyze the hydrogens of the input PDB to determine histidine/cysteine states", action="store_true")
parser.add_argument("--dna",help="Automatically interpret nucleic acids as DNA", action="store_true")
parser.add_argument("--rna",help="Automatically interpret nucleic acids as RNA", action="store_true")
parser.add_argument("--dna_chain",help="chain IDs for DNA", nargs='+')
parser.add_argument("--rna_chain",help="chain IDs for RNA", nargs='+')
parser.add_argument("--termini",help="Add N/C-ter (5'/3'-ter for nucleic acids) to each chain", action="store_true")
parser.add_argument("--nter", "--nterm" , dest="nter", help="Add N-ter/5'-ter to the specified residue number",
                    action="append", type=int, default=[])
parser.add_argument("--cter","--cterm", dest="cter", help="Add C-ter/3'-ter to the specified residue number",
                    action="append", type=int, default=[])
parser.add_argument("--startres", help="Set residue number of first residue", type=int, default=1)
parser.add_argument("--manual",help="""Enables manual mode.
    In automatic mode (default), in case of missing atoms, a number of last-resort fixes are attempted that add pseudo-hydrogens
    at the position of its connected heavy atom. If there too many missing atoms, the residue is discarded and a warning is given.
    In manual mode, last-resort fixes are disabled, and missing atoms are simply printed with XXXXXXX in their coordinates.
    """, action="store_true")
parser.add_argument("--top", "--topfile",dest="topfile",help="Additional topology file in CNS format that contains additional user-defined atom types (e.g. modified amino acids)", action="append",default=[])
parser.add_argument("--patch", help="Provide a chain, a residue number and a patch name to apply (ex: [A 1 Nter], [None 1 None])",
                    dest="patches", nargs=3, action="append",default=[])
#parser.add_argument("--keepchains", help="Keep the chain IDs of the input PDB", default=" ")
parser.add_argument("--no_renum_res", help="renumber residues from 1", action="store_true")
parser.add_argument("--mutate", help="Provide a 2-column residue mutation file",
                    dest="mutatefiles", action="append",default=[])
parser.add_argument("--modres", help="Interpret HETATM records as ATOM if they have a protein backbone", action="store_true")
parser.add_argument("--modbase", help="Interpret HETATM records as ATOM if they have at least three sugar atoms", action="store_true")
parser.add_argument("--batch", help="run aareduce in batch mode. Input and output must be lists of PDBs", action="store_true")
parser.add_argument("--dumppatch", help="Dump all applied patches to a file", action="store_true")
parser.add_argument("--readpatch", help="Read previously applied patches from a file (requires converted input pdb)", action="store_true")

args = parser.parse_args()

if args.rna and args.dna:
    raise ValueError("--dna and --rna are mutually incompatible")

if args.heavy and (args.autorefe or args.refe):
    raise ValueError("--(auto)refe and --heavy are mutually incompatible")

if args.autorefe and args.refe:
    raise ValueError("--autorefe and --refe are mutually incompatible")
if args.autorefe:
    args.refe = args.pdb

if args.readpatch and len(args.patches):
    raise ValueError("--readpatch and explicit patch specification via --patch are mutually incompatible")

for f in args.topfile:
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

if args.batch:
    infiles = read_filelist(args.pdb)
    for f in infiles:
        assert os.path.exists(f), f
    if args.output:
        outfiles = read_filelist(args.output)
    else:
        outfiles = [f.split(".pdb")[0] + "-aa.pdb" for f in infiles]
    for pdb, outfile in zip(infiles, outfiles):
        pdblines, mapping, pdbtop = run_pdbcomplete(pdb, args)
        outf = open(outfile, "w")
        for l in pdblines:
            print(l, file = outf)
        outf.close()
        if args.dumppatch:
            outfilep = os.path.splitext(outfile)[0] + ".patch"
            outf = open(outfilep, "w")
            for i,res in enumerate(pdbtop):
                if len(res.topology.patches):
                    for p in res.topology.patches:
                        outf.write(str(i+args.startres)+' '+p+'\n')
            outf.close()
else:
    outfile = os.path.splitext(args.pdb)[0] + "-aa.pdb"
    if args.output is not None:
        outfile = args.output
    outf = open(outfile, "w")
    pdbs = read_models(args.pdb)
    pp(pdbs)
    m = 0
    for pdb in pdbs:
        m += 1
        pdblines, mapping, pdbtop = run_pdbcomplete(pdb, args)
        if len(pdbs) > 1:
            print('MODEL %i'%m, file = outf)
        for l in pdblines:
            print(l, file = outf)
        if len(pdbs) > 1:
            print('ENDMDL', file = outf)
        if m == 1:
            for v1, v2 in mapping:
                print(v1, v2)
        if args.dumppatch and m == 1:
            outfilep = os.path.splitext(outfile)[0] + ".patch"
            outf2 = open(outfilep, "w")
            for i,res in enumerate(pdbtop):
                if len(res.topology.patches):
                    for p in res.topology.patches:
                        outf2.write(str(i+args.startres)+' '+p+'\n')
            outf2.close()
    outf.close()
