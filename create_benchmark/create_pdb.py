#!/usr/bin/env python3

import os, sys

def write_file(input_dir, output_dir, pdb_name, chain, prot_rna, list_hetatm):
    file_std = open('{}/{}_{}.pdb'.format(output_dir, prot_rna, chain), 'w')

    for l in open('{}/{}'.format(input_dir, pdb_name), 'r'):
        if l.startswith("TER") and l[21] == chain:
            file_std.write(l)
        if l.startswith("ATOM"):
            if l[21] == chain:
                file_std.write(l)
        elif l.startswith("HETATM") and prot_rna != "Prot":
            if l[21] == chain and l[17:20] in list_hetatm:
                file_std.write(l)
        if l.startswith("MODEL        2"):
            break

def write_file_nmr(input_dir, output_dir, pdb_name, chain, prot_rna, nmr_model, list_hetatm):
    file_nmr = open('{}/{}_{}.pdb'.format(output_dir, prot_rna, chain), 'w')
    find = False
    for l in open('{}/{}'.format(input_dir, pdb_name), 'r'):
        if l.startswith("MODEL"):
            if int(l[11:14]) == eval(nmr_model):
                find = True
        if find:
            if l.startswith("ENDMDL"):
                break
            else:
                if l.startswith("TER") and l[21] == chain:
                    file_nmr.write(l)
                if l.startswith("ATOM"):
                    if l[21] == chain:
                        file_nmr.write(l)
                elif l.startswith("HETATM") and prot_rna != "Prot":
                    if l[21] == chain and l[17:20] in list_hetatm:
                        file_std.write(l)

def create_pdb(input_dir, output_dir, pdb_name, rna_chains, nmr_model):

    list_hetatm = []
    with open("/home/amoniot/attract/allatom/rna-mutate.list", 'r') as f:
        lines = f.readlines()
        for line in lines:
            list_hetatm.append(line.split()[0])

    list_created_file = []

    for rna_chain in rna_chains["rna"]:
        list_created_file.append('RNA_{}.pdb'.format(rna_chain))
        if eval(nmr_model) > 0:
            write_file_nmr(input_dir, output_dir, pdb_name, rna_chain, "RNA", nmr_model, list_hetatm)
        else:
            write_file(input_dir, output_dir, pdb_name, rna_chain, "RNA", list_hetatm)

    for prot_chain in rna_chains["protein"]:
        list_created_file.append('Prot_{}.pdb'.format(prot_chain))
        if eval(nmr_model) > 0:
            write_file_nmr(input_dir, output_dir, pdb_name, prot_chain, "Prot", nmr_model, list_hetatm)
        else:
            write_file(input_dir, output_dir, pdb_name, prot_chain, "Prot", list_hetatm)

    return list_created_file
