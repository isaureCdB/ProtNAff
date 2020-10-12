#!/usr/bin/env python3
import json, sys, os
import argparse
import shutil #package for managing files
from multiprocessing import Pool
from itertools import repeat
#import the query tou want to use
from queries.query_ss import query_ss

#import of some useful fonctions
from extract_chains import extract_chains
from create_pdb import create_pdb

'''
Automatization for the creation of trinucleotide fragments.

Steps you have to follow to use this file:
    - First you need to find the query you will use
    - Then import your querry from the package queries
    - Change the query name line XXX
    - Then launch with python 3
'''

###################################################
###################################################
################                   ################
################     ARGPARSE      ################
################                   ################
###################################################
###################################################

parser = argparse.ArgumentParser(description=__doc__,
                          formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("file", type=str, help="The json file on which you want to find ssRNA and dsRNA.")
parser.add_argument("runName", type=str, help="The name of the run.")
parser.add_argument("output", type=str, help="The directory where you want your output.")
args = parser.parse_args()

# Put input in variable names
json_file=args.file #x3dna.json
outpdir = args.output # ~/projets/ssRNA/benchmark/
runName = args.runName

###################################################
###################################################
################                   ################
################     FUNCTIONS     ################
################                   ################
###################################################
###################################################

def create_folder(path, folder):
    """
    This function create a new folder in the path. It checks if it exists or not.
    Input : - path, where to create the folder
            - folder, the name of the new folder
    Output: The new path to this folder
    """
    output = '{}/{}'.format(path, folder)
    if not os.path.exists(output):
        os.mkdir(output)
    return output

def write_file(path, prot_rna, aa_r, list_created_file):
    """
    This function is merging in one file.pdb all chains involved in the
    biological assembly
    Input : - path, the path of the biological assembly
            - prot_rna, if it is a Prot or a RNA
            - aa_r, if files to merge are all atoms or reduced
            - list_created_file, the list of files created by the function create_pdb
    """
    with open("{}/{}-{}.pdb".format(path, prot_rna, aa_r), 'w') as f:
        for created_file in list_created_file:
            if created_file.startswith(prot_rna):
                shutil.copyfileobj(open("{}/{}-{}.pdb".format(path, created_file.split(".")[0], aa_r), 'r'), f)

def delete_files(path, list_created_file):
    """
    This function deletes files which were merged by the write_file function
    Input : - path, the path of the biological assembly
            - list_created_file, the list of files created by the function create_pdb
    """
    for created_file in list_created_file:
        os.remove("{}/{}-aa.pdb".format(path, created_file.split(".")[0]))
        os.remove("{}/{}-aar.pdb".format(path, created_file.split(".")[0]))

def write_fragments_file(path, aa_r, nuclfrag, chain_id):
    """
    This function write file with nucleotides from the query
    Input : - path, the path of the biological assembly
            - aa_r, if files to merge are all atoms or reduced
    """

    chain_id = chain_id.split('_')[-1]
    rnafile = open('{}/RNA_{}_fragments-{}.pdb'.format(path, chain_id, aa_r), 'w')

    for l in open('{}/../RNA-{}.pdb'.format(path, aa_r), 'r'):
        if l.startswith("TER") and l[21] == chain_id:
            rnafile.write(l)
        if l.startswith("ATOM") or l.startswith("HETATM"):
            if int(l[22:26]) in nuclfrag and l[21] == chain_id:
                rnafile.write(l)
    rnafile.close()


def script_for_parallel(pdb_id, tmpPDB, outpdir, runName, js):
    try:

        pdb_info = js[pdb_id]

        #To make the difference between NMR and crystallo complexes
        nmr_model = "0"

        #TODO manage to find the best model in all NMR models
        if pdb_info["Nmodels"] > 1:
            best = "1"
            nmr_model = best

        #Creation of the directory linked to the PDB complex

        dirPDBName = create_folder(outpdir, pdb_id)

        #Creation of the dictionary chains containing for each biological assembly
        #every chains (prot and rna) involved
        chains = extract_chains(pdb_id, tmpPDB, True)

        #If there is not biological assembly from the PDB
        if not chains:
            assembly_bool = False
            #Then download the .pdb file
            os.system('./pdb_download_gz {} {}'.format(pdb_id, tmpPDB))

            #chains involved are those from the json
            protein_list = pdb_info["protchains"]
            rna_list = pdb_info["nachains"]

            #Creation of the dictionary chains
            chains[0] = {"protein":protein_list, "rna":rna_list}
        else:
            assembly_bool = True

        #for every biological assembly
        for assembly in chains.keys():

            #Creation of one folder for each biological assembly
            folder_name = "chains_{}_{}".format('-'.join(chains[assembly]["protein"]), '-'.join(chains[assembly]["rna"]))
            dirBiologicalAssembly = create_folder(dirPDBName, folder_name)

            if assembly_bool:
                pdb_name = '{}.pdb{}'.format(pdb_id, assembly)
            else:
                pdb_name = '{}.pdb'.format(pdb_id)

            #Creation of files.pdb which we will use later, and a list of which
            # files are created.
            list_created_file = create_pdb(tmpPDB, dirBiologicalAssembly, pdb_name, chains[assembly], nmr_model)

            with open("{}/file_to_reduce_{}.list".format(dirBiologicalAssembly, runName), 'w') as f:
                for file in list_created_file:
                    f.write("{}\n".format(file))

            #Creation of all atoms files and reduced files for chains.pdb
            os.system('./aareduce_benchmark {}/ {}'.format(dirBiologicalAssembly, runName))

            #Merging of Prot/RNA_chain-aa/aar.pdb files
            write_file(dirBiologicalAssembly, "RNA", "aa", list_created_file)
            write_file(dirBiologicalAssembly, "RNA", "aar", list_created_file)
            write_file(dirBiologicalAssembly, "Prot", "aa", list_created_file)
            write_file(dirBiologicalAssembly, "Prot", "aar", list_created_file)

            delete_files(dirBiologicalAssembly, list_created_file)

            #Creation of the folder for the run (or query it is the same here)
            dirNameRun = create_folder(dirBiologicalAssembly, runName)

            #Then do the query for every RNA chains in the biological assembly
            for chain_id in chains[assembly]["rna"]:
                length = 3
                nuclfrag = query_ss(pdb_info, "chain_{}".format(chain_id), length)
                # If the query return nothing: continue
                if len(nuclfrag) == 0:
                    print("No structure of interest in {}\n\n".format(pdb_id))
                    continue
                print('{} {} : {}'.format(pdb_id, chain_id, nuclfrag), end=' ')
                print('')

                #Write the pdb name if the query is successful
                with open("{}/query_successful_{}.txt".format(outpdir, runName), 'a') as f:
                    f.write("{}\n".format(pdb_id))

                # Write the file with RNA information
                write_fragments_file(dirNameRun, "aa", nuclfrag, chain_id)
                write_fragments_file(dirNameRun, "aar", nuclfrag, chain_id)

#                 #Creation of fragments for the docking
#                 os.system("extractfrag.sh rna {}".format(dirNameRun))

    #Catch the problem, more or less...
    except Exception as e:
        print("Erreur sur le PDB {}".format(pdb_id))
        with open("{}/pdb_with_error_{}.txt".format(outpdir, runName), 'a') as f:
            f.write("{} {}\n".format(pdb_id, e))

    print("Done for {}\n\n".format(pdb_id))
    return None

###################################################
###################################################
################                   ################
################      SCRIPT       ################
################                   ################
###################################################
###################################################

js = json.load(open(json_file))

# Creation of a tmp/ directory on which every pdb will be downloaded
# This directory is deleted at the end
# If you will use this script often, comment the deletion of this directory
tmpPDB = "{}/tmp/".format(outpdir)
if not os.path.exists(tmpPDB):
    os.mkdir(tmpPDB)

pool = Pool(processes=1)

pool.starmap(script_for_parallel, zip(js.keys(), repeat(tmpPDB), repeat(outpdir), repeat(runName), repeat(js)))

#To remove the tmp directory, comment if you often use the benchmark.
#shutil.rmtree(tmpPDB)
