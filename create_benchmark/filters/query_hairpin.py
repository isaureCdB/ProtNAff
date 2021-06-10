#!/usr/bin/env python3

import json
import argparse

def get_ss(pdb_info, chain_id, interact, nucl_ss):
    """
    INPUT : - pdb_info : information of the pdb on json format
            - chain_id : the name of the RNA chain
            - interact : boolean, True if ss part have to be in contact with protein
            - nucl_ss : the list of nucleotides which are single stranded
    OUTPUT: - a list of single stranded nucleotides, in contact or not
    """
    if interact:
        try:
            interf_prot = pdb_info["interface_protein"]["model_1"][chain_id]
            nucl_interf = list(interf_prot.keys())
            nucl_interf = [int(element.split("_")[1]) for element in nucl_interf]
            nucl_interf_ss = []
            for n in nucl_interf:
                if n in nucl_ss:
                    test = True
                    for key in interf_prot["res_"+str(n)].keys():
                        if interf_prot["res_"+str(n)][key] == 999:
                            test = False
                            break
                    if test:
                        nucl_interf_ss.append(int(n))
        except:
            nucl_interf_ss = []
        return nucl_interf_ss
    else:
        return nucl_ss


def check_ds_binding(pdb_id, chain_id, dir_3dna):
    """
    This function is checking in 3dna outputs the couples of double stranded
    based pairs.
    INPUT : - pdb_id : The name of the pdb_id
            - chain_id : The name of the chain of RNA we are looking in the pdb
            - dir_3dna : the folder which contains 3dna outputs
    OUTPUT: - couples : a list of couples
    """
    chain_id = chain_id.split('_')[1]
    json_file = json.load(open("{}/{}-1-dssr.json".format(dir_3dna, pdb_id)))
    couples = []
    try:
        pairs = json_file["pairs"]

        for pair in pairs:
            nt1 = pair["nt1"]
            nt2 = pair["nt2"]
            num_nt1 = int(nt1.split('.')[3])
            chain_nt1 = nt1.split('.')[1]
            num_nt2 = int(nt2.split('.')[3])
            chain_nt2 = nt2.split('.')[1]
            if chain_nt1 == chain_id and chain_nt2 == chain_id:
                couples.append([num_nt1, num_nt2])
    except:
        pass
    return couples


def query_hairpin(pdb_id, pdb_info, chain_id, dir_3dna, interact, length_ds, length_ss):
    """
    The main function, which is doing a selection of hairpins in the pdb file
    given in input.
    INPUT : - pdb_id : The name of the pdb id
            - pdb_info : information of the pdb in json format
            - chain_id : the name of the RNA chain
            - dir_3dna : the folder which contains 3dna outputs
            - interact : boolean, True if the ss part have to be in contact
            - length_ds : the length of the ds part
            - length_ss : the length of the ss part
    OUTPUT: - results : a list of hairpin [res_start, res_end]
    """

    # Define what is ds and what is ss
    ds_set = set(["D"])
    ss_set = set(["L", "T", "S", "J", "B", "I"])

    nuclfrag = set()
    secondary_structure = pdb_info['ss'][chain_id]

    nucl_ds = [n.split('_')[1] for n in secondary_structure if secondary_structure[n][0] in ds_set]
    nucl_ss = [int(n.split('_')[1]) for n in secondary_structure if secondary_structure[n][0] in ss_set]

    nucl_ss = get_ss(pdb_info, chain_id, interact, nucl_ss)

    # If there is less ss nucleotide than the minimum size of the ss part, then exit
    if len(nucl_ss) < length_ss:
        return nuclfrag

    list_ds = []
    for x in nucl_ds:
        try:
            list_ds.append([int(x), "ds"])
        except:
            pass
    list_ss = []
    for x in nucl_ss:
        try:
            list_ss.append([int(x), "ss"])
        except:
            pass

    list_all = list_ds + list_ss
    list_all = sorted(list_all, key=lambda l:l[0], reverse=False)

    list_hairpin = []
    i = 0

    # main loop, a bit tricky, checking the sequence of nucleotides to find hairpins
    while i < len(list_all):
        if list_all[i][1] == "ds": #Start if the nucleotide is ds
            # Define variables
            change1 = False
            change2 = False
            ds1 = 0
            ds2 = 0
            ss = 1
            end_ss = None
            hairpin = True
            j = i
            # Second loop
            while j < len(list_all) -1 and hairpin == True:
                # if there is no gap between nucleotides
                if list_all[j+1][0] - list_all[j][0] == 1:
                    if change2:
                        ds2 += 1
                        if ds2 >= length_ds:
                            break
                        if list_all[j][1] != list_all[j+1][1]:
                            break
                    else:
                        if change1:
                            if list_all[j][1] != list_all[j+1][1]:
                                change2=True
                                end_ss = j
                            else:
                                ss += 1
                        else:
                            if list_all[j][1] != list_all[j+1][1]:
                                change1=True
                            else:
                                if ds1 >= (length_ds - 1):
                                    i += 1
                                else:
                                    ds1 += 1
                else:
                    hairpin = False
                j += 1
                # print(i, j, change1, change2, ds1, ds2, ss)
            if not(change2 and ds1 >= (length_ds - 1) and ds2 >= length_ds and ss >= length_ss):
                hairpin = False
            if hairpin:
                list_hairpin.append([list_all[i][0],list_all[j][0]])
            if end_ss != None:
                i = end_ss
            else:
                i = j + 1
        else:
            i +=1

    couples = check_ds_binding(pdb_id, chain_id, dir_3dna)
    results= []
    for hairpin in list_hairpin:
        test_hairpin = True
        for i in range(length_ds):
            # If the ds part is not paired to the one in front, it is not a hairpin
            if [hairpin[0] + i, hairpin[1] - i] not in couples:
                test_hairpin = False

        if test_hairpin:
            result = [pdb_info["mapping"][chain_id][str(hairpin[0])], pdb_info["mapping"][chain_id][str(hairpin[1])]]
            results.append([result,pdb_info["sequence"][chain_id][(hairpin[0]-1):hairpin[1]]])
    return results


def check_ds_binding2(pdb_id, chain_id):
    """

    Old version, not useful anymore.

    This function is checking in 3dna outputs the couples of double stranded
    based pairs.
    INPUT : - pdb_id : The name of the pdb_id
            - chain_id : The name of the chain of RNA we are looking in the pdb
    OUTPUT: - couples : a dictionnary of couples
    #TODO : switch from dictionnary to set
    """
    couples = {}
    seuls = {}
    for key in pdb_info["intraNA_hb"][chain_id]:
        nb = key.split('_')[1]
        for key2 in pdb_info["intraNA_hb"][chain_id][key]:
            if key2 == 'base':
                for key3 in pdb_info["intraNA_hb"][chain_id][key]["base"]:
                    if key3 == 'other':
                        for element in pdb_info["intraNA_hb"][chain_id][key]["base"]["other"]:
                            distance = element[0]
                            for key4 in seuls:
                                if distance in seuls[key4]:
                                    couples["{}-{}".format(key, key4)] = 1
                            if key in seuls:
                                seuls[key].append(distance)
                            else:
                                seuls[key] = [distance]
    return couples

def main():
    a = argparse.ArgumentParser(prog="query_hairpin.py")
    a.add_argument("folder_3dna", help="The folder which contains the 3dna outputs.")
    a.add_argument("structures", help="The structures.json file which is created by the create_database.sh of NAfragDB")
    a.add_argument("output", help="The path of the output file")
    a.add_argument("--ss",type=int, help="The length of the single stranded part.")
    a.add_argument("--ds", type=int, help="The length of the double stranded part.")
    a.add_argument("--contact", type=int, help="True if the ss part have to be in contact, False if you do not care.")
    args = a.parse_args()

    dir_3dna = args.folder_3dna
    js = json.load(open(args.structures))

    length_ds = args.ds
    length_ss = args.ss
    contact = args.contact
    if contact == 1:
        contact = True
    else:
        contact = False
    towrite = []
    for pdb_id in js:
        pdb_info = js[pdb_id]
        print(pdb_id)
        # if pdb_id == "1ASY":
        for chain_id in pdb_info['nachains']:
            nuclfrag = query_hairpin(pdb_id, pdb_info, "chain_{}".format(chain_id), dir_3dna, contact, length_ds, length_ss)
            towrite.append([pdb_id, chain_id, nuclfrag])
    with open(args.output, 'w') as f:
        for hairpin in towrite:
            for nucl in hairpin[2]:
                f.write("{} {} {} {}\n".format(hairpin[0], hairpin[1], nucl[0], nucl[1]))

    print("Done")
if __name__ == '__main__':
    main()
