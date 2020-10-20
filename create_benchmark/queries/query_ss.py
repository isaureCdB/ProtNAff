#!/usr/bin/env python3

import json

"""
A query to use in the benchmark script, it selects only ssRNA in interaction with protein.
"""

def query_ss(pdb_info, chain_id, length):
    """
    Input: pdb_info, information of the pdb in json format
           chain_id, the name of the RNA chain
           length, the length of the sequence selected
    Output: nuclfrag, a set of nucleotide of interest
    """

    ss_set = set(["L", "T", "S", "J", "B", "I"])

    nuclfrag = set()
    # print(pdb_info['ss'][chain_id])
    ss = pdb_info['ss'][chain_id]

    interf_prot = pdb_info["interface_protein"]["model_1"][chain_id]
    # print(interf_prot)
    for element in interf_prot:
        element = element.split("_")[1]

    # List of interface nucleotide resid
    nucl_interf = list(interf_prot.keys())
    nucl_interf = [element.split("_")[1] for element in nucl_interf]
    print(nucl_interf)
    # List of single-stranded nucleotides resid
    nucl_ss = [n.split("_")[1] for n in ss if ss[n][0] in ss_set]
    print(nucl_ss)
    # List of single-stranded interface nucleotide resid
    nucl_interf_ss = [int(n) for n in nucl_interf if n in nucl_ss]
    nucl_interf_ss.sort()
    print(nucl_interf_ss)
    if len(nucl_interf_ss) < length:
        return nuclfrag
    for i in range(1, len(nucl_interf_ss) - 1):
        if nucl_interf_ss[i+1] - nucl_interf_ss[i] == 1:
            nuclfrag.add(nucl_interf_ss[i])
            nuclfrag.add(nucl_interf_ss[i]+1)
        elif nucl_interf_ss[i+1] - nucl_interf_ss[i] == 2:
            nuclfrag.add(nucl_interf_ss[i])
            nuclfrag.add(nucl_interf_ss[i]+1)
            nuclfrag.add(nucl_interf_ss[i]+2)
        elif nucl_interf_ss[i+1] - nucl_interf_ss[i] == 3:
            nuclfrag.add(nucl_interf_ss[i])
            nuclfrag.add(nucl_interf_ss[i]+1)
            nuclfrag.add(nucl_interf_ss[i]+2)
            nuclfrag.add(nucl_interf_ss[i]+3)
    if nucl_interf_ss[1] - nucl_interf_ss[0] <= 3:
        nuclfrag.add(nucl_interf_ss[0])
        nuclfrag.add(nucl_interf_ss[0]+1)
        nuclfrag.add(nucl_interf_ss[0]+2)

    nuclfrag = list(nuclfrag)
    nuclfrag.sort()
    nuclfrag = keep_length(nuclfrag, length)

    # nuclfrag_mapped = []
    # for nucl in nuclfrag:
    #     nuclfrag_mapped.append(int(pdb_info["mapping"][chain_id][str(nucl)]))

    return nuclfrag


def keep_length(nuclfrag, length):

    i = 0
    while i < len(nuclfrag) - 1:
        start = i
        j = i
        while j < len(nuclfrag) - 1 and nuclfrag[j+1] - nuclfrag[j] == 1:
            j += 1
        if j - start < length:
            nuclfrag = nuclfrag[:start] + nuclfrag[j+1:]
        else:
            i = j + 1
    return nuclfrag

def main():
    js = json.load(open("/home/amoniot/Documents/libraries/nalib_oct2019_new/structures.json"))
    print(query_ss(js["1B23"], "chain_R", 3))
    # for pdb_id in js:
    #     pdb_info = js[pdb_id]
    #     for chain_id in pdb_info['nachains']:
    #         nuclfrag = query_ss(pdb_info, chain_id, 5)
    #         print("{} {} {}".format(pdb_id, chain_id, nuclfrag))

if __name__ == '__main__':
    main()
