#!/usr/bin/env python3

import json
import argparse
import copy

def non_modified_trinucl(all_frag, nb_frag, js_struc):
    """
    Input:
    Output:
    """
    frag = all_frag[nb_frag]
    pdb_id = frag["structure"]
    chain_id = "chain_" + frag["chain"]
    res_ids = frag["resid"]

    result = False
    if sum(frag["missing_atoms"]) == 0:
        for res in res_ids:
            try:
                # print(res, chain_id, js_struc[pdb_id]["canonized"])
                if res in js_struc[pdb_id]["canonized"][chain_id]:
                    # print("bouh")
                    result = True
                    break
            except:
                pass
    else:
        result = True

    return result


def main():
    a = argparse.ArgumentParser(prog="filter_no_modified.py")
    a.add_argument("structures", help="The structures.json file which is created by the create_database.sh of ProtNAff")
    a.add_argument("fragments", help="The fragments.json file which is created by the create_database.sh of ProtNAff")
    a.add_argument("output", help="The path of the output file")
    args = a.parse_args()

    js_struc = json.load(open(args.structures))
    js_frag = json.load(open(args.fragments))

    copy_js_frag = copy.deepcopy(js_frag)

    for seq in js_frag.keys():
        all_frag = js_frag[seq]
        for nb_frag in all_frag.keys():
            results = non_modified_trinucl(all_frag, nb_frag, js_struc)
            if results:
                print("bouh")
                copy_js_frag[seq].pop(nb_frag)

    with open('{}'.format(args.output), 'w') as fp:
        json.dump(copy_js_frag, fp, indent = 2, sort_keys = True)

if __name__ == '__main__':
    main()
