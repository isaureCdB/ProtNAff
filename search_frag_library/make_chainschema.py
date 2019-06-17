#!/usr/bin/env python3

import json
# TODO: Documentation

'''
    Queries on the fragment library use 3 dictionaries:
    _ the data dictionnary, provided by structures.json (written when creating the database)
    _ the chainschema (see below) that describes the format of the data and how it can be queried
    _ the query "variables" dictionnary (see query.py), that contains description of the subdata we are interested in

    The data dictinonary contains the datakeys. For example "1B7F" is a datakey that describes a PDB code, and "bptype"
    is a datakey that describe the types of base-pairs in that PDB structure.

    chainschema contains schemakeys and under each schemakey, the subschema is stored that describes the
    subdata for the corresponding datakey. For example, for "pdbchains", the subschema is [], which means
    that the subdata data["1B7F"]["pdbchains"] is a list of elementary values (here strings that are chain IDs)

    The data and the chainschema are both multi-level (nested) dictionaries, but a datakey or schemakey is only for data at that level.

    A schemakey can be:
    _ a constant schemakey (e.g. "bptype"), that correspond exactly to the datakey of the same name, of the data at that level.
    _ a variable schemakey (e.g. the PDB code of the structure), starting with a star,
      that correspond to all the datakeys of the data at that level.

    ex;
    The constant schemakey "bptype" describes the data as a dictionnary containing datakey "bptype".
    Its subschema (which is {"?*chain": {"?*res": []} } ) applies only to (the subdata under) the "bptype" datakey of the data.

    The variable schemakey "*pdbcode" describes the data as a dictionnary containing one or more datakeys.
    The subschema of "pdbcode" applies to (the subdata under) each datakey of the data.

    Questionmark distinguish schemakeys that are non-obligatory. If they are not found, the query returns None.

    The following variable schemakeys are being used:
    chain = chain ID
    res = residue index
    part = part of the nucleotide ("ph" for phosphate group, "sugar", "base")
    pos = position in the NA sequence toward a given nucl at position n ("n-1", "n+1"...)

    Variable schemakeys (with *) correspond to querykeys in the query dictionnary named "variables" in query.py
    To query a constant schemakey, the schema contains at the end an entry ("=>", <querykey>) that describes which
    constant schemakey is to be queried.
    For example, if we are interested in data["1B7F"]["interface_protein"]["A"]["111"]["base"]:
    - Level 1: a variable schemakey called "pdbcode".
        query dictionary: {"pdbcode": "1B7F"}
    - Level 2: constant schemakey. ("=>","name") means that the constant schemakey name must be defined under "name"
        query dictionary: {"name": "interface_protein"}
    - Level 3:  a variable schemakey called "chain".
        query dictionary: {"chain": "A"}
    - Level 4:  a variable schemakey called "res".
        query dictionary: {"res": "111"}
    - Level 5:  a variable schemakey called "part".
        query dictionary: {"part": "base"}
    - Level 6: None, meaning that a str/float/int is returned (in this case, a float describing the distance).
    If it was [], a list would be returned, and [[]] for a list-of-lists.
    Final query dict: {
        "pdbcode": "1B7F",
        "name": "interface_protein",
        "chain": "A",
        "res": "111",
        "part": "base",
    }
'''

chainschema = {
        "*pdbcode": [
            ("bptype", {"?*chain": {"?*res": []} } ),
            ("NAprot_hb", {"?*chain": {"?*res": {"?*part": None}} } ),
            ("pdbchains", []),
            ("intraNA_hb", {"?*chain": {"?*res": {"?*part": {"?*pos":None}}} }),
            ("canonized", {"?*chain": [] } ),
            ("mut_sequence", {"*chain": [] }),
            "Nmodels",
            "method",
            ("nachains", []),
            ("mapping", {"*chain": {"*res": None} } ),
            ("interface_hetatoms", {"?*chain": {"?*res": []} }),
            "resolution",
            ("interface_protein", {"*model": {"?*chain": {"?*res": {"?*part": None}}} }),
            ("fragments", {"*model": {"*chain": [[]] }}),
            ("hetnames", {}),
            "sequence",
            ("ss", {"*chain": {"*res": []} }),
            ("missing_atoms", {"?*chain": [] }),
            ("stacking", {"?*chain": {"?*res":{"?*pos":None}} }),
            ("=>", "name")
        ]
}

json.dump(chainschema, open("chainschema.json", "w"), sort_keys = True, indent = 2)
