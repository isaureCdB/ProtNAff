import numpy as np
import json

# "mapping": from numerical string to numerical string
chainschema = {
    "*pdbcode": [
        ("bptype", {"?*chain": {"?*res": []}}),
        ("NAprot_hb", {"?*chain": {"?*res": {"?*part": None}}}),
        ("protchains", []),
        ("intraNA_hb", {"?*chain": {"?*res": {"?*part": {"?*pos": None}}}}),
        ("canonized", {"?*chain": []}),
        "Nmodels",
        "method",
        ("nachains", []),
        ("mapping", {"*chain": {"*res"}}),
        ("interface_hetatoms", {"*model": {"?*chain": {"?*res": []}}}),
        "resolution",
        ("interface_protein", {"*model": {"?*chain": {"?*res": {"?*part": None}}}}),
        ("fragments", {"*model": {"*chain": [[]]}}),
        ("hetnames", {}),
        "sequence",
        ("ss", {"*chain": {"*res": []}}),
        ("missing_atoms", {"?*chain": []}),
        ("stacking", {"?*chain": {"?*res": {"?*pos": None}}}),
        ("=>", "name"),
    ]
}


def query_schema_leaf(schema):
    if schema is None:
        return True
    if schema in ([], [[]], {}):
        return True
    return False


# * = variable instead of string
# ? = optional
def query_schema_star(schema):
    if isinstance(schema, dict):
        assert len(schema.keys()) == 1, schema.keys()
        key = list(schema.keys())[0]
        if key[:2] == "?*":
            return key, True
        else:
            assert key[:1] == "*"
            return key, False
    else:
        return None, False


def query_schema_list(schema):
    if isinstance(schema, list):
        last = schema[-1]
        assert last[0] == "=>", schema
        return last[1]


def _query(data, schema, variables):
    if query_schema_leaf(schema):
        return data
    starkey, optional = query_schema_star(schema)
    if starkey is not None:
        if optional:
            subkey = variables[starkey[2:]]
        else:
            subkey = variables[starkey[1:]]
        try:
            subdata = data[subkey]
        except KeyError:
            if optional:
                return None
            raise KeyError(subkey, starkey)
        subschema = schema[starkey]
        return _query(subdata, subschema, variables)
    listkey = query_schema_list(schema)
    assert listkey is not None, schema
    subkey = variables[listkey]
    subdata = data[subkey]
    for subschema in schema:
        if isinstance(subschema, str):
            if subschema == subkey:
                return subdata
        elif subschema[0] == subkey:
            return _query(subdata, subschema[1], variables)
    else:
        raise Exception("variable not found in schema", subkey, schema)


def query_one_res(chaindata, chainschema, variables):
    """
    Documentation
    """
    return _query(chaindata, chainschema, variables)


def query_one_frag(chaindata, chainschema, frag, name, part=None, pos=None):
    variables = {
        "name": name,
        "model": frag["model"],
        "pdbcode": frag["structure"].decode(),
        "chain": frag["chain"].decode(),
    }
    if part is not None:
        variables["part"] = part
    if pos is not None:
        variables["pos"] = pos
    result = []
    for n in range(3):
        variables.update(
            {
                "res": str(frag["indices"][n]),
            }
        )
        q = query_one_res(chaindata, chainschema, variables)
        result.append(q)

    return result


def query(chaindata, chainschema, frags, func, name, part=None, pos=None):
    result = []
    for frag in frags:
        r = query_one_frag(chaindata, chainschema, frag, name, part)
        rr = func(r)
        result.append(rr)
    return result


'''

"""
Example: Do the single-stranded fragments interact more with the base than with sugar or phosphate
 compared to double-stranded ?
"""

# load fragments from fragments_clust.npy
# load chaindata from structures.json

frag = fragments[0]
variables = {
    "name": "ss",
    "pdbcode": frag["structure"].decode(),
    "chain": frag["chain"].decode(),
    "res": str(frag["indices"][0]),
}
print(query_one_res(chaindata, chainschema, variables))

result = query_one_frag(chaindata, chainschema, frag, "ss")
'''


def is_ss(result):
    ss = "".join([r[0] for r in result]).replace("T", "S")
    if ss == "SSS":
        return True
    return False


def is_ds(result):
    ss = "".join([r[0] for r in result]).replace("T", "S")
    if ss == "DDD":
        return True
    return False


def contact_parts(result):
    if result[1] is not None:
        return True
    return False


"""
ss = ["T", "S", "T"] => True
ss = ["T", "S", "S"] => True
ss [everything else] => False

def ss_true(value):
    if value ==  ["T", "S", "T"]:
        return True
    #...
    return False

make_query("ss", ss_true)
[True, False]

def any_05(value):
    return value[0] < 0.5 or value[1] < 0.5 or value[2] < 0.5

def mean(value):
    return np.mean(value)
make_query("interface_protein", data="ph", mean)
=> [0.6, 0.4, 0.2, ..., N]

# Fragment 2254, "1CVJ", "1", "A", ["5", "6", "7"],
query(name="interface_protein", {"pdbcode": "1CVJ", "model": "1", "chain": "A", "res": "5"}, data="base") => count
query(name="interface_protein", {"pdbcode": "1CVJ", "model": "1", "chain": "A", "res": "6"}, data="base") => count
query(name="interface_protein", {"pdbcode": "1CVJ", "model": "1", "chain": "A", "res": "7"}, data="base") => count

query(name="interface_protein", {"pdbcode": "...", "model": "...", "chain": "...", "res": "..."}, data="base") => count

query(name="ss", {"pdbcode": "...", "model": "...", "chain": "...", "res": "..."}, data="base") => count
d = {
    "name",
    "pdbcode",
    "model",
    "chain",
    "res",
    "data",
    }

apply query to each 3 res of all single stranded fragments => sum up all counts, divide by total number of single stranded
same for double-stranded
"""
