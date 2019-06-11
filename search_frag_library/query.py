#!/usr/bin/env python3
import numpy as np
import json
import sys

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
        assert len(list(schema.keys())) == 1, list(schema.keys())
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
            subkey = variables[starkey[2:]] #select "res" in "?*res"
        else:
            print(starkey[1:])
            print(variables)
            subkey = variables[starkey[1:]] #select "res" in "*res"
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

def query_one_frag(chaindata, chainschema, frag, default, name, result, part=None, pos=None):
    variables = {
        "name": name,
        "model": "model_" + str(frag["model"]),
        "pdbcode": frag["structure"].decode(),
        "chain": "chain_" + frag["chain"].decode(),
    }
    if part is not None:
        variables["part"] = part
    if pos is not None:
        variables["pos"] = pos
    for n in range(3):
        variables.update({
            "res": "res_" + str(frag["indices"][n]),
        })
        q = query_one_res(chaindata, chainschema, variables)
        if q is None:
            q = default
        try:
            result[n] = q
        except Exception as exc:
            msg = "ERROR, query result: %s, expected type: %s\n"
            raise ValueError(msg % (q, result.dtype))


    return result

def query(chaindata, chainschema, frags, default, name, part=None, pos=None):
    default_type = default if isinstance(default, type) else type(default)
    result = np.zeros((len(frags),3), dtype=default_type)
    for nr, frag in enumerate(frags):
        r = query_one_frag(chaindata, chainschema, frag, default, name, result[nr], part, pos)
    return result
