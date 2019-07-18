 #!/usr/bin/env python3

def is_ss(result):
    ss = "".join([r[0] for r in result]).replace("T", "S").replace("L", "S")
    if ss == "SSS":
        return True
    return False

def is_ds(result):
    ss = "".join([r[0] for r in result]).replace("T", "S")
    if ss == "DDD":
        return True
    return False

def get_mid(result):
    if result[1] is not None:
        return True
    return False

def get_smaller(result, cutoff):
    if result[1] is None:
        return False
    if result[1] < cutoff:
        return True
    return False

def get_mid_bool(result):
    if result[1] is not 0:
        return True
    return False

import functools
def set_contact_parts(cutoff):
    return functools.partial(get_smaller, cutoff = cutoff)

def is_close(result):
    if result[1] <= 99:
        return True
    return False

def is_veryclose(result):
    if result[1] <= 3.5:
        return True
    return False
