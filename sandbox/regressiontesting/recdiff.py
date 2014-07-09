# -*- coding: utf-8 -*-

from six.moves import zip, iterkeys

class DiffMarkerType:
    def __init__(self, name):
        self.name = name

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

DiffMissing = DiffMarkerType("<value missing>")
DiffEqual = DiffMarkerType("<equal>")

_default_recdiff_epsilon = 1e-8

def recdiff_dict(data1, data2, epsilon=_default_recdiff_epsilon):
    keys1 = set(iterkeys(data1))
    keys2 = set(iterkeys(data2))
    keys = keys1.intersection(keys2)
    diff = {}
    for k in keys1-keys:
        diff[k] = (data1[k], DiffMissing)
    for k in keys2-keys:
        diff[k] = (DiffMissing, data2[k])
    for k in keys:
        d1 = data1[k]
        d2 = data2[k]
        d = recdiff(d1, d2, epsilon)
        if d is not DiffEqual:
            diff[k] = d
    return diff or DiffEqual

def recdiff(data1, data2, epsilon=_default_recdiff_epsilon):
    if type(data1) != type(data2):
        return (data1, data2)
    elif isinstance(data1, dict):
        return recdiff_dict(data1, data2, epsilon)
    elif isinstance(data1, list):
        diff = [recdiff(d1, d2, epsilon) for (d1,d2) in zip(data1, data2)]
        return DiffEqual if all(d is DiffEqual for d in diff) else diff
    elif isinstance(data1, float):
        diff = data1 - data2
        return DiffEqual if abs(diff) < epsilon else (data1, data2)
    else:
        return DiffEqual if data1 == data2 else (data1, data2)

def _print(line):
    print line

def print_recdiff(diff, epsilon=_default_recdiff_epsilon, indent=0, printer=_print, prekey=""):

    if isinstance(diff, dict):
        for k in sorted(iterkeys(diff)):
             key = str(k)
             if prekey: key = ".".join((prekey, key))
             printer("%s%s: " % ("  "*indent, key))
             print_recdiff(diff[k], epsilon, indent+1, printer, key)

    elif isinstance(diff, list):
        # Limiting this to lists of scalar values!
        for i, d in enumerate(diff):
            if isinstance(d, tuple):
                data1, data2 = d
                printer("%s%d: %s != %s" % ("  "*indent, i, data1, data2))

    elif isinstance(diff, tuple):
        assert len(diff) == 2
        data1, data2 = diff
        data1 = str(data1)
        data2 = str(data2)
        if len(data1) + len(data2) + 2*indent + 4 > 70:
            printer("%s%s" % ("  "*indent, data1))
            printer("%s!=" % ("  "*indent))
            printer("%s%s" % ("  "*indent, data2))
        else:
            printer("%s%s != %s" % ("  "*indent, data1, data2))
