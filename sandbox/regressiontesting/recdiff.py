# -*- coding: utf-8 -*-

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
    keys1 = set(data1.keys())
    keys2 = set(data2.keys())
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
        for k in sorted(diff.keys()):
             key = str(k)
             if prekey: key = ".".join((prekey, key))
             printer("%s%s: " % ("  "*indent, key))
             print_recdiff(diff[k], epsilon, indent+1, printer, key)

    elif isinstance(diff, list):
        for d in diff:
            print_recdiff(d, epsilon, indent+1, printer, prekey)

    elif isinstance(diff, tuple):
        data1, data2 = diff
        printer("%s%s != %s" % ("  "*indent, data1, data2))   

class RecDiffTestCase:#(TestCase):
    def assertEqual(self, a, b):
        if not (a == b):
            print a
            print b
            assert a == b
    
    def assertDiffEqual(self, diff):
        self.assertEqual(diff, DiffEqual)

    def test_recdiff_equal_items(self):
        self.assertDiffEqual(recdiff(1,1))
        self.assertDiffEqual(recdiff(1.1,1.1+1e-7,epsilon=1e-6))
        self.assertDiffEqual(recdiff(1.1,1.1-1e-7,epsilon=1e-6))
        self.assertDiffEqual(recdiff("foo", "foo"))

    def test_recdiff_not_equal_items(self):
        self.assertEqual(recdiff(1,2), (1,2))
        self.assertEqual(recdiff(1.1,1.2+1e-7,epsilon=1e-6), (1.1,1.2+1e-7))
        self.assertEqual(recdiff(1.1,1.2-1e-7,epsilon=1e-6), (1.1,1.2-1e-7))
        self.assertEqual(recdiff("foo", "bar"), ("foo", "bar"))

    def test_recdiff_equal_list(self):
        self.assertDiffEqual(recdiff([1,2], [1,2]))

    def test_recdiff_not_equal_list(self):
        self.assertEqual(recdiff([1,2], [1,3]), [DiffEqual, (2,3)])

    def test_recdiff_equal_dict(self):
        self.assertDiffEqual(recdiff({1:2}, {1:2}))

    def test_recdiff_not_equal_dict(self):
        self.assertEqual(recdiff({1:2,2:3}, {1:3,3:4}), {1:(2,3), 2:(3,DiffMissing), 3:(DiffMissing,4)})

c = RecDiffTestCase()
for k in dir(c):
    if k.startswith("test_"):
        print "Running ", k
        getattr(c,k)()


form1 = {
        "num_coefficients": 2,
        "num_arguments": 2,
        "has_default_cell_integral": 1,
        "cell_integrals": { 0: { "tabulate_tensor_input1": ["data"] } }, }

form2 = {
        "num_coefficients": 1,
        "num_arguments": 3,
        "has_default_cell_integral": 0,
        "cell_integrals": { 0: { "tabulate_tensor_input1": ["data2"] } }, }

diff = recdiff(form1, form2)
print_recdiff(diff)

