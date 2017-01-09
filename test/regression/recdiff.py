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

_default_recdiff_tolerance = 1e-5


def recdiff_dict(data1, data2, tolerance=_default_recdiff_tolerance):
    keys1 = set(data1.keys())
    keys2 = set(data2.keys())
    keys = keys1.intersection(keys2)
    diff = {}
    for k in keys1 - keys:
        diff[k] = (data1[k], DiffMissing)
    for k in keys2 - keys:
        diff[k] = (DiffMissing, data2[k])
    for k in keys:
        d1 = data1[k]
        d2 = data2[k]
        d = recdiff(d1, d2, tolerance)
        if d is not DiffEqual:
            diff[k] = d
    return diff or DiffEqual


def recdiff(data1, data2, tolerance=_default_recdiff_tolerance):
    if isinstance(data1, (float, int)) and isinstance(data2, (float, int)):
        # This approach allows numbers formatted as ints and floats interchangably as long as the values are equal
        delta = abs(data1 - data2)
        avg = (abs(data1) + abs(data2)) / 2.0

        if 0:
            # Using relative comparison, i.e. a tolerance of 1e-2 means one percent error is acceptable
            eps = tolerance * avg
            same = avg < 1e-14 or delta < eps
        else:
            # Using absolute comparison, this is what the old .out comparison does
            same = delta < tolerance

        return DiffEqual if same else (data1, data2)
    elif not isinstance(data1, type(data2)):
        return (data1, data2)
    elif isinstance(data1, dict):
        return recdiff_dict(data1, data2, tolerance)
    elif isinstance(data1, list):
        diff = [recdiff(d1, d2, tolerance) for (d1, d2) in zip(data1, data2)]
        return DiffEqual if all(d is DiffEqual for d in diff) else diff
    else:
        return DiffEqual if data1 == data2 else (data1, data2)


def _print(line):
    print(line)


def print_recdiff(diff, indent=0, printer=_print, prekey=""):

    if isinstance(diff, dict):
        for k in sorted(diff.keys()):
            key = str(k)
            if prekey:
                key = ".".join((prekey, key))
            printer("%s%s: " % ("  " * indent, key))
            print_recdiff(diff[k], indent + 1, printer, key)

    elif isinstance(diff, list):
        # Limiting this to lists of scalar values!
        for i, d in enumerate(diff):
            if isinstance(d, tuple):
                data1, data2 = d
                printer("%s%d: %s != %s" % ("  " * indent, i, data1, data2))

    elif isinstance(diff, tuple):
        assert len(diff) == 2
        data1, data2 = diff
        data1 = str(data1)
        data2 = str(data2)
        if len(data1) + len(data2) + 2 * indent + 4 > 70:
            printer("%s%s" % ("  " * indent, data1))
            printer("%s!=" % ("  " * indent))
            printer("%s%s" % ("  " * indent, data2))
        else:
            printer("%s%s != %s" % ("  " * indent, data1, data2))


# ---------- Unittest code
import unittest
# from recdiff import recdiff, print_recdiff, DiffEqual, DiffMissing


class RecDiffTestCase(unittest.TestCase):

    def assertEqual(self, a, b):
        if not (a == b):
            print(a)
            print(b)
            assert a == b

    def assertDiffEqual(self, diff):
        self.assertEqual(diff, DiffEqual)

    def test_recdiff_equal_items(self):
        self.assertDiffEqual(recdiff(1, 1))
        self.assertDiffEqual(recdiff(0, 0))
        self.assertDiffEqual(recdiff(0, 1e-15))
        self.assertDiffEqual(recdiff(1.1, 1.1 + 1e-7, tolerance=1e-6))
        self.assertDiffEqual(recdiff(1.1, 1.1 - 1e-7, tolerance=1e-6))
        self.assertDiffEqual(recdiff("foo", "foo"))

    def test_recdiff_not_equal_items(self):
        self.assertEqual(recdiff(1, 2), (1, 2))
        self.assertEqual(recdiff(0, 0.0001), (0, 0.0001))
        self.assertEqual(recdiff(0, 1e-13), (0, 1e-13))
        self.assertEqual(recdiff(1.1, 1.2 + 1e-7, tolerance=1e-6), (1.1, 1.2 + 1e-7))
        self.assertEqual(recdiff(1.1, 1.2 - 1e-7, tolerance=1e-6), (1.1, 1.2 - 1e-7))
        self.assertEqual(recdiff("foo", "bar"), ("foo", "bar"))

    def test_recdiff_equal_list(self):
        self.assertDiffEqual(recdiff([1, 2], [1, 2]))

    def test_recdiff_not_equal_list(self):
        self.assertEqual(recdiff([1, 2], [1, 3]), [DiffEqual, (2, 3)])

    def test_recdiff_equal_dict(self):
        self.assertDiffEqual(recdiff({1: 2}, {1: 2}))

    def test_recdiff_not_equal_dict(self):
        self.assertEqual(recdiff({1: 2, 2: 3}, {1: 3, 3: 4}), {1: (2, 3), 2: (3, DiffMissing), 3: (DiffMissing, 4)})

    def test_recdiff_equal_dict_hierarchy(self):
        self.assertDiffEqual(recdiff({1: {2: {3: 4, 5: 6}}}, {1: {2: {3: 4, 5: 6}}}))

    def test_recdiff_not_equal_dict_hierarchy(self):
        self.assertEqual(recdiff({1: {2: {3: 4, 5: 6}}}, {1: {2: {3: 4, 5: 7}}}), {1: {2: {5: (6, 7)}}})

    def test_example(self):
        form1 = {
            "num_coefficients": 2,
            "num_arguments": 2,
            "has_default_cell_integral": 1,
            "cell_integrals": {0: {"tabulate_tensor_input1": ["data"]}},
        }

        form2 = eval("""{
            "num_coefficients": 2,
            "rank": 2,
            "has_default_cell_integral": 0,
            "cell_integrals": { 0: { "tabulate_tensor_input1": ["data2"] } },
        }""")

        actual_diff = recdiff(form1, form2)
        if 0:
            print_recdiff(actual_diff)

        expected_diff = {
            #"num_coefficients": DiffEqual,
            "num_arguments": (2, DiffMissing),
            "rank": (DiffMissing, 2),
            "has_default_cell_integral": (1, 0),
            "cell_integrals": {0: {"tabulate_tensor_input1": [("data", "data2")]}},
        }
        self.assertEqual(actual_diff, expected_diff)


def main(a, b, tolerance=_default_recdiff_tolerance):
    print("Running diff on files %s and %s" % (a, b))
    a = eval(open(a).read())
    b = eval(open(b).read())
    d = recdiff(a, b, float(tolerance))
    print_recdiff(d)


if __name__ == "__main__":
    import sys
    args = sys.argv[1:]
    if not args:  # Hack to be able to use this as a script, TODO: do something nicer
        print("No arguments, running tests.")
        unittest.main()
    else:
        main(*args)
