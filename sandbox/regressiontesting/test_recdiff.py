#!/usr/bin/env python

from unittest import TestCase, main

from recdiff import recdiff, print_recdiff, DiffEqual, DiffMissing

class RecDiffTestCase(TestCase):
    def assertEqual(self, a, b):
        if not (a == b):
            print(a)
            print(b)
            assert a == b
    
    def assertDiffEqual(self, diff):
        self.assertEqual(diff, DiffEqual)

    def test_recdiff_equal_items(self):
        self.assertDiffEqual(recdiff(1, 1))
        self.assertDiffEqual(recdiff(1.1, 1.1+1e-7, epsilon=1e-6))
        self.assertDiffEqual(recdiff(1.1, 1.1-1e-7, epsilon=1e-6))
        self.assertDiffEqual(recdiff("foo", "foo"))

    def test_recdiff_not_equal_items(self):
        self.assertEqual(recdiff(1, 2), (1, 2))
        self.assertEqual(recdiff(1.1, 1.2+1e-7, epsilon=1e-6), (1.1, 1.2+1e-7))
        self.assertEqual(recdiff(1.1, 1.2-1e-7, epsilon=1e-6), (1.1, 1.2-1e-7))
        self.assertEqual(recdiff("foo", "bar"), ("foo", "bar"))

    def test_recdiff_equal_list(self):
        self.assertDiffEqual(recdiff([1, 2], [1, 2]))

    def test_recdiff_not_equal_list(self):
        self.assertEqual(recdiff([1, 2], [1, 3]), [DiffEqual, (2, 3)])

    def test_recdiff_equal_dict(self):
        self.assertDiffEqual(recdiff({1:2}, {1:2}))

    def test_recdiff_not_equal_dict(self):
        self.assertEqual(recdiff({1:2,2:3}, {1:3,3:4}), {1:(2, 3), 2:(3, DiffMissing), 3:(DiffMissing, 4)})

    def test_recdiff_equal_dict_hierarchy(self):
        self.assertDiffEqual(recdiff({1:{2:{3:4,5:6}}}, {1:{2:{3:4,5:6}}}))

    def test_recdiff_not_equal_dict_hierarchy(self):
        self.assertEqual(recdiff({1:{2:{3:4,5:6}}}, {1:{2:{3:4,5:7}}}), {1:{2:{5:(6, 7)}}})

def _test():
    c = RecDiffTestCase()
    for k in dir(c):
        if k.startswith("test_"):
            print("Running ", k)
            getattr(c, k)()

def _example():
    form1 = {
        "num_coefficients": 2,
        "num_arguments": 2,
        "has_default_cell_integral": 1,
        "cell_integrals": { 0: { "tabulate_tensor_input1": ["data"] } },
    }

    form2 = {
        "num_coefficients": 1,
        "num_arguments": 3,
        "has_default_cell_integral": 0,
        "cell_integrals": { 0: { "tabulate_tensor_input1": ["data2"] } },
    }

    diff = recdiff(form1, form2)
    print_recdiff(diff)

if __name__ == "__main__":
    main()

