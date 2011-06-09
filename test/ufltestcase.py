#!/usr/bin/env python
import unittest

class UflTestCase(unittest.TestCase):
    def setUp(self):
        super(UflTestCase, self).setUp()
        #print "UflTestCase.setup"

    def tearDown(self):
        #print "UflTestCase.tearDown"
        super(UflTestCase, self).tearDown()

    ### Asserts available in TestCase from python 2.7:

    def _assertIsInstance(self, obj, cl, msg=None):
        self.assertTrue(isinstance(obj, cl), msg=None)

    def _assertNotIsInstance(self, obj, cl, msg=None):
        self.assertFalse(isinstance(obj, cl), msg=msg)

    def _assertIs(self, obj, cl, msg=None):
        self.assertTrue(obj is cl, msg=None)

    def _assertIsNot(self, obj, cl, msg=None):
        self.assertTrue(obj is not cl, msg=msg)

    def _assertIsNone(self, obj, msg=None):
        self.assertTrue(obj is None, msg=msg)

    ### UFL specific asserts

    def assertSameIndices(self, expr, free_indices, msg=None):
        self.assertEqual(expr.free_indices(), free_indices, msg=msg)

    def assertSameShape(self, expr, shape, msg=None):
        self.assertEqual(expr.shape(), shape, msg=msg)

    def assertSameExprProps(self, expr, shape=None, free_indices=None, terminal=None, msg=None):
        if shape is not None:
            self.assertSameShape(expr, shape, msg=msg)
        if free_indices is not None:
            self.assertSameIndices(expr, free_indices, msg=msg)
        if terminal is not None:
            if terminal:
                self.assertIsInstance(expr, Terminal, msg=msg)
            else:
                self.assertIsInstance(expr, Operator, msg=msg)


# Hack for different versions of python unittest:
for func in ('assertIsInstance', 'assertNotIsInstance', 'assertIs', 'assertIsNot', 'assertIsNone'):
    if not hasattr(UflTestCase, func):
        setattr(UflTestCase, func, getattr(UflTestCase, '_'+func))

def main(*args, **kwargs):
    "Hook to do something before running single file tests."
    return unittest.main(*args, **kwargs)

if __name__ == "__main__":
    print "Not to be run directly."
    print "Call main function from this module"
    print "in modules with test cases."

