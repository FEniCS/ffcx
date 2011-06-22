#!/usr/bin/env python

"""
Tests of generic C++ compilation code.
"""

# These are thin wrappers on top of unittest.TestCase and unittest.main
from ufltestcase import UflTestCase, main
from uflacs.codeutils.cpp_compiler import compile_form
import ufl

class CppFormatterTest(UflTestCase):

    def test_cpp_compilation(self):
        from ufl import cell2D, FiniteElement, Coefficient, dx
        M = Coefficient(FiniteElement("CG", cell2D, 1))**2/2*dx
        code = compile_form(M)
        print code


if __name__ == "__main__":
    main()
