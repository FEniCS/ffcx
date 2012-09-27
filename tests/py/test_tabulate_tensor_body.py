#!/usr/bin/env python
from codegentestcase import CodegenTestCase, unittest

import ufl
from ufl import *
#from ufl.common import product

from uflacs.backends.cpp2.compiler import compile_expression

class test_tabulate_tensor_body(CodegenTestCase):
    '''TODO: Fill in the blanks!

    HEADER:
    /**
    Unit tests of tabulate tensor body code structure.
    */
    #include <iostream>
    #include "mock_cells.h"
    using std::cout;
    using std::endl;
    '''

    def compile_expression0(self, expr):
        code = ""
        return code

    def compile_expression(self, expr):
        code = compile_expression(expr, "")
        return code

    def compile_integral(self, integral):
        # TODO: Handle measure type etc.
        if isinstance(integral, Form):
            integral, = integral.integrals()
        expr = integral.integrand()
        code = compile_expression(expr, "")
        return code

    def test_evaluation_of_cell_volume_on_interval(self):
        """Test that...

        PRE:
        mock_interval c;
        double A[1] = { 0.0 };

        POST:
        ASSERT_EQ(A[0], 0.0);
        """
        cell = interval
        x = cell.x
        expr = x
        integral = expr*dP
        code = "" #self.compile_integral(integral)
        self.emit_test(code)

if __name__ == "__main__":
    unittest.main()
