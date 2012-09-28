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

    def test_tabulate_tensor_interval_point(self):
        """Test code generation of body of the ufc function:

        void tabulate_tensor(
            double* A,
            const double * const * w,
            const cell& c,
            const double * x) const;

        PRE:
        double A[1];
        memset(A, 0, sizeof(A));

        double w[1][2] = { { 2.0, 3.0 } };

        mock_interval c;
        c.coordinates[0][0] = 0.1;
        c.coordinates[1][0] = 0.2;

        double x[1] = { 0.15 };

        POST:
        // TODO
        """
        code = "// TODO"
        self.emit_test(code)

    def test_tabulate_tensor_interval_cell(self):
        """Test code generation of body of the ufc function:

        void tabulate_tensor(
            double* A,
            const double * const * w,
            const cell& c) const;

        PRE:
        double A[1];
        memset(A, 0, sizeof(A));
        double w[1][2] = { { 2.0, 3.0 } };
        mock_interval c;
        c.coordinates[0][0] = 0.1;
        c.coordinates[1][0] = 0.2;

        POST:
        // TODO
        """
        code = "// TODO"
        self.emit_test(code)

    def test_tabulate_tensor_interval_facet(self):
        """Test code generation of body of the ufc function:

        void tabulate_tensor(
            double* A,
            const double * const * w,
            const cell& c,
            unsigned int facet) const;

        PRE:
        double A[1];
        memset(A, 0, sizeof(A));

        double w[1][2] = { { 2.0, 3.0 } };

        mock_interval c;
        c.coordinates[0][0] = 0.1;
        c.coordinates[1][0] = 0.2;

        unsigned int facet = 0;

        POST:
        // TODO
        """
        code = "// TODO"
        self.emit_test(code)

    def test_tabulate_tensor_interval_facet(self):
        """Test code generation of body of the ufc function:

        void tabulate_tensor(
            double* A,
            const double * const * w,
            const cell& c0,
            const cell& c1,
            unsigned int facet0,
            unsigned int facet1) const;

        PRE:
        double A[2];
        memset(A, 0, sizeof(A));

        double w[1][4] = { { 2.0, 3.0,
                             4.0, 5.0 } };

        mock_interval c0;
        c0.coordinates[0][0] = 0.1;
        c0.coordinates[1][0] = 0.2;
        mock_interval c1;
        c1.coordinates[0][0] = 0.2;
        c1.coordinates[1][0] = 0.3;

        unsigned int facet0 = 1;
        unsigned int facet1 = 0;

        POST:
        // TODO
        """
        code = "// TODO"
        self.emit_test(code)

if __name__ == "__main__":
    unittest.main()
