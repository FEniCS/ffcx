#!/usr/bin/env python
from codegentestcase import CodegenTestCase, unittest

import ufl
from ufl import *
#from ufl.common import product

from uflacs.backends.cpp2.compiler import compile_expression
""" FIXME: Update to these:
    /// UFC tabulate_tensor signatures:

    /// Tabulate the tensor for the contribution from a local cell
    virtual void tabulate_tensor(double* A,
                                 const double * const * w,
                                 const double* vertex_coordinates) const = 0;

    /// Tabulate the tensor for the contribution from a local exterior facet
    virtual void tabulate_tensor(double* A,
                                 const double * const * w,
                                 const double* vertex_coordinates,
                                 std::size_t facet) const = 0;

    /// Tabulate the tensor for the contribution from a local interior facet
    virtual void tabulate_tensor(double* A,
                                 const double * const * w,
                                 const double* vertex_coordinates_0,
                                 const double* vertex_coordinates_1,
                                 std::size_t facet0,
                                 std::size_t facet1) const = 0;

    /// Tabulate the tensor for the contribution from the local vertex
    virtual void tabulate_tensor(double* A,
                                 const double * const * w,
                                 const double* vertex_coordinates,
                                 std::size_t vertex) const = 0;

    /// Tabulate the tensor for the contributions from a set of points
    virtual void tabulate_tensor(double* A,
                                 const double * const * w,
                                 const double* vertex_coordinates,
                                 std::size_t num_points,
                                 const double * const * points) const = 0;
*/
"""
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

    def compile_tabulate_tensor_body(self, integral):
        # TODO: Handle measure type etc.
        if isinstance(integral, Form):
            integral, = integral.integrals()
        expr = integral.integrand()
        code = compile_expression(expr, "")
        return code

    def test_interval_tabten_x_given(self):
        """Test code generation of body of the ufc function:

        void tabulate_tensor(
            double* A,
            const double * const * w,
            const double * vertex_coordinates) const;

        PRE:
        mock_cell mc;
        mc.fill_reference_interval_1d();
        double * vertex_coordinates = mc.vertex_coordinates;
        double A[1] = { 0.0 };
        double x[1] = { 1.2 };

        POST:
        ASSERT_EQ(A[0], 1.2);
        """
        cell = interval

        x = cell.x[0]
        expr = x
        integral = expr*dP

        code = self.compile_tabulate_tensor_body(integral)
        self.emit_test(code)

    def test_interval_tabten_dg0_given(self):
        """Test code generation of body of the ufc function:

        void tabulate_tensor(
            double* A,
            const double * const * w,
            const cell& c,
            const double * x) const;

        with mock values for DG0/Real coefficients in w[][].

        PRE:
        double A[1];
        memset(A, 0, sizeof(A));

        double w[2][2] = { { 1.2, 0.0 }, // second value here is unused
                           { 2.0, 3.0 } };

        mock_cell mc;
        mc.fill_reference_interval_1d();
        double * vertex_coordinates = mc.vertex_coordinates;
        vertex_coordinates[0*mc.geometric_dimension + 0] = 0.1;
        vertex_coordinates[1*mc.geometric_dimension + 0] = 0.2;

        double x[1] = { 0.15 };

        POST:
        ASSERT_EQ(A[0], 0.15*1.2*(2.0+3.0));
        """
        cell = interval
        x = cell.x[0]

        V = VectorElement("DG", cell, 0, dim=2)
        w0 = Constant(cell, count=0)
        w1 = Coefficient(V, count=1)

        expr = x*w0*(w1[0] + w1[1])

        integral = expr*dP
        code = self.compile_tabulate_tensor_body(integral)
        self.emit_test(code)

    def xtest_interval_tabten_geometry_mappings(self):
        """Test code generation of body of the ufc function:

        void tabulate_tensor(
            double* A,
            const double * const * w,
            const cell& c) const;

        PRE:
        double A[1];
        memset(A, 0, sizeof(A));

        mock_cell mc;
        mc.fill_reference_interval_1d();
        double * vertex_coordinates = mc.vertex_coordinates;
        vertex_coordinates[0*mc.geometric_dimension + 0] = 0.2;
        vertex_coordinates[1*mc.geometric_dimension + 0] = 0.1;

        POST:
        ASSERT_EQ(J[0], -0.1);
        ASSERT_EQ(Jinv[0], -10.0);
        ASSERT_EQ(detJ, -0.1);
        ASSERT_EQ(volume, 0.1);
        ASSERT_EQ(D, 0.1);

        ASSERT_EQ(A[0], 0.03*2);
        """
        cell = interval
        x = cell.x[0]
        J = cell.J[0,0]
        Jinv = cell.Jinv[0,0]
        detJ = cell.detJ
        xi = cell.xi[0]
        T = cell.volume

        expr = 3*T

        integral = expr*dx
        code = self.compile_tabulate_tensor_body(integral)
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

        mock_cell mc;
        mc.fill_reference_interval_1d();
        double * vertex_coordinates = mc.vertex_coordinates;
        vertex_coordinates[0*mc.geometric_dimension + 0] = 0.1;
        vertex_coordinates[1*mc.geometric_dimension + 0] = 0.2;

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

        mock_cell mc0;
        mock_cell mc1;
        mc0.fill_reference_interval_1d();
        mc1.fill_reference_interval_1d();
        double * vertex_coordinates0 = mc0.vertex_coordinates;
        vertex_coordinates0[0*mc0.geometric_dimension + 0] = 0.1;
        vertex_coordinates0[1*mc0.geometric_dimension + 0] = 0.2;
        double * vertex_coordinates1 = mc1.vertex_coordinates;
        vertex_coordinates1[0*mc1.geometric_dimension + 0] = 0.2;
        vertex_coordinates1[1*mc1.geometric_dimension + 0] = 0.3;

        unsigned int facet0 = 1;
        unsigned int facet1 = 0;

        POST:
        // TODO
        """
        code = "// TODO"
        self.emit_test(code)

if __name__ == "__main__":
    unittest.main()
