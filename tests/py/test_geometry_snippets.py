#!/usr/bin/env python
from codegentestcase import CodegenTestCase, unittest

import ufl
from uflacs.codeutils.format_code_structure import format_code_structure
from uflacs.geometry.generate_geometry_snippets import (
    generate_array_definition_snippets,
    generate_z_Axpy_snippets,
    generate_z_Axmy_snippets,
    generate_jacobian_snippets,
    generate_jacobian_determinants_snippets,
    generate_jacobian_inverse_snippets,
    generate_x_from_xi_snippets,
    generate_xi_from_x_snippets,
    )
# TODO: Add tests of the rest of the snippets code


class test_geometry_snippets(CodegenTestCase):
    '''Geometry snippets based on ufc_geometry.h.

    TODO: Fill in the blanks!

    HEADER:
    /**
    Unit tests of generated geometry snippet code.
    */
    #include <iostream>
    //#include <ufc.h>
    using std::cout;
    using std::endl;

    // Utility function for debugging values in tests
    /* Copy these into a POST: section for quick display of variable values:
        disp("vc", vertex_coordinates, mc.num_vertices, mc.geometric_dimension);
        disp("J", J, mc.geometric_dimension, mc.topological_dimension);
        disp("K", K, mc.topological_dimension, mc.geometric_dimension);
        disp("x", x, 1, mc.geometric_dimension);
        disp("xi", xi, 1, mc.topological_dimension);
    */
    void disp(const char * name, const double * values, std::size_t m, std::size_t n)
    {
        cout << name << ":" << endl;
        for (std::size_t i = 0; i < m; ++i)
        {
            for (std::size_t j = 0; j < n; ++j)
            {
                cout << values[i*n+j] << " ";
            }
            cout << endl;
        }
    }

    #include "mock_cells.h"
    #include <ufc_geometry.h>
    using namespace uflacs;
    '''

    def test_compute_jacobian_interval_1d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval(1);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('interval', 1)
        r = None
        snippets = generate_jacobian_snippets(cell, r)
        code = format_code_structure(snippets)
        self.emit_test(code)

    def test_compute_jacobian_interval_2d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval(2);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('interval', 2)
        r = None
        snippets = generate_jacobian_snippets(cell, r)
        code = format_code_structure(snippets)
        self.emit_test(code)

    def test_compute_jacobian_interval_3d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval(3);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('interval', 3)
        r = None
        snippets = generate_jacobian_snippets(cell, r)
        code = format_code_structure(snippets)
        self.emit_test(code)

    def test_compute_jacobian_triangle_2d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_triangle(2);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('triangle', 2)
        r = None
        snippets = generate_jacobian_snippets(cell, r)
        code = format_code_structure(snippets)
        self.emit_test(code)

    def test_compute_jacobian_triangle_3d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_triangle(3);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('triangle', 3)
        r = None
        snippets = generate_jacobian_snippets(cell, r)
        code = format_code_structure(snippets)
        self.emit_test(code)

    def test_compute_jacobian_tetrahedron_3d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_tetrahedron(3);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('tetrahedron', 3)
        r = None
        snippets = generate_jacobian_snippets(cell, r)
        code = format_code_structure(snippets)
        self.emit_test(code)

    # ...............................................................

    def test_compute_jacobian_inverse_interval_1d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval(1);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('interval', 1)
        r = None
        snippets = []
        snippets.append( generate_jacobian_snippets(cell, r) )
        snippets.append( generate_jacobian_determinants_snippets(cell, r) )
        snippets.append( generate_jacobian_inverse_snippets(cell, r) )
        code = format_code_structure(snippets)
        self.emit_test(code)

    def test_compute_jacobian_inverse_interval_2d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval(2);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('interval', 2)
        r = None
        snippets = []
        snippets.append( generate_jacobian_snippets(cell, r) )
        snippets.append( generate_jacobian_determinants_snippets(cell, r) )
        snippets.append( generate_jacobian_inverse_snippets(cell, r) )
        code = format_code_structure(snippets)
        self.emit_test(code)

    def test_compute_jacobian_inverse_interval_3d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval(3);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('interval', 3)
        r = None
        snippets = []
        snippets.append( generate_jacobian_snippets(cell, r) )
        snippets.append( generate_jacobian_determinants_snippets(cell, r) )
        snippets.append( generate_jacobian_inverse_snippets(cell, r) )
        code = format_code_structure(snippets)
        self.emit_test(code)

    def test_compute_jacobian_inverse_triangle_2d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_triangle(2);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('triangle', 2)
        r = None
        snippets = []
        snippets.append( generate_jacobian_snippets(cell, r) )
        snippets.append( generate_jacobian_determinants_snippets(cell, r) )
        snippets.append( generate_jacobian_inverse_snippets(cell, r) )
        code = format_code_structure(snippets)
        self.emit_test(code)

    def test_compute_jacobian_inverse_triangle_3d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_triangle(3);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('triangle', 3)
        r = None
        snippets = []
        snippets.append( generate_jacobian_snippets(cell, r) )
        snippets.append( generate_jacobian_determinants_snippets(cell, r) )
        snippets.append( generate_jacobian_inverse_snippets(cell, r) )
        code = format_code_structure(snippets)
        self.emit_test(code)

    def test_compute_jacobian_inverse_tetrahedron_3d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_tetrahedron(3);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('tetrahedron', 3)
        r = None
        snippets = []
        snippets.append( generate_jacobian_snippets(cell, r) )
        snippets.append( generate_jacobian_determinants_snippets(cell, r) )
        snippets.append( generate_jacobian_inverse_snippets(cell, r) )
        code = format_code_structure(snippets)
        self.emit_test(code)

    # ...............................................................

    def test_compute_x_from_xi_interval_1d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval(1);
        mc.scale(0.3); mc.translate(0.1);
        double * vertex_coordinates = mc.vertex_coordinates;
        double xi[1] = { 0.2 };

        POST:
        ASSERT_DOUBLE_EQ(0.2*0.3+0.1, x[0]);
        """
        cell = ufl.Cell('interval', 1)
        r = None
        snippets = []
        snippets.append( generate_jacobian_snippets(cell, r) )
        snippets.append( generate_x_from_xi_snippets(cell, r) )
        code = format_code_structure(snippets)
        self.emit_test(code)

    def test_compute_x_from_xi_interval_2d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval(2);
        mc.scale(0.3); mc.translate(0.1, 0.4);
        double * vertex_coordinates = mc.vertex_coordinates;
        double xi[1] = { 0.2 };

        POST:
        ASSERT_DOUBLE_EQ(0.2*0.3+0.1, x[0]);
        ASSERT_DOUBLE_EQ(0.4, x[1]);
        """
        cell = ufl.Cell('interval', 2)
        r = None
        snippets = []
        snippets.append( generate_jacobian_snippets(cell, r) )
        snippets.append( generate_x_from_xi_snippets(cell, r) )
        code = format_code_structure(snippets)
        self.emit_test(code)

    def test_compute_x_from_xi_interval_3d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval(3);
        mc.scale(0.3); mc.translate(0.1, 0.4, 0.5);

        double * vertex_coordinates = mc.vertex_coordinates;
        double xi[1] = { 0.2 };

        POST:
        ASSERT_DOUBLE_EQ(0.2*0.3+0.1, x[0]);
        ASSERT_DOUBLE_EQ(0.4, x[1]);
        ASSERT_DOUBLE_EQ(0.5, x[2]);
        """
        cell = ufl.Cell('interval', 3)
        r = None
        snippets = []
        snippets.append( generate_jacobian_snippets(cell, r) )
        snippets.append( generate_x_from_xi_snippets(cell, r) )
        code = format_code_structure(snippets)
        self.emit_test(code)

    def test_compute_x_from_xi_triangle_2d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_triangle(2);
        mc.scale(0.3); mc.translate(0.1, 0.4);

        double * vertex_coordinates = mc.vertex_coordinates;
        double xi[2] = { 0.2, 0.6 };

        POST:
        ASSERT_DOUBLE_EQ(0.2*0.3+0.1, x[0]);
        ASSERT_DOUBLE_EQ(0.6*0.3+0.4, x[1]);
        """
        cell = ufl.Cell('triangle', 2)
        r = None
        snippets = []
        snippets.append( generate_jacobian_snippets(cell, r) )
        snippets.append( generate_x_from_xi_snippets(cell, r) )
        code = format_code_structure(snippets)
        self.emit_test(code)

    def test_compute_x_from_xi_triangle_3d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_triangle(3);
        mc.scale(0.3); mc.translate(0.1, 0.4, 0.5);

        double * vertex_coordinates = mc.vertex_coordinates;
        double xi[2] = { 0.2, 0.6 };

        POST:
        ASSERT_DOUBLE_EQ(0.2*0.3+0.1, x[0]);
        ASSERT_DOUBLE_EQ(0.6*0.3+0.4, x[1]);
        ASSERT_DOUBLE_EQ(0.5, x[2]);
        """
        cell = ufl.Cell('triangle', 3)
        r = None
        snippets = []
        snippets.append( generate_jacobian_snippets(cell, r) )
        snippets.append( generate_x_from_xi_snippets(cell, r) )
        code = format_code_structure(snippets)
        self.emit_test(code)

    def test_compute_x_from_xi_tetrahedron_3d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_tetrahedron(3);
        mc.scale(0.3); mc.translate(0.1, 0.4, 0.5);

        double * vertex_coordinates = mc.vertex_coordinates;
        double xi[3] = { 0.2, 0.6, 0.9 };

        POST:
        ASSERT_DOUBLE_EQ(0.2*0.3+0.1, x[0]);
        ASSERT_DOUBLE_EQ(0.6*0.3+0.4, x[1]);
        ASSERT_DOUBLE_EQ(0.9*0.3+0.5, x[2]);
        """
        cell = ufl.Cell('tetrahedron', 3)
        r = None
        snippets = []
        snippets.append( generate_jacobian_snippets(cell, r) )
        snippets.append( generate_x_from_xi_snippets(cell, r) )
        code = format_code_structure(snippets)
        self.emit_test(code)

    # ...............................................................

    def test_compute_xi_from_x_interval_1d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval(1);
        mc.scale(0.3); mc.translate(0.1);

        double * vertex_coordinates = mc.vertex_coordinates;
        double x[1] = { 0.2 };

        POST:
        ASSERT_DOUBLE_EQ(1.0/0.3*(0.2-0.1), xi[0]);
        """
        cell = ufl.Cell('interval', 1)
        r = None
        snippets = []
        snippets.append( generate_jacobian_snippets(cell, r) )
        snippets.append( generate_jacobian_determinants_snippets(cell, r) )
        snippets.append( generate_jacobian_inverse_snippets(cell, r) )
        snippets.append( generate_xi_from_x_snippets(cell, r) )
        code = format_code_structure(snippets)
        self.emit_test(code)

    def test_compute_xi_from_x_interval_2d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval(2);
        mc.scale(0.3); mc.translate(0.1, 0.4);

        double * vertex_coordinates = mc.vertex_coordinates;
        double x[2] = { 0.2, 0.6 };

        POST:
        ASSERT_DOUBLE_EQ(1.0/0.3*(0.2-0.1), xi[0]);
        """
        cell = ufl.Cell('interval', 2)
        r = None
        snippets = []
        snippets.append( generate_jacobian_snippets(cell, r) )
        snippets.append( generate_jacobian_determinants_snippets(cell, r) )
        snippets.append( generate_jacobian_inverse_snippets(cell, r) )
        snippets.append( generate_xi_from_x_snippets(cell, r) )
        code = format_code_structure(snippets)
        self.emit_test(code)

    def test_compute_xi_from_x_interval_3d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval(3);
        mc.scale(0.3); mc.translate(0.1, 0.4, 0.5);

        double * vertex_coordinates = mc.vertex_coordinates;
        double x[3] = { 0.2, 0.6, 0.9 };

        POST:
        ASSERT_DOUBLE_EQ(1.0/0.3*(0.2-0.1), xi[0]);
        """
        cell = ufl.Cell('interval', 3)
        r = None
        snippets = []
        snippets.append( generate_jacobian_snippets(cell, r) )
        snippets.append( generate_jacobian_determinants_snippets(cell, r) )
        snippets.append( generate_jacobian_inverse_snippets(cell, r) )
        snippets.append( generate_xi_from_x_snippets(cell, r) )
        code = format_code_structure(snippets)
        self.emit_test(code)

    def test_compute_xi_from_x_triangle_2d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_triangle(2);
        mc.scale(0.3); mc.translate(0.1, 0.4);

        double * vertex_coordinates = mc.vertex_coordinates;
        double x[2] = { 0.2, 0.6 };

        POST:
        ASSERT_DOUBLE_EQ(1.0/0.3*(0.2-0.1), xi[0]);
        ASSERT_DOUBLE_EQ(1.0/0.3*(0.6-0.4), xi[1]);
        """
        cell = ufl.Cell('triangle', 2)
        r = None
        snippets = []
        snippets.append( generate_jacobian_snippets(cell, r) )
        snippets.append( generate_jacobian_determinants_snippets(cell, r) )
        snippets.append( generate_jacobian_inverse_snippets(cell, r) )
        snippets.append( generate_xi_from_x_snippets(cell, r) )
        code = format_code_structure(snippets)
        self.emit_test(code)

    def test_compute_xi_from_x_triangle_3d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_triangle(3);
        mc.scale(0.3); mc.translate(0.1, 0.4, 0.5);

        double * vertex_coordinates = mc.vertex_coordinates;
        double x[3] = { 0.2, 0.6, 0.9 };

        POST:
        ASSERT_DOUBLE_EQ(1.0/0.3*(0.2-0.1), xi[0]);
        ASSERT_DOUBLE_EQ(1.0/0.3*(0.6-0.4), xi[1]);
        """
        cell = ufl.Cell('triangle', 3)
        r = None
        snippets = []
        snippets.append( generate_jacobian_snippets(cell, r) )
        snippets.append( generate_jacobian_determinants_snippets(cell, r) )
        snippets.append( generate_jacobian_inverse_snippets(cell, r) )
        snippets.append( generate_xi_from_x_snippets(cell, r) )
        code = format_code_structure(snippets)
        self.emit_test(code)

    def test_compute_xi_from_x_tetrahedron_3d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_tetrahedron(3);
        mc.scale(0.3); mc.translate(0.1, 0.4, 0.5);

        double * vertex_coordinates = mc.vertex_coordinates;
        double x[3] = { 0.2, 0.6, 0.9 };

        POST:
        ASSERT_DOUBLE_EQ(1.0/0.3*(0.2-0.1), xi[0]);
        ASSERT_DOUBLE_EQ(1.0/0.3*(0.6-0.4), xi[1]);
        ASSERT_DOUBLE_EQ(1.0/0.3*(0.9-0.5), xi[2]);
        """
        cell = ufl.Cell('tetrahedron', 3)
        r = None
        snippets = []
        snippets.append( generate_jacobian_snippets(cell, r) )
        snippets.append( generate_jacobian_determinants_snippets(cell, r) )
        snippets.append( generate_jacobian_inverse_snippets(cell, r) )
        snippets.append( generate_xi_from_x_snippets(cell, r) )
        code = format_code_structure(snippets)
        self.emit_test(code)

    # ...............................................................

    def xtest_computate_foo_on_bar_Nd(self):
        """
        PRE:
        POST:
        """
        code = ""
        self.emit_test(code)

if __name__ == "__main__":
    unittest.main()
