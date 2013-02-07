#!/usr/bin/env python
from codegentestcase import CodegenTestCase, unittest

class test_mock_cells(CodegenTestCase):
    '''These tests check that the mock implementations
    of ufc cells are properly constructed. The mock cells
    are not actual subclasses of ufc::cell, instead using
    sufficiently large fixed size arrays to avoid inconvenient
    memory allocation code in the unit tests. The tests here
    can be seen as a demonstration of how to use the mock cells.

    HEADER:
    /**
    Unit tests of mock implementations of ufc cell.
    */
    #include <iostream>
    #include <ufc.h>
    using std::cout;
    using std::endl;

    #include "mock_cells.h"
    using namespace uflacs;
    '''
    def test_mock_interval(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval_1d();
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_EQ(vertex_coordinates[0*mc.geometric_dimension + 0], 0.0);
        ASSERT_EQ(vertex_coordinates[1*mc.geometric_dimension + 0], 2.0);
        """
        code = """
        mc.scale(2.0);
        """
        self.emit_test(code)

    def test_mock_triangle(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_triangle_2d();
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_EQ(vertex_coordinates[0*mc.geometric_dimension + 0], 0.0);
        ASSERT_EQ(vertex_coordinates[0*mc.geometric_dimension + 1], 0.0);
        ASSERT_EQ(vertex_coordinates[1*mc.geometric_dimension + 0], 2.0);
        ASSERT_EQ(vertex_coordinates[1*mc.geometric_dimension + 1], 0.0);
        ASSERT_EQ(vertex_coordinates[2*mc.geometric_dimension + 0], 0.0);
        ASSERT_EQ(vertex_coordinates[2*mc.geometric_dimension + 1], 3.0);
        """
        code = """
        double factors[2] = { 2.0, 3.0 };
        mc.scale(factors);
        """
        self.emit_test(code)

    def test_mock_tetrahedron(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_tetrahedron_3d();
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_EQ(vertex_coordinates[0*mc.geometric_dimension + 0], 2.0);
        ASSERT_EQ(vertex_coordinates[0*mc.geometric_dimension + 1], 3.0);
        ASSERT_EQ(vertex_coordinates[0*mc.geometric_dimension + 2], 4.0);
        ASSERT_EQ(vertex_coordinates[1*mc.geometric_dimension + 0], 3.0);
        ASSERT_EQ(vertex_coordinates[1*mc.geometric_dimension + 1], 3.0);
        ASSERT_EQ(vertex_coordinates[1*mc.geometric_dimension + 2], 4.0);
        ASSERT_EQ(vertex_coordinates[2*mc.geometric_dimension + 0], 2.0);
        ASSERT_EQ(vertex_coordinates[2*mc.geometric_dimension + 1], 4.0);
        ASSERT_EQ(vertex_coordinates[2*mc.geometric_dimension + 2], 4.0);
        ASSERT_EQ(vertex_coordinates[3*mc.geometric_dimension + 0], 2.0);
        ASSERT_EQ(vertex_coordinates[3*mc.geometric_dimension + 1], 3.0);
        ASSERT_EQ(vertex_coordinates[3*mc.geometric_dimension + 2], 5.0);
        """
        code = """
        double offset[3] = { 2.0, 3.0, 4.0 };
        mc.translate(offset);
        """
        self.emit_test(code)

    def test_mock_quadrilateral(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_quadrilateral_2d();
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_EQ(vertex_coordinates[0*mc.geometric_dimension + 0], 2.0 + 2.0*0.0 + 3.0*0.0);
        ASSERT_EQ(vertex_coordinates[0*mc.geometric_dimension + 1], 3.0 + 4.0*0.0 + 5.0*0.0);
        ASSERT_EQ(vertex_coordinates[1*mc.geometric_dimension + 0], 2.0 + 2.0*1.0 + 3.0*0.0);
        ASSERT_EQ(vertex_coordinates[1*mc.geometric_dimension + 1], 3.0 + 4.0*1.0 + 5.0*0.0);
        ASSERT_EQ(vertex_coordinates[2*mc.geometric_dimension + 0], 2.0 + 2.0*1.0 + 3.0*1.0);
        ASSERT_EQ(vertex_coordinates[2*mc.geometric_dimension + 1], 3.0 + 4.0*1.0 + 5.0*1.0);
        ASSERT_EQ(vertex_coordinates[3*mc.geometric_dimension + 0], 2.0 + 2.0*0.0 + 3.0*1.0);
        ASSERT_EQ(vertex_coordinates[3*mc.geometric_dimension + 1], 3.0 + 4.0*0.0 + 5.0*1.0);
        """
        code = """
        double G[4] = { 2.0, 3.0,
                        4.0, 5.0 };
        double x[2] = { 2.0, 3.0 };
        mc.affine_map(G, x);
        """
        self.emit_test(code)

    def test_mock_hexahedron(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_hexahedron_3d();
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_EQ(vertex_coordinates[0*mc.geometric_dimension + 0], 5.0 * (0.0 + 2.0));
        ASSERT_EQ(vertex_coordinates[0*mc.geometric_dimension + 1], 6.0 * (0.0 + 3.0));
        ASSERT_EQ(vertex_coordinates[0*mc.geometric_dimension + 2], 7.0 * (0.0 + 4.0));
        ASSERT_EQ(vertex_coordinates[1*mc.geometric_dimension + 0], 5.0 * (1.0 + 2.0));
        ASSERT_EQ(vertex_coordinates[1*mc.geometric_dimension + 1], 6.0 * (0.0 + 3.0));
        ASSERT_EQ(vertex_coordinates[1*mc.geometric_dimension + 2], 7.0 * (0.0 + 4.0));
        ASSERT_EQ(vertex_coordinates[2*mc.geometric_dimension + 0], 5.0 * (1.0 + 2.0));
        ASSERT_EQ(vertex_coordinates[2*mc.geometric_dimension + 1], 6.0 * (1.0 + 3.0));
        ASSERT_EQ(vertex_coordinates[2*mc.geometric_dimension + 2], 7.0 * (0.0 + 4.0));
        // ...
        ASSERT_EQ(vertex_coordinates[7*mc.geometric_dimension + 0], 5.0 * (0.0 + 2.0));
        ASSERT_EQ(vertex_coordinates[7*mc.geometric_dimension + 1], 6.0 * (1.0 + 3.0));
        ASSERT_EQ(vertex_coordinates[7*mc.geometric_dimension + 2], 7.0 * (1.0 + 4.0));
        """
        code = """
        double offset[3] = { 2.0, 3.0, 4.0 };
        double factors[3] = { 5.0, 6.0, 7.0 };
        mc.translate(offset);
        mc.scale(factors);
        """
        self.emit_test(code)

if __name__ == "__main__":
    unittest.main()
