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
        mock_interval c;

        POST:
        #define cc c.coordinates
        ASSERT_EQ(cc[0][0], 0.0);
        ASSERT_EQ(cc[1][0], 2.0);
        """
        code = """
        scale_cell(c, 2.0);
        """
        self.emit_test(code)

    def test_mock_triangle(self):
        """
        PRE:
        mock_triangle c;

        POST:
        #define cc c.coordinates
        ASSERT_EQ(cc[0][0], 0.0);
        ASSERT_EQ(cc[0][1], 0.0);
        ASSERT_EQ(cc[1][0], 2.0);
        ASSERT_EQ(cc[1][1], 0.0);
        ASSERT_EQ(cc[2][0], 0.0);
        ASSERT_EQ(cc[2][1], 3.0);
        """
        code = """
        double factors[2] = { 2.0, 3.0 };
        scale_cell(c, factors);
        """
        self.emit_test(code)

    def test_mock_tetrahedron(self):
        """
        PRE:
        mock_tetrahedron c;

        POST:
        #define cc c.coordinates
        ASSERT_EQ(cc[0][0], 2.0);
        ASSERT_EQ(cc[0][1], 3.0);
        ASSERT_EQ(cc[0][2], 4.0);
        ASSERT_EQ(cc[1][0], 3.0);
        ASSERT_EQ(cc[1][1], 3.0);
        ASSERT_EQ(cc[1][2], 4.0);
        ASSERT_EQ(cc[2][0], 2.0);
        ASSERT_EQ(cc[2][1], 4.0);
        ASSERT_EQ(cc[2][2], 4.0);
        ASSERT_EQ(cc[3][0], 2.0);
        ASSERT_EQ(cc[3][1], 3.0);
        ASSERT_EQ(cc[3][2], 5.0);
        """
        code = """
        double offset[3] = { 2.0, 3.0, 4.0 };
        translate_cell(c, offset);
        """
        self.emit_test(code)

    def test_mock_quadrilateral(self):
        """
        PRE:
        mock_quadrilateral c;

        POST:
        #define cc c.coordinates
        ASSERT_EQ(cc[0][0], 2.0 + 2.0*0.0 + 3.0*0.0);
        ASSERT_EQ(cc[0][1], 3.0 + 4.0*0.0 + 5.0*0.0);
        ASSERT_EQ(cc[1][0], 2.0 + 2.0*1.0 + 3.0*0.0);
        ASSERT_EQ(cc[1][1], 3.0 + 4.0*1.0 + 5.0*0.0);
        ASSERT_EQ(cc[2][0], 2.0 + 2.0*1.0 + 3.0*1.0);
        ASSERT_EQ(cc[2][1], 3.0 + 4.0*1.0 + 5.0*1.0);
        ASSERT_EQ(cc[3][0], 2.0 + 2.0*0.0 + 3.0*1.0);
        ASSERT_EQ(cc[3][1], 3.0 + 4.0*0.0 + 5.0*1.0);
        """
        code = """
        double G[4] = { 2.0, 3.0,
                        4.0, 5.0 };
        double x[2] = { 2.0, 3.0 };
        map_cell(c, G, x);
        """
        self.emit_test(code)

    def test_mock_hexahedron(self):
        """
        PRE:
        mock_hexahedron c;

        POST:
        #define cc c.coordinates
        ASSERT_EQ(cc[0][0], 5.0 * (0.0 + 2.0));
        ASSERT_EQ(cc[0][1], 6.0 * (0.0 + 3.0));
        ASSERT_EQ(cc[0][2], 7.0 * (0.0 + 4.0));
        ASSERT_EQ(cc[1][0], 5.0 * (1.0 + 2.0));
        ASSERT_EQ(cc[1][1], 6.0 * (0.0 + 3.0));
        ASSERT_EQ(cc[1][2], 7.0 * (0.0 + 4.0));
        ASSERT_EQ(cc[2][0], 5.0 * (1.0 + 2.0));
        ASSERT_EQ(cc[2][1], 6.0 * (1.0 + 3.0));
        ASSERT_EQ(cc[2][2], 7.0 * (0.0 + 4.0));
        // ...
        ASSERT_EQ(cc[7][0], 5.0 * (0.0 + 2.0));
        ASSERT_EQ(cc[7][1], 6.0 * (1.0 + 3.0));
        ASSERT_EQ(cc[7][2], 7.0 * (1.0 + 4.0));
        """
        code = """
        double offset[3] = { 2.0, 3.0, 4.0 };
        double factors[3] = { 5.0, 6.0, 7.0 };
        translate_cell(c, offset);
        scale_cell(c, factors);
        """
        self.emit_test(code)

if __name__ == "__main__":
    unittest.main()
