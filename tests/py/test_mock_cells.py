#!/usr/bin/env python
from codegentestcase import CodegenTestCase, unittest

class test_mock_cells(CodegenTestCase):
    '''Test mock implementations of ufc cell.

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
        mock_interval      s1;

        POST:
        ASSERT_EQ(s1.coordinates[0][0], 0.0);
        ASSERT_EQ(s1.coordinates[1][0], 2.0);
        """
        code = """
        scale_cell(s1, 2.0);
        """
        self.emit_test(code)

    def test_mock_triangle(self):
        """
        PRE:
        mock_triangle      s2;

        POST:
        ASSERT_EQ(s2.coordinates[0][0], 0.0);
        ASSERT_EQ(s2.coordinates[0][1], 0.0);
        ASSERT_EQ(s2.coordinates[1][0], 2.0);
        ASSERT_EQ(s2.coordinates[1][1], 0.0);
        ASSERT_EQ(s2.coordinates[2][0], 0.0);
        ASSERT_EQ(s2.coordinates[2][1], 3.0);
        """
        code = """
        double factors2[2] = { 2.0, 3.0 };
        scale_cell(s2, factors2);
        """
        self.emit_test(code)

    def test_mock_tetrahedron(self):
        """
        PRE:
        mock_tetrahedron   s3;

        POST:
        ASSERT_EQ(s3.coordinates[0][0], 2.0);
        ASSERT_EQ(s3.coordinates[0][1], 3.0);
        ASSERT_EQ(s3.coordinates[0][2], 4.0);
        ASSERT_EQ(s3.coordinates[1][0], 3.0);
        ASSERT_EQ(s3.coordinates[1][1], 3.0);
        ASSERT_EQ(s3.coordinates[1][2], 4.0);
        ASSERT_EQ(s3.coordinates[2][0], 2.0);
        ASSERT_EQ(s3.coordinates[2][1], 4.0);
        ASSERT_EQ(s3.coordinates[2][2], 4.0);
        ASSERT_EQ(s3.coordinates[3][0], 2.0);
        ASSERT_EQ(s3.coordinates[3][1], 3.0);
        ASSERT_EQ(s3.coordinates[3][2], 5.0);
        """
        code = """
        double offset[3] = { 2.0, 3.0, 4.0 };
        translate_cell(s3, offset);
        """
        self.emit_test(code)

    def test_mock_quadrilateral(self):
        """
        PRE:
        mock_quadrilateral q2;

        POST:
        ASSERT_EQ(q2.coordinates[0][0], 2.0 + 2.0*0.0 + 3.0*0.0);
        ASSERT_EQ(q2.coordinates[0][1], 3.0 + 4.0*0.0 + 5.0*0.0);
        ASSERT_EQ(q2.coordinates[1][0], 2.0 + 2.0*1.0 + 3.0*0.0);
        ASSERT_EQ(q2.coordinates[1][1], 3.0 + 4.0*1.0 + 5.0*0.0);
        ASSERT_EQ(q2.coordinates[2][0], 2.0 + 2.0*1.0 + 3.0*1.0);
        ASSERT_EQ(q2.coordinates[2][1], 3.0 + 4.0*1.0 + 5.0*1.0);
        ASSERT_EQ(q2.coordinates[3][0], 2.0 + 2.0*0.0 + 3.0*1.0);
        ASSERT_EQ(q2.coordinates[3][1], 3.0 + 4.0*0.0 + 5.0*1.0);
        """
        code = """
        double G[4] = { 2.0, 3.0,
                        4.0, 5.0 };
        double x[2] = { 2.0, 3.0 };
        map_cell(q2, G, x);
        """
        self.emit_test(code)

    def test_mock_hexahedron(self):
        """
        PRE:
        mock_hexahedron    q3;

        POST:
        ASSERT_EQ(q3.coordinates[0][0], 5.0 * (0.0 + 2.0));
        ASSERT_EQ(q3.coordinates[0][1], 6.0 * (0.0 + 3.0));
        ASSERT_EQ(q3.coordinates[0][2], 7.0 * (0.0 + 4.0));
        ASSERT_EQ(q3.coordinates[1][0], 5.0 * (1.0 + 2.0));
        ASSERT_EQ(q3.coordinates[1][1], 6.0 * (0.0 + 3.0));
        ASSERT_EQ(q3.coordinates[1][2], 7.0 * (0.0 + 4.0));
        ASSERT_EQ(q3.coordinates[2][0], 5.0 * (1.0 + 2.0));
        ASSERT_EQ(q3.coordinates[2][1], 6.0 * (1.0 + 3.0));
        ASSERT_EQ(q3.coordinates[2][2], 7.0 * (0.0 + 4.0));
        // ...
        ASSERT_EQ(q3.coordinates[7][0], 5.0 * (0.0 + 2.0));
        ASSERT_EQ(q3.coordinates[7][1], 6.0 * (1.0 + 3.0));
        ASSERT_EQ(q3.coordinates[7][2], 7.0 * (1.0 + 4.0));
        """
        code = """
        double offset[3] = { 2.0, 3.0, 4.0 };
        double factors[3] = { 5.0, 6.0, 7.0 };
        translate_cell(q3, offset);
        scale_cell(q3, factors);
        """
        self.emit_test(code)

if __name__ == "__main__":
    unittest.main()
