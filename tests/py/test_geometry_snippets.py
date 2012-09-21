#!/usr/bin/env python
from codegentestcase import CodegenTestCase, unittest

class test_geometry_snippets(CodegenTestCase):
    '''Geometry snippets based on ufc cell.

    TODO: Fill in the blanks!

    HEADER:
    /**
    Unit tests of generated geometry snippet code.
    */
    #include <iostream>
    #include <ufc.h>
    using std::cout;
    using std::endl;

    #include "mock_cells.h"
    using namespace uflacs;
    '''

    def test_mock_cells(self):
        """
        PRE:
        mock_interval      s1;
        mock_triangle      s2;
        mock_tetrahedron   s3;
        mock_quadrilateral q2;
        mock_hexahedron    q3;

        POST:
        ASSERT_EQ(s1.coordinates[0][0], 0.0);
        ASSERT_EQ(s1.coordinates[1][0], 2.0);

        ASSERT_EQ(s2.coordinates[0][0], 0.0);
        ASSERT_EQ(s2.coordinates[0][1], 0.0);
        ASSERT_EQ(s2.coordinates[1][0], 2.0);
        ASSERT_EQ(s2.coordinates[1][1], 0.0);
        ASSERT_EQ(s2.coordinates[2][0], 0.0);
        ASSERT_EQ(s2.coordinates[2][1], 3.0);

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

        ASSERT_EQ(q2.coordinates[0][0], 2.0 + 2.0*0.0 + 3.0*0.0);
        ASSERT_EQ(q2.coordinates[0][1], 3.0 + 4.0*0.0 + 5.0*0.0);
        ASSERT_EQ(q2.coordinates[1][0], 2.0 + 2.0*1.0 + 3.0*0.0);
        ASSERT_EQ(q2.coordinates[1][1], 3.0 + 4.0*1.0 + 5.0*0.0);
        ASSERT_EQ(q2.coordinates[2][0], 2.0 + 2.0*1.0 + 3.0*1.0);
        ASSERT_EQ(q2.coordinates[2][1], 3.0 + 4.0*1.0 + 5.0*1.0);
        ASSERT_EQ(q2.coordinates[3][0], 2.0 + 2.0*0.0 + 3.0*1.0);
        ASSERT_EQ(q2.coordinates[3][1], 3.0 + 4.0*0.0 + 5.0*1.0);

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
        scale_cell(s1, 2.0);

        double factors2[2] = { 2.0, 3.0 };
        scale_cell(s2, factors2);

        double offset[3] = { 2.0, 3.0, 4.0 };
        translate_cell(s3, offset);

        double G2[4] = { 2.0, 3.0,
                         4.0, 5.0 };
        double x2[2] = { 2.0, 3.0 };
        map_cell(q2, G2, x2);

        double factors3[3] = { 5.0, 6.0, 7.0 };
        translate_cell(q3, offset);
        scale_cell(q3, factors3);
        """
        self.emit_test(code)

    def test_extraction_of_v0(self):
        """
        PRE:
        mock_triangle c;
        c.coordinates[0][0] = 0.1;
        c.coordinates[0][1] = 0.2;

        POST:
        ASSERT_EQ(v0[0], 0.1);
        ASSERT_EQ(v0[1], 0.2);
        """
        code = """
        double v0[2];
        v0[0] = c.coordinates[0][0];
        v0[1] = c.coordinates[0][1];
        """
        self.emit_test(code)

    def test_computation_of_detG(self):
        """
        PRE:
        POST:
        """
        code = ""
        self.emit_test(code)

    def test_mapping_from_x_to_xi(self):
        """Test that foobar.

        PRE:
        double x[3] = { 0.0, 0.1, 0.2 };
        double vertices[3][3] = { { 1.0, 2.0, 3.0 },
                                  { 2.0, 2.0, 2.0 },
                                  { 2.0, 2.0, 2.0 } };

        POST:
        ASSERT_EQ(xi[0], 1.0+0.0);
        ASSERT_EQ(xi[1], 2.0+0.1);
        ASSERT_EQ(xi[2], 3.0+0.2);
        """

        code = """
        double *v0 = &vertices[0][0];

        double G[9];
        G[0] = 1.0;
        G[1] = 0.0;
        G[2] = 0.0;
        G[3] = 0.0;
        G[4] = 1.0;
        G[5] = 0.0;
        G[6] = 0.0;
        G[7] = 0.0;
        G[8] = 1.0;

        double xi[3];
        xi[0] = G[0]*x[0] + G[1]*x[1] + G[2]*x[2] + v0[0];
        xi[1] = G[3]*x[0] + G[4]*x[1] + G[5]*x[2] + v0[1];
        xi[2] = G[6]*x[0] + G[7]*x[1] + G[8]*x[2] + v0[2];
        """

        self.emit_test(code)

    def test_mapping_from_xi_to_x(self):
        """
        PRE:
        POST:
        """
        code = ""
        self.emit_test(code)

    def test_computation_of_cell_volume(self):
        """
        PRE:
        POST:
        """
        code = ""
        self.emit_test(code)

    def test_computation_of_cell_surface_area(self):
        """
        PRE:
        POST:
        """
        code = ""
        self.emit_test(code)

    def test_computation_of_facet_area(self):
        """
        PRE:
        POST:
        """
        code = ""
        self.emit_test(code)

    def test_computation_of_facet_normal(self):
        """
        PRE:
        POST:
        """
        code = ""
        self.emit_test(code)

if __name__ == "__main__":
    unittest.main()
