#!/usr/bin/env python
from codegentestcase import CodegenTestCase, unittest

from uflacs.geometry.cellcodegen import (
    IntervalGeometryCG,
    TriangleGeometryCG,
    TetrahedronGeometryCG,
    QuadrilateralGeometryCG,
    HexahedronGeometryCG,
    )

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

    def test_interval_computation_of_geometry_mapping(self):
        """
        PRE:
        mock_interval c;
        c.coordinates[0][0] = 0.2;
        c.coordinates[1][0] = 0.1;

        POST:
        ASSERT_EQ(v0[0], 0.2);
        ASSERT_EQ(G[0], -0.1);
        ASSERT_EQ(detG, -0.1);
        ASSERT_EQ(detGsign, -1.0);
        ASSERT_EQ(absdetG, 0.1);
        ASSERT_EQ(Ginv[0], -1.0/0.1);
        """
        ccg = IntervalGeometryCG()
        snippets = []
        snippets.append(ccg.gen_v0_code())
        snippets.append(ccg.gen_G_code())
        snippets.append(ccg.gen_detG_code())
        snippets.append(ccg.gen_absdetG_code())
        snippets.append(ccg.gen_Ginv_code())
        code = '\n'.join(snippets)
        self.emit_test(code)

    def test_interval_mapping_between_x_and_xi(self):
        pass

    def test_interval_mapping_from_x_to_xi(self):
        pass

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
