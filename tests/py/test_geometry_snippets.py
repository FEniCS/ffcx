#!/usr/bin/env python
from codegentestcase import CodegenTestCase, unittest

from uflacs.geometry import (
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

    def test_computation_of_geometry_mapping_on_interval(self):
        """
        PRE:
        mock_interval c;
        c.coordinates[0][0] = 0.2;
        c.coordinates[1][0] = 0.1;

        POST:
        ASSERT_EQ(0.2, v0[0]);
        ASSERT_EQ(-0.1, J[0]);
        ASSERT_EQ(-0.1, detJ);
        ASSERT_EQ(-1.0, signdetJ);
        ASSERT_EQ(0.1, absdetJ);
        ASSERT_EQ(-1.0/0.1, Jinv[0]);
        """
        ccg = IntervalGeometryCG()
        snippets = []
        snippets.append(ccg.v0_code())
        snippets.append(ccg.J_code())
        snippets.append(ccg.detJ_code())
        snippets.append(ccg.signdetJ_code())
        snippets.append(ccg.absdetJ_code())
        snippets.append(ccg.Jinv_code())
        code = '\n'.join(snippets)
        self.emit_test(code)

    def test_computation_of_geometry_mapping_on_restricted_interval(self):
        """
        PRE:
        mock_interval cr;
        cr.coordinates[0][0] = 0.2;
        cr.coordinates[1][0] = 0.1;

        POST:
        ASSERT_EQ(0.2, v0r[0]);
        ASSERT_EQ(-0.1, Jr[0]);
        ASSERT_EQ(-0.1, detJr);
        ASSERT_EQ(-1.0, signdetJr);
        ASSERT_EQ(0.1, absdetJr);
        ASSERT_EQ(-1.0/0.1, Jinvr[0]);
        """
        # Using a custom restriction name 'r' for the test, any string can be inserted instead
        ccg = IntervalGeometryCG(restriction='r')
        snippets = []
        snippets.append(ccg.v0_code())
        snippets.append(ccg.J_code())
        snippets.append(ccg.detJ_code())
        snippets.append(ccg.signdetJ_code())
        snippets.append(ccg.absdetJ_code())
        snippets.append(ccg.Jinv_code())
        code = '\n'.join(snippets)
        self.emit_test(code)

    def test_mapping_from_xi_to_x_on_interval(self):
        """
        PRE:
        mock_interval c;
        c.coordinates[0][0] = 0.2;
        c.coordinates[1][0] = 0.1;
        double xi[1] = { 0.5 };

        POST:
        ASSERT_DOUBLE_EQ(0.15, x[0]);
        """
        ccg = IntervalGeometryCG()
        snippets = []
        snippets.append(ccg.v0_code())
        snippets.append(ccg.J_code())
        snippets.append(ccg.x_from_xi_code())
        code = '\n'.join(snippets)
        self.emit_test(code)

    def test_mapping_from_x_to_xi_on_interval(self):
        """
        PRE:
        mock_interval c;
        c.coordinates[0][0] = 0.2;
        c.coordinates[1][0] = 0.1;
        double x[1] = { 0.15 };

        POST:
        ASSERT_DOUBLE_EQ(0.5, xi[0]);
        """
        ccg = IntervalGeometryCG()
        snippets = []
        snippets.append(ccg.v0_code())
        snippets.append(ccg.J_code())
        snippets.append(ccg.Jinv_code())
        snippets.append(ccg.xi_from_x_code())
        code = '\n'.join(snippets)
        self.emit_test(code)

    def xtest_computation_of_cell_volume_on_interval(self):
        """
        PRE:
        POST:
        """
        code = ""
        self.emit_test(code)

    def xtest_computation_of_cell_surface_area_on_interval(self):
        """
        PRE:
        POST:
        """
        code = ""
        self.emit_test(code)

    def xtest_computation_of_facet_area_on_interval(self):
        """
        PRE:
        POST:
        """
        code = ""
        self.emit_test(code)

    def xtest_computation_of_facet_normal_on_interval(self):
        """
        PRE:
        POST:
        """
        code = ""
        self.emit_test(code)

    def xtest_mapping_from_x_to_xi_on_tetrahedron(self):
        """
        double xi[3];
        xi[0] = G[0]*x[0] + G[1]*x[1] + G[2]*x[2] + v0[0];
        xi[1] = G[3]*x[0] + G[4]*x[1] + G[5]*x[2] + v0[1];
        xi[2] = G[6]*x[0] + G[7]*x[1] + G[8]*x[2] + v0[2];
        """
        code = ""
        self.emit_test(code)

    def xtest_ffc_codesnippets(self):
        """

        PRE:

        POST:

        """
        from ffc.codesnippets import facet_normal
        keys = { 'restriction': '',
                 'direction': '' }
        code = facet_normal[3] % keys
        self.emit_test(code)

if __name__ == "__main__":
    unittest.main()
