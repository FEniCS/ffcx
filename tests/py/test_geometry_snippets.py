#!/usr/bin/env python
from codegentestcase import CodegenTestCase, unittest

import ufl
from uflacs.geometry import (
    IntervalGeometryCG,
    TriangleGeometryCG,
    TetrahedronGeometryCG,
    QuadrilateralGeometryCG,
    HexahedronGeometryCG,
    )

def generate_jacobian_snippets(cell, restriction):
    decl = 'double J%s[%d*%d];' % (restriction, cell.geometric_dimension(), cell.topological_dimension())
    comp = 'compute_jacobian_%s_%dd(J%s, vertex_coordinates%s);' % (
        cell.cellname(), cell.geometric_dimension(), restriction, restriction)
    return [decl, comp]

def generate_jacobian_inverse_snippets(cell, restriction):
    decl = ['double det%s;' % restriction,
            'double K%s[%d*%d];' % (restriction, cell.geometric_dimension(), cell.topological_dimension())]
    comp = 'compute_jacobian_inverse_%s_%dd(K%s, det%s, vertex_coordinates%s);' % (
        cell.cellname(), cell.geometric_dimension(), restriction, restriction, restriction)
    return decl + [comp]

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

    #include "mock_cells.h"
    #include "common/ufc_geometry.h"
    using namespace uflacs;
    '''

    def test_compute_jacobian_interval_1d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval(1);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('interval', 1)
        snippets = generate_jacobian_snippets(cell, '')
        code = '\n'.join(snippets)
        self.emit_test(code)

    def test_compute_jacobian_interval_2d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval(2);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('interval', 2)
        snippets = generate_jacobian_snippets(cell, '')
        code = '\n'.join(snippets)
        self.emit_test(code)

    def test_compute_jacobian_interval_3d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval(3);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('interval', 3)
        snippets = generate_jacobian_snippets(cell, '')
        code = '\n'.join(snippets)
        self.emit_test(code)

    def test_compute_jacobian_triangle_2d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_triangle(2);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('triangle', 2)
        snippets = generate_jacobian_snippets(cell, '')
        code = '\n'.join(snippets)
        self.emit_test(code)

    def test_compute_jacobian_triangle_3d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_triangle(3);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('triangle', 3)
        snippets = generate_jacobian_snippets(cell, '')
        code = '\n'.join(snippets)
        self.emit_test(code)

    def test_compute_jacobian_tetrahedron_3d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_tetrahedron(3);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('tetrahedron', 3)
        snippets = generate_jacobian_snippets(cell, '')
        code = '\n'.join(snippets)
        self.emit_test(code)

    def test_compute_jacobian_inverse_interval_1d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval(1);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('interval', 1)
        snippets = generate_jacobian_inverse_snippets(cell, '')
        code = '\n'.join(snippets)
        self.emit_test(code)

    def test_compute_jacobian_inverse_interval_2d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval(2);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('interval', 2)
        snippets = generate_jacobian_inverse_snippets(cell, '')
        code = '\n'.join(snippets)
        self.emit_test(code)

    def test_compute_jacobian_inverse_interval_3d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval(3);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('interval', 3)
        snippets = generate_jacobian_inverse_snippets(cell, '')
        code = '\n'.join(snippets)
        self.emit_test(code)

    def test_compute_jacobian_inverse_triangle_2d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_triangle(2);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('triangle', 2)
        snippets = generate_jacobian_inverse_snippets(cell, '')
        code = '\n'.join(snippets)
        self.emit_test(code)

    def test_compute_jacobian_inverse_triangle_3d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_triangle(3);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('triangle', 3)
        snippets = generate_jacobian_inverse_snippets(cell, '')
        code = '\n'.join(snippets)
        self.emit_test(code)

    def test_compute_jacobian_inverse_tetrahedron_3d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_tetrahedron(3);
        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('tetrahedron', 3)
        snippets = generate_jacobian_inverse_snippets(cell, '')
        code = '\n'.join(snippets)
        self.emit_test(code)


    def test_computation_of_geometry_mapping_on_interval(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval(1);
        double * vertex_coordinates = mc.vertex_coordinates;
        cout << mc.geometric_dimension << endl;
        vertex_coordinates[0*mc.geometric_dimension + 0] = 0.2;
        vertex_coordinates[1*mc.geometric_dimension + 0] = 0.1;

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
        mock_cell mc;
        mc.fill_reference_interval(1);
        double * vertex_coordinatesr = mc.vertex_coordinates;
        vertex_coordinatesr[0*mc.geometric_dimension + 0] = 0.2;
        vertex_coordinatesr[1*mc.geometric_dimension + 0] = 0.1;

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
        mock_cell mc;
        mc.fill_reference_interval(1);
        double * vertex_coordinates = mc.vertex_coordinates;
        vertex_coordinates[0*mc.geometric_dimension + 0] = 0.2;
        vertex_coordinates[1*mc.geometric_dimension + 0] = 0.1;
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
        mock_cell mc;
        mc.fill_reference_interval(1);
        double * vertex_coordinates = mc.vertex_coordinates;
        vertex_coordinates[0*mc.geometric_dimension + 0] = 0.2;
        vertex_coordinates[1*mc.geometric_dimension + 0] = 0.1;
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
