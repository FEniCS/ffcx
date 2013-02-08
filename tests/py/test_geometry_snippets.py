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


# ...............................................................

dependencies = {
   'J': ('vertex_coordinates',),
   'detJ': ('J',),
   'K': ('J','detJ'),
   'x': ('xi', 'J', 'vertex_coordinates'),
   'xi': ('x', 'K', 'vertex_coordinates'),
  }

# ...............................................................

def generate_jacobian_snippets(cell, restriction):
    decl = 'double J%s[%d*%d];' % (restriction, cell.geometric_dimension(), cell.topological_dimension())
    comp = 'compute_jacobian_%s_%dd(J%s, vertex_coordinates%s);' % (
        cell.cellname(), cell.geometric_dimension(), restriction, restriction)
    return [decl, comp]

def generate_jacobian_inverse_snippets(cell, restriction):
    decl = ['double det%s;' % restriction,
            'double K%s[%d*%d];' % (restriction, cell.geometric_dimension(), cell.topological_dimension())]
    comp = 'compute_jacobian_inverse_%s_%dd(K%s, det%s, J%s);' % (
        cell.cellname(), cell.geometric_dimension(), restriction, restriction, restriction)
    return decl + [comp]

# ...............................................................

def generate_array_definition_snippets(name, expressions, d):
    "Generate combined definition and declaration of name[] = expressions[] dimension d."
    decl = ['double %s[%d] = {' % (name, d)]
    decl += [e+',' for e in expressions[:-1]]
    decl += [expressions[-1]]
    decl += ['    };']
    return decl

def generate_z_Axpy_snippets(name_z, name_A, name_x, name_y, zd, xd):
    fmt_A = { (i,j): '%s[%d*%d+%d]' % (name_A, i, xd, j) for i in xrange(zd) for j in xrange(xd) }
    fmt_x = ['%s[%d]' % (name_x, j) for j in xrange(xd)]
    fmt_y = ['%s[%d]' % (name_y, i) for i in xrange(zd)]
    fmt_Ax = [' + '.join('%s * %s' % (fmt_A[(i,j)], fmt_x[j]) for j in xrange(xd)) for i in xrange(zd)]
    fmt_z = ['%s + %s' % (fmt_Ax[i], fmt_y[i]) for i in xrange(zd)]
    return generate_array_definition_snippets(name_z, fmt_z, zd)

def generate_z_Axmy_snippets(name_z, name_A, name_x, name_y, zd, xd):
    "Generate combined definition and declaration of z = A (x - y) with dimensions zd,xd."
    fmt_A = { (i,j): '%s[%d*%d+%d]' % (name_A, i, xd, j) for i in xrange(zd) for j in xrange(xd) }
    fmt_xmy = ['(%s[%d] - %s[%d])' % (name_x, j, name_y, j) for j in xrange(xd)]
    fmt_z = [' + '.join('%s * %s' % (fmt_A[(i,j)], fmt_xmy[j]) for j in xrange(xd)) for i in xrange(zd)]
    return generate_array_definition_snippets(name_z, fmt_z, zd)

# ...............................................................

def generate_x_from_xi_snippets(cell, restriction):
    "Generate combined definition and declaration of x = J xi + v."
    gd = cell.geometric_dimension()
    td = cell.topological_dimension()
    name_A = "J%s" % restriction
    name_x = "xi%s" % restriction
    name_y = "vertex_coordinates%s" % restriction
    name_z = "x%s" % restriction
    return generate_z_Axpy_snippets(name_z, name_A, name_x, name_y, gd, td)

def generate_xi_from_x_snippets(cell, restriction):
    "Generate combined definition and declaration of xi = K (x - v)."
    gd = cell.geometric_dimension()
    td = cell.topological_dimension()
    name_A = "K%s" % restriction
    name_x = "x%s" % restriction
    name_y = "vertex_coordinates%s" % restriction
    name_z = "xi%s" % restriction
    return generate_z_Axmy_snippets(name_z, name_A, name_x, name_y, gd, td)


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
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
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
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
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
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
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
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
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
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
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
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('tetrahedron', 3)
        snippets = generate_jacobian_snippets(cell, '')
        code = '\n'.join(snippets)
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
        snippets = generate_jacobian_snippets(cell, '')
        snippets += generate_jacobian_inverse_snippets(cell, '')
        code = '\n'.join(snippets)
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
        snippets = generate_jacobian_snippets(cell, '')
        snippets += generate_jacobian_inverse_snippets(cell, '')
        code = '\n'.join(snippets)
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
        snippets = generate_jacobian_snippets(cell, '')
        snippets += generate_jacobian_inverse_snippets(cell, '')
        code = '\n'.join(snippets)
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
        snippets = generate_jacobian_snippets(cell, '')
        snippets += generate_jacobian_inverse_snippets(cell, '')
        code = '\n'.join(snippets)
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
        snippets = generate_jacobian_snippets(cell, '')
        snippets += generate_jacobian_inverse_snippets(cell, '')
        code = '\n'.join(snippets)
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
        snippets = generate_jacobian_snippets(cell, '')
        snippets += generate_jacobian_inverse_snippets(cell, '')
        code = '\n'.join(snippets)
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
        snippets = generate_jacobian_snippets(cell, '')
        snippets += generate_x_from_xi_snippets(cell, '')
        code = '\n'.join(snippets)
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
        snippets = generate_jacobian_snippets(cell, '')
        snippets += generate_x_from_xi_snippets(cell, '')
        code = '\n'.join(snippets)
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
        snippets = generate_jacobian_snippets(cell, '')
        snippets += generate_x_from_xi_snippets(cell, '')
        code = '\n'.join(snippets)
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
        snippets = generate_jacobian_snippets(cell, '')
        snippets += generate_x_from_xi_snippets(cell, '')
        code = '\n'.join(snippets)
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
        cout << endl
             << "J:" << endl
             << J[0] << " " << J[1] << endl
             << J[2] << " " << J[3] << endl
             << J[4] << " " << J[5] << endl
             << "xi:" << endl
             << xi[0] << " " << xi[1] << endl
             << "vc:" << endl
             << vertex_coordinates[0] << " " << vertex_coordinates[1] << " " << vertex_coordinates[2] << endl
             << endl;
        ASSERT_DOUBLE_EQ(0.2*0.3+0.1, x[0]);
        ASSERT_DOUBLE_EQ(0.6*0.3+0.4, x[1]);
        ASSERT_DOUBLE_EQ(0.5, x[3]);
        """
        cell = ufl.Cell('triangle', 3)
        snippets = generate_jacobian_snippets(cell, '')
        snippets += generate_x_from_xi_snippets(cell, '')
        code = '\n'.join(snippets)
        print code
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
        ASSERT_DOUBLE_EQ(0.9*0.3+0.5, x[3]);
        """
        cell = ufl.Cell('tetrahedron', 3)
        snippets = generate_jacobian_snippets(cell, '')
        snippets += generate_x_from_xi_snippets(cell, '')
        code = '\n'.join(snippets)
        self.emit_test(code)

    # ...............................................................

    def test_compute_xi_from_x_interval_1d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval(1);
        double * vertex_coordinates = mc.vertex_coordinates;
        double x[1] = { 0.2 };

        POST:
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('interval', 1)
        snippets = generate_jacobian_snippets(cell, '')
        snippets += generate_jacobian_inverse_snippets(cell, '')
        snippets += generate_xi_from_x_snippets(cell, '')
        code = '\n'.join(snippets)
        self.emit_test(code)

    def test_compute_xi_from_x_interval_2d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval(2);
        double * vertex_coordinates = mc.vertex_coordinates;
        double x[2] = { 0.2, 0.6 };

        POST:
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('interval', 2)
        snippets = generate_jacobian_snippets(cell, '')
        snippets += generate_jacobian_inverse_snippets(cell, '')
        snippets += generate_xi_from_x_snippets(cell, '')
        code = '\n'.join(snippets)
        self.emit_test(code)

    def test_compute_xi_from_x_interval_3d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval(3);
        double * vertex_coordinates = mc.vertex_coordinates;
        double x[3] = { 0.2, 0.6, 0.9 };

        POST:
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('interval', 3)
        snippets = generate_jacobian_snippets(cell, '')
        snippets += generate_jacobian_inverse_snippets(cell, '')
        snippets += generate_xi_from_x_snippets(cell, '')
        code = '\n'.join(snippets)
        self.emit_test(code)

    def test_compute_xi_from_x_triangle_2d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_triangle(2);
        double * vertex_coordinates = mc.vertex_coordinates;
        double x[2] = { 0.2, 0.6 };

        POST:
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('triangle', 2)
        snippets = generate_jacobian_snippets(cell, '')
        snippets += generate_jacobian_inverse_snippets(cell, '')
        snippets += generate_xi_from_x_snippets(cell, '')
        code = '\n'.join(snippets)
        self.emit_test(code)

    def test_compute_xi_from_x_triangle_3d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_triangle(3);
        double * vertex_coordinates = mc.vertex_coordinates;
        double x[3] = { 0.2, 0.6, 0.9 };

        POST:
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('triangle', 3)
        snippets = generate_jacobian_snippets(cell, '')
        snippets += generate_jacobian_inverse_snippets(cell, '')
        snippets += generate_xi_from_x_snippets(cell, '')
        code = '\n'.join(snippets)
        self.emit_test(code)

    def test_compute_xi_from_x_tetrahedron_3d(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_tetrahedron(3);
        double * vertex_coordinates = mc.vertex_coordinates;
        double x[3] = { 0.2, 0.6, 0.9 };

        POST:
        ASSERT_DOUBLE_EQ(1, 1); // FIXME
        """
        cell = ufl.Cell('tetrahedron', 3)
        snippets = generate_jacobian_snippets(cell, '')
        snippets += generate_jacobian_inverse_snippets(cell, '')
        snippets += generate_xi_from_x_snippets(cell, '')
        code = '\n'.join(snippets)
        self.emit_test(code)

    # ...............................................................

    def test_computation_of_geometry_mapping_on_interval(self):
        """
        PRE:
        mock_cell mc;
        mc.fill_reference_interval(1);
        mc.scale(-0.1);
        mc.translate(0.2);

        double * vertex_coordinates = mc.vertex_coordinates;

        POST:
        ASSERT_DOUBLE_EQ(0.2, v0[0]);
        ASSERT_DOUBLE_EQ(-0.1, J[0]);
        ASSERT_DOUBLE_EQ(-0.1, detJ);
        ASSERT_DOUBLE_EQ(-1.0, signdetJ);
        ASSERT_DOUBLE_EQ(0.1, absdetJ);
        ASSERT_DOUBLE_EQ(-1.0/0.1, Jinv[0]);
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
        mc.scale(-0.1);
        mc.translate(0.2);

        double * vertex_coordinatesr = mc.vertex_coordinates;

        POST:
        ASSERT_DOUBLE_EQ(0.2, v0r[0]);
        ASSERT_DOUBLE_EQ(-0.1, Jr[0]);
        ASSERT_DOUBLE_EQ(-0.1, detJr);
        ASSERT_DOUBLE_EQ(-1.0, signdetJr);
        ASSERT_DOUBLE_EQ(0.1, absdetJr);
        ASSERT_DOUBLE_EQ(-1.0/0.1, Jinvr[0]);
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
        mc.scale(-0.1);
        mc.translate(0.2);

        double * vertex_coordinates = mc.vertex_coordinates;
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
        mc.scale(-0.1);
        mc.translate(0.2);

        double * vertex_coordinates = mc.vertex_coordinates;
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
