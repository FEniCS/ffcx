#!/usr/bin/env py.test

"""
Unit tests of generated geometry snippet code.
"""

import ufl
from uflacs.codeutils.format_code import format_code
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

debugging = '''
/* Copy these into a """

    post = """ section for quick display of variable values:
    disp("vc", vertex_coordinates, mc.num_vertices, mc.geometric_dimension);
    disp("J", J, mc.geometric_dimension, mc.topological_dimension);
    disp("K", K, mc.topological_dimension, mc.geometric_dimension);
    disp("x", x, 1, mc.geometric_dimension);
    disp("xi", xi, 1, mc.topological_dimension);
*/
'''

def test_compute_jacobian_interval_1d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_interval(1);
    double * vertex_coordinates = mc.vertex_coordinates;
    """

    post = """
    ASSERT_DOUBLE_EQ(1, 1); // FIXME
    """

    cell = ufl.Cell('interval', 1)
    r = None
    snippets = generate_jacobian_snippets(cell, r)
    code = format_code(snippets)

    gtest.add(pre + code + post)

def test_compute_jacobian_interval_2d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_interval(2);
    double * vertex_coordinates = mc.vertex_coordinates;
    """

    post = """
    ASSERT_DOUBLE_EQ(1, 1); // FIXME
    """

    cell = ufl.Cell('interval', 2)
    r = None
    snippets = generate_jacobian_snippets(cell, r)
    code = format_code(snippets)

    gtest.add(pre + code + post)

def test_compute_jacobian_interval_3d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_interval(3);
    double * vertex_coordinates = mc.vertex_coordinates;
    """

    post = """
    ASSERT_DOUBLE_EQ(1, 1); // FIXME
    """

    cell = ufl.Cell('interval', 3)
    r = None
    snippets = generate_jacobian_snippets(cell, r)
    code = format_code(snippets)

    gtest.add(pre + code + post)

def test_compute_jacobian_triangle_2d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_triangle(2);
    double * vertex_coordinates = mc.vertex_coordinates;
    """

    post = """
    ASSERT_DOUBLE_EQ(1, 1); // FIXME
    """

    cell = ufl.Cell('triangle', 2)
    r = None
    snippets = generate_jacobian_snippets(cell, r)
    code = format_code(snippets)

    gtest.add(pre + code + post)

def test_compute_jacobian_triangle_3d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_triangle(3);
    double * vertex_coordinates = mc.vertex_coordinates;
    """

    post = """
    ASSERT_DOUBLE_EQ(1, 1); // FIXME
    """

    cell = ufl.Cell('triangle', 3)
    r = None
    snippets = generate_jacobian_snippets(cell, r)
    code = format_code(snippets)

    gtest.add(pre + code + post)

def test_compute_jacobian_tetrahedron_3d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_tetrahedron(3);
    double * vertex_coordinates = mc.vertex_coordinates;
    """

    post = """
    ASSERT_DOUBLE_EQ(1, 1); // FIXME
    """

    cell = ufl.Cell('tetrahedron', 3)
    r = None
    snippets = generate_jacobian_snippets(cell, r)
    code = format_code(snippets)

    gtest.add(pre + code + post)

# ...............................................................

def test_compute_jacobian_inverse_interval_1d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_interval(1);
    double * vertex_coordinates = mc.vertex_coordinates;
    """

    post = """
    ASSERT_DOUBLE_EQ(1, 1); // FIXME
    """

    cell = ufl.Cell('interval', 1)
    r = None
    snippets = []
    snippets.append( generate_jacobian_snippets(cell, r) )
    snippets.append( generate_jacobian_determinants_snippets(cell, r) )
    snippets.append( generate_jacobian_inverse_snippets(cell, r) )
    code = format_code(snippets)

    gtest.add(pre + code + post)

def test_compute_jacobian_inverse_interval_2d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_interval(2);
    double * vertex_coordinates = mc.vertex_coordinates;
    """

    post = """
    ASSERT_DOUBLE_EQ(1, 1); // FIXME
    """

    cell = ufl.Cell('interval', 2)
    r = None
    snippets = []
    snippets.append( generate_jacobian_snippets(cell, r) )
    snippets.append( generate_jacobian_determinants_snippets(cell, r) )
    snippets.append( generate_jacobian_inverse_snippets(cell, r) )
    code = format_code(snippets)

    gtest.add(pre + code + post)

def test_compute_jacobian_inverse_interval_3d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_interval(3);
    double * vertex_coordinates = mc.vertex_coordinates;
    """

    post = """
    ASSERT_DOUBLE_EQ(1, 1); // FIXME
    """

    cell = ufl.Cell('interval', 3)
    r = None
    snippets = []
    snippets.append( generate_jacobian_snippets(cell, r) )
    snippets.append( generate_jacobian_determinants_snippets(cell, r) )
    snippets.append( generate_jacobian_inverse_snippets(cell, r) )
    code = format_code(snippets)

    gtest.add(pre + code + post)

def test_compute_jacobian_inverse_triangle_2d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_triangle(2);
    double * vertex_coordinates = mc.vertex_coordinates;
    """

    post = """
    ASSERT_DOUBLE_EQ(1, 1); // FIXME
    """

    cell = ufl.Cell('triangle', 2)
    r = None
    snippets = []
    snippets.append( generate_jacobian_snippets(cell, r) )
    snippets.append( generate_jacobian_determinants_snippets(cell, r) )
    snippets.append( generate_jacobian_inverse_snippets(cell, r) )
    code = format_code(snippets)

    gtest.add(pre + code + post)

def test_compute_jacobian_inverse_triangle_3d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_triangle(3);
    double * vertex_coordinates = mc.vertex_coordinates;
    """

    post = """
    ASSERT_DOUBLE_EQ(1, 1); // FIXME
    """

    cell = ufl.Cell('triangle', 3)
    r = None
    snippets = []
    snippets.append( generate_jacobian_snippets(cell, r) )
    snippets.append( generate_jacobian_determinants_snippets(cell, r) )
    snippets.append( generate_jacobian_inverse_snippets(cell, r) )
    code = format_code(snippets)

    gtest.add(pre + code + post)

def test_compute_jacobian_inverse_tetrahedron_3d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_tetrahedron(3);
    double * vertex_coordinates = mc.vertex_coordinates;
    """

    post = """
    ASSERT_DOUBLE_EQ(1, 1); // FIXME
    """

    cell = ufl.Cell('tetrahedron', 3)
    r = None
    snippets = []
    snippets.append( generate_jacobian_snippets(cell, r) )
    snippets.append( generate_jacobian_determinants_snippets(cell, r) )
    snippets.append( generate_jacobian_inverse_snippets(cell, r) )
    code = format_code(snippets)

    gtest.add(pre + code + post)

# ...............................................................

def test_compute_x_from_xi_interval_1d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_interval(1);
    mc.scale(0.3); mc.translate(0.1);
    double * vertex_coordinates = mc.vertex_coordinates;
    double xi[1] = { 0.2 };
    """

    post = """
    ASSERT_DOUBLE_EQ(0.2*0.3+0.1, x[0]);
    """

    cell = ufl.Cell('interval', 1)
    r = None
    snippets = []
    snippets.append( generate_jacobian_snippets(cell, r) )
    snippets.append( generate_x_from_xi_snippets(cell, r) )
    code = format_code(snippets)

    gtest.add(pre + code + post)

def test_compute_x_from_xi_interval_2d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_interval(2);
    mc.scale(0.3); mc.translate(0.1, 0.4);
    double * vertex_coordinates = mc.vertex_coordinates;
    double xi[1] = { 0.2 };
    """

    post = """
    ASSERT_DOUBLE_EQ(0.2*0.3+0.1, x[0]);
    ASSERT_DOUBLE_EQ(0.4, x[1]);
    """

    cell = ufl.Cell('interval', 2)
    r = None
    snippets = []
    snippets.append( generate_jacobian_snippets(cell, r) )
    snippets.append( generate_x_from_xi_snippets(cell, r) )
    code = format_code(snippets)

    gtest.add(pre + code + post)

def test_compute_x_from_xi_interval_3d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_interval(3);
    mc.scale(0.3); mc.translate(0.1, 0.4, 0.5);

    double * vertex_coordinates = mc.vertex_coordinates;
    double xi[1] = { 0.2 };
    """

    post = """
    ASSERT_DOUBLE_EQ(0.2*0.3+0.1, x[0]);
    ASSERT_DOUBLE_EQ(0.4, x[1]);
    ASSERT_DOUBLE_EQ(0.5, x[2]);
    """

    cell = ufl.Cell('interval', 3)
    r = None
    snippets = []
    snippets.append( generate_jacobian_snippets(cell, r) )
    snippets.append( generate_x_from_xi_snippets(cell, r) )
    code = format_code(snippets)

    gtest.add(pre + code + post)

def test_compute_x_from_xi_triangle_2d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_triangle(2);
    mc.scale(0.3); mc.translate(0.1, 0.4);

    double * vertex_coordinates = mc.vertex_coordinates;
    double xi[2] = { 0.2, 0.6 };
    """

    post = """
    ASSERT_DOUBLE_EQ(0.2*0.3+0.1, x[0]);
    ASSERT_DOUBLE_EQ(0.6*0.3+0.4, x[1]);
    """

    cell = ufl.Cell('triangle', 2)
    r = None
    snippets = []
    snippets.append( generate_jacobian_snippets(cell, r) )
    snippets.append( generate_x_from_xi_snippets(cell, r) )
    code = format_code(snippets)

    gtest.add(pre + code + post)

def test_compute_x_from_xi_triangle_3d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_triangle(3);
    mc.scale(0.3); mc.translate(0.1, 0.4, 0.5);

    double * vertex_coordinates = mc.vertex_coordinates;
    double xi[2] = { 0.2, 0.6 };
    """

    post = """
    ASSERT_DOUBLE_EQ(0.2*0.3+0.1, x[0]);
    ASSERT_DOUBLE_EQ(0.6*0.3+0.4, x[1]);
    ASSERT_DOUBLE_EQ(0.5, x[2]);
    """

    cell = ufl.Cell('triangle', 3)
    r = None
    snippets = []
    snippets.append( generate_jacobian_snippets(cell, r) )
    snippets.append( generate_x_from_xi_snippets(cell, r) )
    code = format_code(snippets)

    gtest.add(pre + code + post)

def test_compute_x_from_xi_tetrahedron_3d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_tetrahedron(3);
    mc.scale(0.3); mc.translate(0.1, 0.4, 0.5);

    double * vertex_coordinates = mc.vertex_coordinates;
    double xi[3] = { 0.2, 0.6, 0.9 };
    """

    post = """
    ASSERT_DOUBLE_EQ(0.2*0.3+0.1, x[0]);
    ASSERT_DOUBLE_EQ(0.6*0.3+0.4, x[1]);
    ASSERT_DOUBLE_EQ(0.9*0.3+0.5, x[2]);
    """

    cell = ufl.Cell('tetrahedron', 3)
    r = None
    snippets = []
    snippets.append( generate_jacobian_snippets(cell, r) )
    snippets.append( generate_x_from_xi_snippets(cell, r) )
    code = format_code(snippets)

    gtest.add(pre + code + post)

# ...............................................................

def test_compute_xi_from_x_interval_1d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_interval(1);
    mc.scale(0.3); mc.translate(0.1);

    double * vertex_coordinates = mc.vertex_coordinates;
    double x[1] = { 0.2 };
    """

    post = """
    ASSERT_DOUBLE_EQ(1.0/0.3*(0.2-0.1), xi[0]);
    """

    cell = ufl.Cell('interval', 1)
    r = None
    snippets = []
    snippets.append( generate_jacobian_snippets(cell, r) )
    snippets.append( generate_jacobian_determinants_snippets(cell, r) )
    snippets.append( generate_jacobian_inverse_snippets(cell, r) )
    snippets.append( generate_xi_from_x_snippets(cell, r) )
    code = format_code(snippets)

    gtest.add(pre + code + post)

def test_compute_xi_from_x_interval_2d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_interval(2);
    mc.scale(0.3); mc.translate(0.1, 0.4);

    double * vertex_coordinates = mc.vertex_coordinates;
    double x[2] = { 0.2, 0.6 };
    """

    post = """
    ASSERT_DOUBLE_EQ(1.0/0.3*(0.2-0.1), xi[0]);
    """

    cell = ufl.Cell('interval', 2)
    r = None
    snippets = []
    snippets.append( generate_jacobian_snippets(cell, r) )
    snippets.append( generate_jacobian_determinants_snippets(cell, r) )
    snippets.append( generate_jacobian_inverse_snippets(cell, r) )
    snippets.append( generate_xi_from_x_snippets(cell, r) )
    code = format_code(snippets)

    gtest.add(pre + code + post)

def test_compute_xi_from_x_interval_3d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_interval(3);
    mc.scale(0.3); mc.translate(0.1, 0.4, 0.5);

    double * vertex_coordinates = mc.vertex_coordinates;
    double x[3] = { 0.2, 0.6, 0.9 };
    """

    post = """
    ASSERT_DOUBLE_EQ(1.0/0.3*(0.2-0.1), xi[0]);
    """

    cell = ufl.Cell('interval', 3)
    r = None
    snippets = []
    snippets.append( generate_jacobian_snippets(cell, r) )
    snippets.append( generate_jacobian_determinants_snippets(cell, r) )
    snippets.append( generate_jacobian_inverse_snippets(cell, r) )
    snippets.append( generate_xi_from_x_snippets(cell, r) )
    code = format_code(snippets)

    gtest.add(pre + code + post)

def test_compute_xi_from_x_triangle_2d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_triangle(2);
    mc.scale(0.3); mc.translate(0.1, 0.4);

    double * vertex_coordinates = mc.vertex_coordinates;
    double x[2] = { 0.2, 0.6 };
    """

    post = """
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
    code = format_code(snippets)

    gtest.add(pre + code + post)

def test_compute_xi_from_x_triangle_3d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_triangle(3);
    mc.scale(0.3); mc.translate(0.1, 0.4, 0.5);

    double * vertex_coordinates = mc.vertex_coordinates;
    double x[3] = { 0.2, 0.6, 0.9 };
    """

    post = """
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
    code = format_code(snippets)

    gtest.add(pre + code + post)

def test_compute_xi_from_x_tetrahedron_3d(gtest):
    """
    """

    pre = """
    mock_cell mc;
    mc.fill_reference_tetrahedron(3);
    mc.scale(0.3); mc.translate(0.1, 0.4, 0.5);

    double * vertex_coordinates = mc.vertex_coordinates;
    double x[3] = { 0.2, 0.6, 0.9 };
    """

    post = """
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
    code = format_code(snippets)

    gtest.add(pre + code + post)

