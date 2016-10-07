#!/usr/bin/env py.test
# -*- coding: utf-8 -*-

"""
These tests check that the mock implementations
of ufc cells are properly constructed. The mock cells
are not actual subclasses of ufc::cell, instead using
sufficiently large fixed size arrays to avoid inconvenient
memory allocation code in the unit tests. The tests here
can be seen as a demonstration of how to use the mock cells.
"""

def test_mock_interval(gtest):
    pre = """
    mock_cell mc;
    mc.fill_reference_interval(1);
    double * coordinate_dofs = mc.coordinate_dofs;
    """

    post = """
    ASSERT_EQ(coordinate_dofs[0*mc.geometric_dimension + 0], 0.0);
    ASSERT_EQ(coordinate_dofs[1*mc.geometric_dimension + 0], 2.0);
    """

    code = """
    mc.scale(2.0);
    """

    gtest.add(pre + code + post)

def test_mock_triangle(gtest):
    pre = """
    mock_cell mc;
    mc.fill_reference_triangle(2);
    double * coordinate_dofs = mc.coordinate_dofs;
    """

    post = """
    ASSERT_EQ(coordinate_dofs[0*mc.geometric_dimension + 0], 0.0);
    ASSERT_EQ(coordinate_dofs[0*mc.geometric_dimension + 1], 0.0);
    ASSERT_EQ(coordinate_dofs[1*mc.geometric_dimension + 0], 2.0);
    ASSERT_EQ(coordinate_dofs[1*mc.geometric_dimension + 1], 0.0);
    ASSERT_EQ(coordinate_dofs[2*mc.geometric_dimension + 0], 0.0);
    ASSERT_EQ(coordinate_dofs[2*mc.geometric_dimension + 1], 3.0);
    """

    code = """
    double factors[2] = { 2.0, 3.0 };
    mc.scale(factors);
    """

    gtest.add(pre + code + post)

def test_mock_tetrahedron(gtest):
    pre = """
    mock_cell mc;
    mc.fill_reference_tetrahedron(3);
    double * coordinate_dofs = mc.coordinate_dofs;
    """

    post = """
    ASSERT_EQ(coordinate_dofs[0*mc.geometric_dimension + 0], 2.0);
    ASSERT_EQ(coordinate_dofs[0*mc.geometric_dimension + 1], 3.0);
    ASSERT_EQ(coordinate_dofs[0*mc.geometric_dimension + 2], 4.0);
    ASSERT_EQ(coordinate_dofs[1*mc.geometric_dimension + 0], 3.0);
    ASSERT_EQ(coordinate_dofs[1*mc.geometric_dimension + 1], 3.0);
    ASSERT_EQ(coordinate_dofs[1*mc.geometric_dimension + 2], 4.0);
    ASSERT_EQ(coordinate_dofs[2*mc.geometric_dimension + 0], 2.0);
    ASSERT_EQ(coordinate_dofs[2*mc.geometric_dimension + 1], 4.0);
    ASSERT_EQ(coordinate_dofs[2*mc.geometric_dimension + 2], 4.0);
    ASSERT_EQ(coordinate_dofs[3*mc.geometric_dimension + 0], 2.0);
    ASSERT_EQ(coordinate_dofs[3*mc.geometric_dimension + 1], 3.0);
    ASSERT_EQ(coordinate_dofs[3*mc.geometric_dimension + 2], 5.0);
    """

    code = """
    double offset[3] = { 2.0, 3.0, 4.0 };
    mc.translate(offset);
    """

    gtest.add(pre + code + post)

def test_mock_quadrilateral(gtest):
    pre = """
    mock_cell mc;
    mc.fill_reference_quadrilateral(2);
    double * coordinate_dofs = mc.coordinate_dofs;
    """

    post = """
    ASSERT_EQ(coordinate_dofs[0*mc.geometric_dimension + 0], 2.0 + 2.0*0.0 + 3.0*0.0);
    ASSERT_EQ(coordinate_dofs[0*mc.geometric_dimension + 1], 3.0 + 4.0*0.0 + 5.0*0.0);
    ASSERT_EQ(coordinate_dofs[1*mc.geometric_dimension + 0], 2.0 + 2.0*1.0 + 3.0*0.0);
    ASSERT_EQ(coordinate_dofs[1*mc.geometric_dimension + 1], 3.0 + 4.0*1.0 + 5.0*0.0);
    ASSERT_EQ(coordinate_dofs[2*mc.geometric_dimension + 0], 2.0 + 2.0*1.0 + 3.0*1.0);
    ASSERT_EQ(coordinate_dofs[2*mc.geometric_dimension + 1], 3.0 + 4.0*1.0 + 5.0*1.0);
    ASSERT_EQ(coordinate_dofs[3*mc.geometric_dimension + 0], 2.0 + 2.0*0.0 + 3.0*1.0);
    ASSERT_EQ(coordinate_dofs[3*mc.geometric_dimension + 1], 3.0 + 4.0*0.0 + 5.0*1.0);
    """

    code = """
    double G[4] = { 2.0, 3.0,
                    4.0, 5.0 };
    double x[2] = { 2.0, 3.0 };
    mc.affine_map(G, x);
    """

    gtest.add(pre + code + post)

def test_mock_hexahedron(gtest):
    pre = """
    mock_cell mc;
    mc.fill_reference_hexahedron(3);
    double * coordinate_dofs = mc.coordinate_dofs;
    """

    post = """
    ASSERT_EQ(coordinate_dofs[0*mc.geometric_dimension + 0], 5.0 * (0.0 + 2.0));
    ASSERT_EQ(coordinate_dofs[0*mc.geometric_dimension + 1], 6.0 * (0.0 + 3.0));
    ASSERT_EQ(coordinate_dofs[0*mc.geometric_dimension + 2], 7.0 * (0.0 + 4.0));
    ASSERT_EQ(coordinate_dofs[1*mc.geometric_dimension + 0], 5.0 * (1.0 + 2.0));
    ASSERT_EQ(coordinate_dofs[1*mc.geometric_dimension + 1], 6.0 * (0.0 + 3.0));
    ASSERT_EQ(coordinate_dofs[1*mc.geometric_dimension + 2], 7.0 * (0.0 + 4.0));
    ASSERT_EQ(coordinate_dofs[2*mc.geometric_dimension + 0], 5.0 * (1.0 + 2.0));
    ASSERT_EQ(coordinate_dofs[2*mc.geometric_dimension + 1], 6.0 * (1.0 + 3.0));
    ASSERT_EQ(coordinate_dofs[2*mc.geometric_dimension + 2], 7.0 * (0.0 + 4.0));
    // ...
    ASSERT_EQ(coordinate_dofs[7*mc.geometric_dimension + 0], 5.0 * (0.0 + 2.0));
    ASSERT_EQ(coordinate_dofs[7*mc.geometric_dimension + 1], 6.0 * (1.0 + 3.0));
    ASSERT_EQ(coordinate_dofs[7*mc.geometric_dimension + 2], 7.0 * (1.0 + 4.0));
    """

    code = """
    double offset[3] = { 2.0, 3.0, 4.0 };
    double factors[3] = { 5.0, 6.0, 7.0 };
    mc.translate(offset);
    mc.scale(factors);
    """

    gtest.add(pre + code + post)

