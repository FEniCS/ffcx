#!/usr/bin/env py.test
# -*- coding: utf-8 -*-

import ufl
from ufl import *
#from ufl import product

from ffc.uflacs.backends.toy.toy_compiler import compile_expression

def compile_expression0(expr):
    code = ""
    return code

def compile_expression1(expr):
    code = compile_expression(expr, "")
    return code

def compile_tabulate_tensor_body(integral):
    # TODO: Handle measure type etc.
    if isinstance(integral, Form):
        integral, = integral.integrals()
    expr = integral.integrand()
    code = compile_expression(expr, "")
    return code

def test_interval_tabten_x_given(gtest):
    """Test code generation of body of the ufc function:

    void tabulate_tensor(
        double* A,
        const double * const * w,
        const double * coordinate_dofs,
        std::size_t num_points,
        const double * const * points
        ) const;

    """

    pre = """
    mock_cell mc;
    mc.fill_reference_interval(1);

    double * coordinate_dofs = mc.coordinate_dofs;
    double A[1] = { 0.0 };
    std::size_t num_points = 1;
    double points[1] = { 1.2 };

    // "Fetching" a point
    double * x = points;
    """

    post = """
    ASSERT_EQ(A[0], 1.2);
    """

    cell = interval

    x = SpatialCoordinate(cell)[0]
    expr = x
    integral = expr*dP

    code = compile_tabulate_tensor_body(integral)

    gtest.add(pre + code + post)

def test_interval_tabten_dg0_given(gtest):
    """Test code generation of body of the ufc function:

    void tabulate_tensor(
        double* A,
        const double * const * w,
        const double * coordinate_dofs,
        std::size_t num_points,
        const double * const * points
        ) const;

    with mock values for DG0/Real coefficients in w[][].

    """

    pre = """
    mock_cell mc;
    mc.fill_reference_interval(1);
    mc.coordinate_dofs[0*mc.geometric_dimension + 0] = 0.1;
    mc.coordinate_dofs[1*mc.geometric_dimension + 0] = 0.2;

    double A[1];
    memset(A, 0, sizeof(A));

    double w[2][2] = { { 1.2, 0.0 }, // second value here is unused
                       { 2.0, 3.0 } };

    double * coordinate_dofs = mc.coordinate_dofs;
    std::size_t num_points = 1;
    double points[1] = { 0.15 };

    // "Fetching" a point
    double * x = points;
    """

    post = """
    ASSERT_EQ(A[0], 0.15*1.2*(2.0+3.0));
    """

    cell = interval
    x = SpatialCoordinate(cell)[0]

    V = VectorElement("DG", cell, 0, dim=2)
    w0 = Constant(cell, count=0)
    w1 = Coefficient(V, count=1)

    expr = x*w0*(w1[0] + w1[1])

    integral = expr*dP

    code = compile_tabulate_tensor_body(integral)

    gtest.add(pre + code + post)

def xtest_interval_tabten_geometry_mappings(gtest):
    """Test code generation of body of the ufc function:

    void tabulate_tensor(
        double* A,
        const double * const * w,
        const double * coordinate_dofs
        ) const;
    """
    pass

def xtest_interval_geometry_expressions(gtest):
    """Test code generation of geometry expressions on an interval.

    """

    pre = """
    mock_cell mc;
    mc.fill_reference_interval(1);
    mc.coordinate_dofs[0*mc.geometric_dimension + 0] = 0.2;
    mc.coordinate_dofs[1*mc.geometric_dimension + 0] = 0.1;

    double * coordinate_dofs = mc.coordinate_dofs;
    double A[1];
    memset(A, 0, sizeof(A));

    """

    post = """
    // Check that geometric quantities are declared and computed correctly
    ASSERT_EQ(x[0], TODO);
    ASSERT_EQ(xi[0], TODO);
    ASSERT_EQ(J[0], -0.1);
    ASSERT_EQ(K[0], -10.0);
    ASSERT_EQ(detJ, -0.1);
    ASSERT_EQ(volume, 0.1);
    ASSERT_EQ(circumradius, TODO);

    // Check that geometric quantities have been placed in the output array
    std::size_t gd = mc.geometric_dimension;
    std::size_t td = mc.topological_dimension;
    double * AA = A;
    ASSERT_EQ(AA[0], x[0]); AA += gd;
    ASSERT_EQ(AA[0], xi[0]); AA += td;
    ASSERT_EQ(AA[0], J[0]); AA += gd * td;
    ASSERT_EQ(AA[0], K[0]); AA += td * gd;
    ASSERT_EQ(AA[0], detJ); AA += 1;
    ASSERT_EQ(AA[0], volume); AA += 1;
    ASSERT_EQ(AA[0], circumradius); AA += 1;
    """
    cell = interval

    gd = cell.geometric_dimension()
    td = cell.topological_dimension()

    values = []
    values.extend(SpatialCoordinate(cell)[i] for i in range(gd))
    values.extend(LocalCoordinate(cell)[i] for i in range(td))
    values.extend(GeometryJacobian(cell)[i, j] for i in range(gd) for i in range(td))
    values.extend(InverseGeometryJacobian(cell)[i, j] for i in range(td) for i in range(gd))
    values.append(GeometryJacobianDeterminant(cell))
    values.append(CellVolume(cell))
    values.append(Circumradius(cell))

    expr = as_vector(values)

    code = compile_expression1(expr)

    gtest.add(pre + code + post)

def test_tabulate_tensor_interval_facet(gtest):
    """Test code generation of body of the ufc function:

    void tabulate_tensor(
        double* A,
        const double * const * w,
        const double * coordinate_dofs,
        unsigned int facet) const;

    """

    pre = """
    mock_cell mc;
    mc.fill_reference_interval(1);
    mc.coordinate_dofs[0*mc.geometric_dimension + 0] = 0.1;
    mc.coordinate_dofs[1*mc.geometric_dimension + 0] = 0.2;

    double A[1];
    memset(A, 0, sizeof(A));

    double w[1][2] = { { 2.0, 3.0 } };

    double * coordinate_dofs = mc.coordinate_dofs;

    unsigned int facet = 0;
    """

    post = """
    // TODO
    """

    code = "// TODO"

    gtest.add(pre + code + post)

def test_tabulate_tensor_interval_interior_facet(gtest):
    """Test code generation of body of the ufc function:

    void tabulate_tensor(
        double* A,
        const double * const * w,
        const double * coordinate_dofs0,
        const double * coordinate_dofs1,
        unsigned int facet0,
        unsigned int facet1) const;
    """

    pre = """
    mock_cell mc0;
    mc0.fill_reference_interval(1);
    mc0.coordinate_dofs[0*mc0.geometric_dimension + 0] = 0.1;
    mc0.coordinate_dofs[1*mc0.geometric_dimension + 0] = 0.2;

    mock_cell mc1;
    mc1.fill_reference_interval(1);
    mc1.coordinate_dofs[0*mc1.geometric_dimension + 0] = 0.2;
    mc1.coordinate_dofs[1*mc1.geometric_dimension + 0] = 0.3;

    double A[2];
    memset(A, 0, sizeof(A));

    double w[1][4] = { { 2.0, 3.0,
                         4.0, 5.0 } };

    double * coordinate_dofs0 = mc0.coordinate_dofs;
    double * coordinate_dofs1 = mc1.coordinate_dofs;

    unsigned int facet0 = 1;
    unsigned int facet1 = 0;
    """

    post = """
    // TODO
    """

    code = "// TODO"

    gtest.add(pre + code + post)
