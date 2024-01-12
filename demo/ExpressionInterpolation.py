# Copyright (C) 2022 Jørgen S. Dokken
#
# This file is part of FFCx.
#
# SPDX-License-Identifier:     LGPL-3.0-or-later
#
# Defines an Expression which evaluates the several different functions at
# a set of interpolation points

import basix
import basix.ufl
from ufl import Coefficient, FunctionSpace, Mesh, grad

# Define mesh
cell = "triangle"
v_el = basix.ufl.element("Lagrange", cell, 1, shape=(2, ))
mesh = Mesh(v_el)

# Define mixed function space
el = basix.ufl.element("P", cell, 2)
el_int = basix.ufl.element("Discontinuous Lagrange", cell, 1, shape=(2, ))
me = basix.ufl.mixed_element([el, el_int])

V = FunctionSpace(mesh, me)
u = Coefficient(V)

# Define expressions on each sub-space
du0 = grad(u[0])
du1 = grad(u[1])

# Define an expression using quadrature elements
q_rule = "gauss_jacobi"
q_degree = 3
q_el = basix.ufl.quadrature_element(cell, scheme=q_rule, degree=q_degree)
Q = FunctionSpace(mesh, q_el)
q = Coefficient(Q)
powq = 3 * q**2

# Extract basix cell type
b_cell = basix.cell.string_to_type(cell)

# Find quadrature points for quadrature element
b_rule = basix.quadrature.string_to_type(q_rule)
quadrature_points, _ = basix.quadrature.make_quadrature(b_cell, q_degree, rule=b_rule)

# Get interpolation points for output space
family = basix.finite_element.string_to_family("Lagrange", cell)
b_element = basix.create_element(family, b_cell, 4, basix.LagrangeVariant.gll_warped, discontinuous=True)
interpolation_points = b_element.points


# Create expressions that can be used for interpolation
expressions = [(du0, X), (du1, X), (powq, quadrature_points)]
