# Copyright (C) 2022 Jørgen S. Dokken
#
# This file is part of FFCx.
#
# SPDX-License-Identifier:     LGPL-3.0-or-later
#
# Defines an Expression which evaluates the several different functions at
# a set of interpolation points

import basix
from ufl import (Coefficient, FiniteElement, FunctionSpace, Mesh, MixedElement,
                 VectorElement, grad, triangle)

# Define mesh
cell = triangle
mesh = Mesh(VectorElement("Lagrange", cell, 1))

# Define mixed function space
me = MixedElement([FiniteElement("CG", cell, 2), VectorElement("Discontinuous Lagrange", cell, 1)])
V = FunctionSpace(mesh, me)
u = Coefficient(V)

# Define expressions on each sub-space
du0 = grad(u[0])
du1 = grad(u[1])

# Define an expression using quadrature elements
q_rule = "gauss_jacobi"
q_degree = 3
q_el = FiniteElement("Quadrature", cell, q_degree, quad_scheme=q_rule)
Q = FunctionSpace(mesh, q_el)
q = Coefficient(Q)
powq = 3 * q**2

# Extract basix cell type
b_cell = basix.cell.string_to_type(cell.cellname())

# Find quadrature points for quadrature element
b_rule = basix.quadrature.string_to_type(q_el.quadrature_scheme())
quadrature_points, _ = basix.make_quadrature(b_rule, b_cell, q_el.degree())

# Get interpolation points for output space
family = basix.finite_element.string_to_family("Lagrange", cell.cellname())
output_degree = 4
b_element = basix.create_element(family, b_cell, output_degree, basix.LagrangeVariant.gll_warped, True)
X = b_element.points

# Create expressions that can be used for interpolation
expressions = [(du0, X), (du1, X), (powq, quadrature_points)]
