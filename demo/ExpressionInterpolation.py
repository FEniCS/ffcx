# Copyright (C) 2022 Jørgen S. Dokken
#
# This file is part of FFCx.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# Defines an Expression which evaluates the several different functions at
# a set of interpolation points

import basix
from ufl import (Coefficient, FiniteElement, FunctionSpace, Mesh, MixedElement,
                 VectorElement, grad, triangle)

# Define mesh
cell = triangle
v_el = VectorElement("Lagrange", cell, 1)
mesh = Mesh(v_el)

# Define mixed function space
el = FiniteElement("CG", cell, 2)
el_int = VectorElement("Discontinuous Lagrange", cell, 1)
me = MixedElement([el, el_int])
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
b_rule = basix.quadrature.string_to_type(q_rule)
quadrature_points, _ = basix.make_quadrature(b_rule, b_cell, q_degree)

# Get interpolation points for output space
family = basix.finite_element.string_to_family("Lagrange", cell.cellname())
b_element = basix.create_element(family, b_cell, 4, basix.LagrangeVariant.gll_warped, True)
interpolation_points = b_element.points

# Create expressions that can be used for interpolation
expressions = [(du0, interpolation_points), (du1, interpolation_points),
               (powq, quadrature_points)]
