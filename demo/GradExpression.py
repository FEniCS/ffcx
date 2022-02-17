# Copyright (C) 2021 Michal Habera
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
# Defines an Expression which evaluates gradient of an Argument
# at predefined set of points.
from ufl import (Coefficient, Constant, FiniteElement, FunctionSpace, Mesh,
                 TrialFunction, VectorElement, grad, triangle)

element = FiniteElement("Lagrange", triangle, 2)
coordinate_el = VectorElement("Lagrange", triangle, 1)
mesh = Mesh(coordinate_el)

V = FunctionSpace(mesh, element)

u = TrialFunction(V)
f = Coefficient(V)
c = Constant(mesh)

grad_u = c * f * grad(u)

expressions = [(grad_u, [[0.0, 0.1], [0.2, 0.3]]),
               (c * f * u, [[0.0, 1.0], [1.0, 0.0]])]
