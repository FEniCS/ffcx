# Copyright (C) 2007 Anders Logg
#
# This file is part of UFL.
#
# UFL is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFL. If not, see <http://www.gnu.org/licenses/>.
#
# Specification of a system F(v, u) = 0 and extraction of
# the bilinear and linear forms a and L for the left- and
# right-hand sides:
#
#   F(v, u) = a(v, u) - L(v) = 0
#
# The example below demonstrates the specification of the
# linear system for a cG(1)/Crank-Nicholson time step for
# the heat equation.
#
# The below formulation is equivalent to writing
#
#  a = v*u*dx + 0.5*k*dot(grad(v), grad(u))*dx
#  L = v*u0*dx - 0.5*k*dot(grad(v), grad(u0))*dx
#
# but instead of manually shuffling terms not including
# the unknown u to the right-hand side, all terms may
# be listed on one line and left- and right-hand sides
# extracted by lhs() and rhs().
from ufl import (Coefficient, FiniteElement, TestFunction, TrialFunction, dot,
                 dx, grad, lhs, rhs, triangle)

element = FiniteElement("Lagrange", triangle, 1)

k = 0.1

v = TestFunction(element)
u = TrialFunction(element)
u0 = Coefficient(element)

F = v * (u - u0) * dx + k * dot(grad(v), grad(0.5 * (u0 + u))) * dx

a = lhs(F)
L = rhs(F)
