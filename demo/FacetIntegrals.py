# Copyright (C) 2009-2010 Anders Logg
#
# This file is part of FFCx.
#
# FFCx is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFCx is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFCx. If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2009-03-20
# Last changed: 2011-03-08
#
# Simple example of a form defined over exterior and interior facets.
from ufl import (FacetNormal, FiniteElement, TestFunction, TrialFunction, avg,
                 ds, dS, grad, inner, jump, triangle)

element = FiniteElement("Discontinuous Lagrange", triangle, 1)

u = TrialFunction(element)
v = TestFunction(element)

n = FacetNormal(triangle)

a = u * v * ds \
    + u('+') * v('-') * dS \
    + inner(jump(u, n), avg(grad(v))) * dS \
    + inner(avg(grad(u)), jump(v, n)) * dS
