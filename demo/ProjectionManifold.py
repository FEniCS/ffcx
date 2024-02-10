# Copyright (C) 2012 Marie E. Rognes and David Ham
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
"""Projection manifold demo.

This demo illustrates use of finite element spaces defined over
simplicies embedded in higher dimensions.
"""

import basix.ufl
from ufl import FunctionSpace, Mesh, TestFunctions, TrialFunctions, div, dx, inner

# Define element over this domain
V = basix.ufl.element("RT", "triangle", 1)
Q = basix.ufl.element("DG", "triangle", 0)
element = basix.ufl.mixed_element([V, Q])
domain = Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(3,)))
space = FunctionSpace(domain, element)

(u, p) = TrialFunctions(space)
(v, q) = TestFunctions(space)

a = (inner(u, v) + inner(div(u), q) + inner(p, div(v))) * dx
