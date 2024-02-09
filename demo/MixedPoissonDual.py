# Copyright (C) 2014 Jan Blechta
#
# This file is part of FFCx.
#
# DOLFINx is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFINx is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFINx. If not, see <http://www.gnu.org/licenses/>.
"""Mixed Poisson dual demo.

The bilinear form a(u, v) and linear form L(v) for a two-field
(mixed) formulation of Poisson's equation.
"""

import basix.ufl
from ufl import Coefficient, FunctionSpace, Mesh, TestFunctions, TrialFunctions, ds, dx, grad, inner

DRT = basix.ufl.element("Discontinuous RT", "triangle", 2)
P = basix.ufl.element("P", "triangle", 3)
W = basix.ufl.mixed_element([DRT, P])
domain = Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
space = FunctionSpace(domain, W)

(sigma, u) = TrialFunctions(space)
(tau, v) = TestFunctions(space)

P1 = basix.ufl.element("P", "triangle", 1)
space = FunctionSpace(domain, P1)
f = Coefficient(space)
g = Coefficient(space)

a = (inner(sigma, tau) + inner(grad(u), tau) + inner(sigma, grad(v))) * dx
L = -inner(f, v) * dx - inner(g, v) * ds
