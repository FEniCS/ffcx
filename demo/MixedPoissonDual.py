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
#
# First added:  2014-01-29
# Last changed: 2014-01-29
#
# The bilinear form a(u, v) and linear form L(v) for a two-field
# (mixed) formulation of Poisson's equation
from basix.ufl import mixed_element, element
from ufl import Coefficient, TestFunctions, TrialFunctions, dot, ds, dx, grad

DRT = element("Discontinuous RT", "triangle", 2)
P = element("P", "triangle", 3)
W = mixed_element([DRT, P])

(sigma, u) = TrialFunctions(W)
(tau, v) = TestFunctions(W)

P1 = element("P", "triangle", 1)
f = Coefficient(P1)
g = Coefficient(P1)

a = (dot(sigma, tau) + dot(grad(u), tau) + dot(sigma, grad(v))) * dx
L = - f * v * dx - g * v * ds
