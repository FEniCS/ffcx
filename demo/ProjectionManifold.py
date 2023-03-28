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
#
# This demo illustrates use of finite element spaces defined over
# simplicies embedded in higher dimensions
import basix.ufl
from ufl import TestFunctions, TrialFunctions, div, dx, inner

# Define element over this domain
V = basix.ufl.element("RT", "triangle", 1, gdim=3)
Q = basix.ufl.element("DG", "triangle", 0, gdim=3)
element = basix.ufl.mixed_element([V, Q])

(u, p) = TrialFunctions(element)
(v, q) = TestFunctions(element)

a = (inner(u, v) + div(u) * q + div(v) * p) * dx
