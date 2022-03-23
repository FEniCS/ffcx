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
from ufl import (Cell, FiniteElement, TestFunctions, TrialFunctions, div, dx,
                 inner)

# Define interval embedded in 3D:
domain = Cell("triangle", geometric_dimension=3)

# Define element over this domain
V = FiniteElement("RT", domain, 1)
Q = FiniteElement("DG", domain, 0)
element = V * Q

(u, p) = TrialFunctions(element)
(v, q) = TestFunctions(element)

a = (inner(u, v) + div(u) * q + div(v) * p) * dx
