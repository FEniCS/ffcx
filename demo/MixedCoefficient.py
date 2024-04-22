# Copyright (C) 2016 Miklós Homolya
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
"""Mixed coefficient demo."""

import basix.ufl
from ufl import Coefficients, FunctionSpace, Mesh, dot, dS, dx

DG = basix.ufl.element("DG", "triangle", 0, shape=(2,))
CG = basix.ufl.element("Lagrange", "triangle", 2)
RT = basix.ufl.element("RT", "triangle", 3)

element = basix.ufl.mixed_element([DG, CG, RT])
domain = Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
space = FunctionSpace(domain, element)

f, g, h = Coefficients(space)

forms = [dot(f("+"), h("-")) * dS + g * dx]
