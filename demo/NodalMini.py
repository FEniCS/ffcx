# Copyright (C) 2018 Jan Blechta
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
from ufl import (FiniteElement, MixedElement, NodalEnrichedElement,
                 TestFunctions, TrialFunctions, VectorElement, div, dx, grad,
                 inner, triangle)

cell = triangle
P1 = FiniteElement("P", cell, 1)
B3 = FiniteElement("B", cell, 3)
V = VectorElement(NodalEnrichedElement(P1, B3))
element = MixedElement(V, P1)

u, p = TrialFunctions(element)
v, q = TestFunctions(element)

a = (inner(grad(u), grad(v)) - div(v) * p - div(u) * q) * dx
