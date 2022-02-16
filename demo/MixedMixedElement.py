# Copyright (C) 2007 Anders Logg
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
# A mixed element of mixed elements
#
# Compile this form with FFCx: ffcx MixedMixedElement.ufl
from ufl import FiniteElement, TestFunction, VectorElement, dx, triangle

cell = triangle

DG = VectorElement("DG", cell, 0)
CG = FiniteElement("Lagrange", cell, 2)
RT = FiniteElement("RT", cell, 3)

element = DG * CG * RT

v = TestFunction(element)
a = v[0] * dx
