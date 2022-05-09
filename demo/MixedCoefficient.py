# -*- coding: utf-8 -*-

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
#
# Mixed coefficient.
from ufl import (Coefficients, FiniteElement, MixedElement, VectorElement, dot,
                 dS, dx, triangle)

cell = triangle

DG = VectorElement("DG", cell, 0)
CG = FiniteElement("Lagrange", cell, 2)
RT = FiniteElement("RT", cell, 3)

element = MixedElement(DG, CG, RT)

f, g, h = Coefficients(element)

forms = [dot(f('+'), h('-')) * dS + g * dx]
