# Copyright (C) 2021 Igor Baratta
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
# The bilinear form for a mass matrix.
from ufl import (FiniteElement, TestFunction, TrialFunction, dx, inner,
                 tetrahedron)

element = FiniteElement("DG", tetrahedron, 0)

v = TestFunction(element)
u = TrialFunction(element)

a = inner(u, v) * dx
