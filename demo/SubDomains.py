# Copyright (C) 2008 Anders Logg and Kristian B. Oelgaard
#
# This file is part of UFL.
#
# UFL is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFL. If not, see <http://www.gnu.org/licenses/>.
#
# This simple example illustrates how forms can be defined on different sub domains.
# It is supported for all three integral types.
from ufl import (FiniteElement, TestFunction, TrialFunction, ds, dS, dx,
                 tetrahedron)

element = FiniteElement("CG", tetrahedron, 1)

v = TestFunction(element)
u = TrialFunction(element)

a = v * u * dx(0) + 10.0 * v * u * dx(1) + v * u * ds(0) + 2.0 * v * u * ds(1)\
    + v('+') * u('+') * dS(0) + 4.3 * v('+') * u('+') * dS(1)
