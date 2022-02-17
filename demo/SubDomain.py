# Copyright (C) 2008 Anders Logg
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
# This example illustrates how to define a form over a
# given subdomain of a mesh, in this case a functional.
from ufl import (Coefficient, FiniteElement, TestFunction, TrialFunction, ds,
                 dx, tetrahedron)

element = FiniteElement("CG", tetrahedron, 1)

v = TestFunction(element)
u = TrialFunction(element)
f = Coefficient(element)

M = f * dx(2) + f * ds(5)
