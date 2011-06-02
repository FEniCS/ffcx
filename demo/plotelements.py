"This program demonstrates how to plot finite elements with FFC."

# Copyright (C) 2010 Anders Logg
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2010-12-07
# Last changed: 2010-12-10

from ufl import *
from ffc import *

#element = FiniteElement("Argyris", triangle, 5)

#element = FiniteElement("Arnold-Winther", triangle)

#element = FiniteElement("Brezzi-Douglas-Marini", triangle, 3)
element = FiniteElement("Brezzi-Douglas-Marini", tetrahedron, 3)

#element = FiniteElement("Crouzeix-Raviart", triangle, 1)
#element = FiniteElement("Crouzeix-Raviart", tetrahedron, 1)

#element = FiniteElement("Discontinuous Lagrange", triangle, 3)
#element = FiniteElement("Discontinuous Lagrange", tetrahedron, 3)

#element = FiniteElement("Hermite", triangle)
#element = FiniteElement("Hermite", tetrahedron)

#element = FiniteElement("Lagrange", triangle, 3)
#element = FiniteElement("Lagrange", tetrahedron, 3)

#element = FiniteElement("Mardal-Tai-Winther", triangle)

#element = FiniteElement("Morley", triangle)

#element = FiniteElement("Nedelec 1st kind H(curl)", triangle, 3)
#element = FiniteElement("Nedelec 1st kind H(curl)", tetrahedron, 3)

#element = FiniteElement("Nedelec 2nd kind H(curl)", triangle, 3)
#element = FiniteElement("Nedelec 2nd kind H(curl)", tetrahedron, 1)

#element = FiniteElement("Raviart-Thomas", triangle, 3)
#element = FiniteElement("Raviart-Thomas", tetrahedron, 3)

plot(element)
#plot(element, rotate=False)
#plot("notation")
