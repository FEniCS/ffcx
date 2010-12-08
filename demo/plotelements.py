"This program demonstrates how to plot finite elements with FFC."

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2010-12-07"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

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
#element = FiniteElement("Nedelec 2nd kind H(curl)", tetrahedron, 3)

#element = FiniteElement("Raviart-Thomas", triangle, 3)
#element = FiniteElement("Raviart-Thomas", tetrahedron, 3)

plot(element)
