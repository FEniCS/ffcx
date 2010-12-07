"This program demonstrates how to plot finite elements with FFC."

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2010-12-07"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from ufl import *
from ffc import *

#element = FiniteElement("CG", triangle, 5)
#element = FiniteElement("CG", tetrahedron, 5)

#element = FiniteElement("CR", triangle, 1)
#element = FiniteElement("CR", tetrahedron, 1)

#element = FiniteElement("BDM", triangle, 2)
#element = FiniteElement("BDM", tetrahedron, 2)

#element = FiniteElement("RT", triangle, 2)
#element = FiniteElement("BDM", tetrahedron, 2)

element = FiniteElement("Nedelec 1st kind H(curl)", triangle, 3)
element = FiniteElement("Nedelec 1st kind H(curl)", tetrahedron, 3)

plot(element)
