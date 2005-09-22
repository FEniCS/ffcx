__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2005-09-16 -- 2005-09-20"
__copyright__ = "Copyright (c) 2005 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
import Numeric

# FFC common modules
from ffc.common.debug import *

# FFC compiler modules
import algebra
import finiteelement
from dofmap import *
from pointmap import *
from vertexeval import *

def BasisFunctions(element):
    "Create tuple of BasisFunction from given MixedElement."
    if not isinstance(element, MixedElement):
        raise RuntimeError, "Basis function tuple must be created from mixed element."
    # Create basis function fox mixed element
    vector = algebra.BasisFunction(element)
    # Pick components/subvectors of the mixed basis function
    subvectors = []
    offset = 0
    for e in element.elements:
        if e.rank() == 0:
            subvector = vector[offset]
            offset += 1
        elif e.rank() == 1:
            subvector = [vector[i] for i in range(offset, offset + e.tensordim(0))]
            offset += e.tensordim(0)
        else:
            raise RuntimeError, "Mixed elements can only be created from scalar or vector-valued elements."
        subvectors += [subvector]
    return tuple(subvectors)

class MixedElement:
    """A MixedElement represents a vector-valued finite element
    created by the tensor product of a list of FiniteElements.

    Attributes:
        elements        - a list of finite elements
        mixed_basis     - a list of mixed basis functions
        mixed_degree    - maximum degree of basis functions
        mixed_spacedim  - number of basis functions
        mixed shapedim  - number of shape dimensions
        mixed_tensordim - number of components
    """

    def __init__(self, elements):
        "Create MixedElement from a list of elements."

        # Create list of elements
        if not isinstance(elements, list):
            self.elements = [elements]
        else:
            self.elements = elements

        # Check that we have at least one element
        if not len(self.elements) > 0:
            raise FormError, "Mixed finite element must contain at least one element."

        # Compute degree
        self.mixed_degree = self.__compute_degree()
        # Compute number of basis functions
        self.mixed_spacedim = self.__compute_spacedim()
        # Compute number of shape dimensions
        self.mixed_shapedim = self.__compute_shapedim()
        # Compute number of components
        self.mixed_tensordim = self.__compute_tensordim()

        # Create dof map
        self.dofmap = DofMap(self.elements)

        # Create point map
        self.pointmap = PointMap(self.elements)

        # Create vertex evaluation
        self.vertexeval = VertexEval(self.elements)

    def basis(self):
        "Return basis of finite element space."
        raise RuntimeError, "Basis cannot be accessed explicitly for mixed element."

    def degree(self):
        "Return degree of polynomial basis."
        return self.mixed_degree

    def shape(self):
        "Return shape used for element."
        return self.elements[0].shape()

    def spacedim(self):
        "Return dimension of finite element space."
        return self.mixed_spacedim

    def shapedim(self):
        "Return dimension of of shape."
        return self.mixed_shapedim

    def rank(self):
        "Return rank of basis functions."
        if len(self.elements) == 1 and self.elements[0].rank() == 0:
            return 0
        else:
            return 1

    def tensordim(self, i):
        "Return size of given dimension."
        if self.rank() == 0:
            raise RuntimeError, "Cannot compute tensor dimension of scalar mixed element."
        elif not i == 0:
            raise RuntimeError, "Illegal tensor dimension for vector-valued mixed element."
        return self.mixed_tensordim

    def tabulate(self, order, points):
        """Return tabulated values of derivatives up to given order of
        basis functions at given points."""
        # Special case: only one element
        if len(self.elements) == 1:
            return elements[0].tabulate(order, points)
        # Iterate over elements and build mixed table from element tables.
        # This is a bit nasty, so it needs some thought...
        mixed_table = []
        offset = 0
        for i in range(len(self.elements)):
            # Get current element and table
            element = self.elements[i]
            table = element.tabulate(order, points)
            # Iterate over the components corresponding to the current element
            if element.rank() == 0:
                component_table = self.__compute_component_table(table, offset)
                mixed_table.append(component_table)
            else:
                for i in range(element.tensordim(0)):
                    component_table = self.__compute_component_table(table[i], offset)
                    mixed_table.append(component_table)
            # Add to offset, the number of the first basis function for the current element
            offset += element.spacedim()

        return mixed_table

#    def __create_basis(self):
#        "Create basis for mixed element."
#        raise RuntimeError, "Not implemented."

    def __compute_degree(self):
        "Compute maximum degree."
        return max([element.degree() for element in self.elements])

    def __compute_spacedim(self):
        "Compute number of basis functions."
        return sum([element.spacedim() for element in self.elements])

    def __compute_shapedim(self):
        "Compute number of shape dimensions."
        # Check that all elements are defined on the same shape
        for i in range(len(self.elements) - 1):
            e0 = self.elements[i]
            e1 = self.elements[i + 1]
            if not e0.shape() == e1.shape():
                raise FormError, ((e0, e1), "Elements defined on different shapes.")
        return self.elements[0].shapedim()
    
    def __compute_tensordim(self):
        "Compute number of components."
        sum = 0
        for element in self.elements:
            if element.rank() == 0:
                sum += 1
            elif element.rank() == 1:
                sum += element.tensordim(0)
            else:
                raise RuntimeError, "Mixed elements can only be created from scalar or vector-valued elements."
        return sum

    def __compute_component_table(self, table, offset):
        "Compute subtable for given component."
        component_table = []
        # Iterate over derivative orders
        for dorder in range(len(table)):
            component_table.append({})
            # Iterate over derivative tuples
            derivative_dictionary = {}
            for dtuple in table[dorder]:
                element_subtable = table[dorder][dtuple]
                num_points = Numeric.shape(element_subtable)[1]
                mixed_subtable = Numeric.zeros((self.mixed_spacedim, num_points), Numeric.Float)
                # Iterate over element basis functions and fill in non-zero values
                for i in range(len(element_subtable)):
                    mixed_subtable[offset + i] = element_subtable[i]
                # Add to dictionary
                component_table[dorder][dtuple] = mixed_subtable
        return component_table

    def __add__(self, other):
        "Create mixed element."
        if isinstance(other, finiteelement.FiniteElement):
            return MixedElement(self.elements + [other])
        elif isinstance(other, MixedElement):
            return MixedElement(self.elements + other.elements)
        else:
            raise RuntimeError, "Unable to create mixed element from given object: " + str(other)

    def __repr__(self):
        "Print nicely formatted representation of MixedElement."
        return "Mixed finite element: " + str(self.elements)
