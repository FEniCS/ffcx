_author__ = "Anders Logg (logg@simula.no)"
__date__ = "2005-09-16 -- 2007-03-20"
__copyright__ = "Copyright (C) 2005-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Modified by Garth N. Wells 2006

# Python modules
import numpy

# FFC common modules
from ffc.common.debug import *

# FFC compiler.language modules
from ffc.compiler.language.algebra import *

# FFC FEM modules
from finiteelement import *

class MixedElement:
    """A MixedElement represents a vector-valued finite element
    created by the tensor product of a list of FiniteElements.

    Attributes:
        type_str               - string name of element type
        shape_str              - string name of element shape
        elements               - a list of finite elements
        mixed_basis            - a list of mixed basis functions
        mixed_degree           - maximum degree of basis functions
        mixed_space_dimension  - number of basis functions
        mixed cell_dimension         - number of shape dimensions
        mixed_tensordim        - number of components
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

        # Initialize data
        self.type_str = "mixed"
        self.shape_str = elements[0].shape_str

        # Compute degree
        self.mixed_degree = self.__compute_degree()
        # Compute number of basis functions
        self.mixed_space_dimension = self.__compute_space_dimension()
        # Compute number of shape dimensions
        self.mixed_cell_dimension = self.__compute_cell_dimension()
        # Compute number of components
        self.mixed_tensordim = self.__compute_tensordim()

        # FIXME: Temporary fix, need to figure this out. Don't really use it at this time.
        self.mappings = [element.mapping for element in elements]
        # self.mapping = [element.mapping for element in elements]

        # FIXME: Not implemented
        self.__entity_dofs = {}

    def signature(self):
        "Return a string identifying the finite element"
        return "Mixed finite element: [%s]" % ", ".join([element.signature() for element in self.elements])

    def basis(self):
        "Return basis of finite element space"
        raise RuntimeError, "Basis cannot be accessed explicitly for a mixed element."

    def degree(self):
        "Return degree of polynomial basis."
        return self.mixed_degree

    def cell_shape(self):
        "Return the cell shape"
        return self.elements[0].cell_shape()
    
    def facet_shape(self):
        "Return shape of facet."
        return self.elements[0].facet_shape()

    def space_dimension(self):
        "Return dimension of finite element space."
        return self.mixed_space_dimension

    def cell_dimension(self):
        "Return dimension of of shape."
        return self.mixed_cell_dimension

    def value_rank(self):
        "Return the rank of the value space"
        if len(self.elements) == 1 and self.elements[0].value_rank() == 0:
            return 0
        else:
            return 1

    def value_dimension(self, i):
        "Return the dimension of the value space for axis i"
        if self.value_rank() == 0:
            raise RuntimeError, "Cannot compute tensor dimension of scalar mixed element."
        elif not i == 0:
            raise RuntimeError, "Illegal tensor dimension for vector-valued mixed element."
        return self.mixed_tensordim

    def num_sub_elements(self):
        "Return the number of sub elements"
        return len(self.elements)

    def sub_element(self, i):
        "Return sub element i"
        return self.elements[i]

    def num_facets(self):
        "Return number of facets for shape of element."
        return self.elements[0].num_facets()

    def tabulate(self, order, points, facet = None):
        """Tabulate values on mixed element by appropriately reordering
        the tabulated values for the sub elements."""

        # Special case: only one element
        if len(self.elements) == 1:
            return elements[0].tabulate(order, points, facet)

        # Iterate over sub elements and build mixed table from element tables
        mixed_table = []
        offset = 0
        for i in range(len(self.elements)):
            # Get current element and table
            element = self.elements[i]
            table = element.tabulate(order, points, facet)
            # Iterate over the components corresponding to the current element
            if element.value_rank() == 0:
                component_table = self.__compute_component_table(table, offset)
                mixed_table.append(component_table)
            else:
                for i in range(element.value_dimension(0)):
                    component_table = self.__compute_component_table(table[i], offset)
                    mixed_table.append(component_table)
            # Add to offset, the number of the first basis function for the current element
            offset += element.space_dimension()

        return mixed_table

    def entity_dofs(self):
        """Return a dictionary mapping the mesh entities of the
        reference cell to the degrees of freedom associated with
        the entity"""
        return self.__entity_dofs

    def __compute_degree(self):
        "Compute maximum degree."
        return max([element.degree() for element in self.elements])

    def __compute_space_dimension(self):
        "Compute number of basis functions."
        return sum([element.space_dimension() for element in self.elements])

    def __compute_cell_dimension(self):
        "Compute number of shape dimensions."
        # Check that all elements are defined on the same shape
        for i in range(len(self.elements) - 1):
            e0 = self.elements[i]
            e1 = self.elements[i + 1]
            if not e0.cell_shape() == e1.cell_shape():
                raise FormError, ((e0, e1), "Elements defined on different shapes.")
        return self.elements[0].cell_dimension()
    
    def __compute_tensordim(self):
        "Compute number of components."
        sum = 0
        for element in self.elements:
            if element.value_rank() == 0:
                sum += 1
            elif element.value_rank() == 1:
                sum += element.value_dimension(0)
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
                num_points = numpy.shape(element_subtable)[1]
                mixed_subtable = numpy.zeros((self.mixed_space_dimension, num_points), dtype = numpy.float)
                # Iterate over element basis functions and fill in non-zero values
                for i in range(len(element_subtable)):
                    mixed_subtable[offset + i] = element_subtable[i]
                # Add to dictionary
                component_table[dorder][dtuple] = mixed_subtable
        return component_table

    def __add__(self, other):
        "Create mixed element"
        if isinstance(other, FiniteElement):
            return MixedElement(self.elements + [other])
        elif isinstance(other, MixedElement):
            return MixedElement(self.elements + other.elements)
        else:
            raise RuntimeError, "Unable to create mixed element from given object: " + str(other)

    def __repr__(self):
        "Pretty print"
        return "Mixed finite element: " + str(self.elements)

def BasisFunctions(element, functiontype = BasisFunction):
    "Create tuple of BasisFunctions from given MixedElement."
    if not isinstance(element, MixedElement):
        raise RuntimeError, "Basis function tuple must be created from mixed element."
    # Create basis function for mixed element
    vector = functiontype(element)
    # Pick components/subvectors of the mixed basis function
    subvectors = []
    offset = 0
    for e in element.elements:
        if e.value_rank() == 0:
            subvector = vector.pick_component_default(offset)
            offset += 1
        elif e.value_rank() == 1:
            if e.mapping == "Piola":
                subvector = [vector.pick_component_piola(k) for k in range(0, e.value_dimension(0))]
            else:
                subvector = [vector.pick_component_default(i) for i in range(offset, offset + e.value_dimension(0))]
            offset += e.value_dimension(0)
        else:
            raise RuntimeError, "Mixed elements can only be created from scalar or vector-valued elements."
        subvectors += [subvector]
    return tuple(subvectors)

def TestFunctions(element):
    "Create tuple of TestFunctions from given MixedElement."
    return BasisFunctions(element, TestFunction)

def TrialFunctions(element):
    "Create tuple of TrialFunctions from given MixedElement."
    return BasisFunctions(element, TrialFunction)

def Functions(element):
    "Create tuple of Functions from given MixedElement."
    if not isinstance(element, MixedElement):
        raise RuntimeError, "Function tuple must be created from mixed element."
    # Create function fox mixed element
    vector = Function(element)
    # Pick components/subvectors of the mixed basis function
    subvectors = []
    offset = 0
    for e in element.elements:
        if e.value_rank() == 0:
            subvector = vector[offset]
            offset += 1
        elif e.value_rank() == 1:
            subvector = [vector[i] for i in range(offset, offset + e.value_dimension(0))]
            offset += e.value_dimension(0)
        else:
            raise RuntimeError, "Mixed elements can only be created from scalar or vector-valued elements."
        subvectors += [subvector]
    return tuple(subvectors)
