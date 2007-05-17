_author__ = "Anders Logg (logg@simula.no)"
__date__ = "2005-09-16 -- 2007-03-29"
__copyright__ = "Copyright (C) 2005-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Modified by Garth N. Wells 2006
# Modified by Marie E. Rognes (meg@math.uio.no) 2007

# Python modules
import numpy

# FFC common modules
from ffc.common.debug import *
from ffc.common.utils import *

# FFC language modules
from ffc.compiler.language.index import Index

class MixedElement:
    """A MixedElement represents a finite element defined as a tensor
    product of finite elements. It is represented as a list of finite
    elements (mixed or simple) and may thus be recursively defined in
    terms of other mixed elements."""
    
    def __init__(self, elements):
        "Create MixedElement from a list of elements."

        # Make sure we get a list of elements
        if not isinstance(elements, list) or not len(elements) > 1:
            raise FormError, "Mixed finite element must be created from a list of at least two elements."

        # Save list of elements
        self.__elements = elements

    def family(self):
        "Return a string indentifying the finite element family"
        return "Mixed"

    def signature(self):
        "Return a string identifying the finite element"
        return "Mixed finite element: [%s]" % ", ".join([element.signature() for element in self.__elements])

    def cell_shape(self):
        "Return the cell shape"
        return pick_first([element.cell_shape() for element in self.__elements])

    def space_dimension(self):
        "Return the dimension of the finite element function space"
        return sum([element.space_dimension() for element in self.__elements])

    def value_rank(self):
        "Return the rank of the value space"
        return 1

    def value_dimension(self, i):
        "Return the dimension of the value space for axis i"
        return sum([element.value_dimension(i) for element in self.__elements])

    def num_sub_elements(self):
        "Return the number of sub elements"
        return len(self.__elements)

    def sub_element(self, i):
        "Return sub element i"
        return self.__elements[i]

    def degree(self):
        "Return degree of polynomial basis"
        return max([element.degree() for element in self.__elements])

    def value_mapping(self, component):
        """Return the type of mapping associated with the given
        component of the element. """
        i = Index(component)
        if i.type == i.FIXED:
            (sub_element, offset) = self.value_offset(i)
            return sub_element.value_mapping(Index(i.index - offset))
        else:
            if component.range:
                return pick_first([self.value_mapping(i)
                                   for i in component.range])
            else:
                return pick_first([self.value_mapping(i)
                                   for i in range(self.value_dimension(0))])

    def space_mapping(self, i):
        """Return the type of mapping associated with the i'th basis
        function of the element"""
        (sub_element, offset) = self.space_offset(i)
        return sub_element.space_mapping(i - offset)
        
    def value_offset(self, component):
        """Given an absolute component (index), return the associated
        subelement and relative position of the component""" 
        # Does not yet work with nested mixed elements
        adjustment = 0
        for element in self.__elements:
            value_dim = element.value_dimension(0)
            if (adjustment + value_dim) > component.index:
                return (element, adjustment)
            else:
                adjustment += value_dim
        raise RuntimeError("Component does not match value dimension")

    def space_offset(self, i):
        """Given an absolute basis_no (i), return the associated
        subelement and offset"""
        # Does not yet work with nested mixed elements
        adjustment = 0
        for element in self.__elements:
            space_dim = element.space_dimension()
            if (adjustment + space_dim) > i:
                return (element, adjustment)
            else:
                adjustment += space_dim
        raise RuntimeError("Basis number does not match space dimension")

    def cell_dimension(self):
        "Return dimension of shape"
        return pick_first([element.cell_dimension() for element in self.__elements])

    def facet_shape(self):
        "Return shape of facet"
        return pick_first([element.facet_shape() for element in self.__elements])

    def num_facets(self):
        "Return number of facets for shape of element"
        return pick_first([element.num_facets() for element in self.__elements])

    def entity_dofs(self):
        """Return the mapping from entities to dofs. Note that we
        unnest the possibly recursively nested entity_dofs here to
        generate just a list of entity dofs for basic elements."""
        return [entity_dofs for element in self.__elements for entity_dofs in element.entity_dofs()]

    def basis(self):
        "Return basis of finite element space"
        raise RuntimeError, "Basis cannot be accessed explicitly for a mixed element."

    def tabulate(self, order, points, facet = None):
        """Tabulate values on mixed element by appropriately reordering
        the tabulated values for the sub elements."""

        # Special case: only one element
        if len(self.__elements) == 1:
            return elements[0].tabulate(order, points, facet)

        # Iterate over sub elements and build mixed table from element tables
        mixed_table = []
        offset = 0
        for i in range(len(self.__elements)):
            # Get current element and table
            element = self.__elements[i]
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

    def basis_elements(self):
        "Returns a list of all basis elements"
        return self.__extract_elements(self)

    def __compute_mixed_entity_dofs(self, elements):
        "Compute mixed entity dofs as a list of entity dof mappings"
        mixed_entity_dofs = []
        for element in elements:
            if isinstance(element.entity_dofs(), list):
                mixed_entity_dofs += element.entity_dofs()
            else:
                mixed_entity_dofs += [element.entity_dofs()]
        return mixed_entity_dofs

    def __compute_component_table(self, table, offset):
        "Compute subtable for given component"
        component_table = []
        # Iterate over derivative orders
        for dorder in range(len(table)):
            component_table.append({})
            # Iterate over derivative tuples
            derivative_dictionary = {}
            for dtuple in table[dorder]:
                element_subtable = table[dorder][dtuple]
                num_points = numpy.shape(element_subtable)[1]
                mixed_subtable = numpy.zeros((self.space_dimension(), num_points), dtype = numpy.float)
                # Iterate over element basis functions and fill in non-zero values
                for i in range(len(element_subtable)):
                    mixed_subtable[offset + i] = element_subtable[i]
                # Add to dictionary
                component_table[dorder][dtuple] = mixed_subtable
        return component_table

    def __extract_elements(self, element):
        """This function extracts the basis elements recursively from vector elements and mixed elements.
        Example, the following mixed element:

        element1 = FiniteElement("Lagrange", "triangle", 1)
        element2 = VectorElement("Lagrange", "triangle", 2)

        element  = element2 + element1, has the structure:
        mixed-element[mixed-element[Lagrange order 2, Lagrange order 2], Lagrange order 1]

        This function returns the list of basis elements:
        elements = [Lagrange order 2, Lagrange order 2, Lagrange order 1]"""

        elements = []

        # Import here to avoid cyclic dependency
        from finiteelement import FiniteElement

        # If the element is not mixed (a basis element, add to list)        
        if isinstance(element, FiniteElement):
            elements += [element]
        # Else call this function again for each subelement
        else:
            for i in range(element.num_sub_elements()):
                elements += self.__extract_elements(element.sub_element(i))

        return elements

    def __add__(self, other):
        "Create mixed element"
        return MixedElement([self, other])

    def __repr__(self):
        "Pretty print"
        return "Mixed finite element: " + str(self.__elements)
