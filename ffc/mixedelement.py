_author__ = "Anders Logg (logg@simula.no)"
__date__ = "2005-09-16"
__copyright__ = "Copyright (C) 2005-2007 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Garth N. Wells 2006-2009
# Modified by Marie E. Rognes (meg@math.uio.no) 2007
# Modified by Kristian B. Oelgaard 2009
# Last changed: 2009-12-16

# Python modules.
import numpy

# UFL modules.
from ufl.classes import FiniteElementBase

# FFC modules.
from log import error
from dofrepresentation import DofRepresentation

class MixedElement(FiniteElementBase):
    """A MixedElement represents a finite element defined as a tensor
    product of finite elements. It is represented as a list of finite
    elements (mixed or simple) and may thus be recursively defined in
    terms of other mixed elements."""

    # TODO: KBO: change the argument list to get the ufl_element, modify
    # create_element accordingly
    def __init__(self, elements, ufl_str, value_shape, domain=None):
        "Create MixedElement from a list of elements."

        # Initialise base class
        degree = max(e.degree() for e in elements)
        FiniteElementBase.__init__(self, "Mixed", elements[0].cell(),
                                   degree, value_shape)
        # Save UFL string representation
        self._repr = ufl_str

        # Make sure we get a list of elements
        if not isinstance(elements, list):
            error(elements, "Mixed finite element must be created from a list of at least two elements.")

        # Save list of elements
        self._elements = elements

        # Save domain, no need to do any checks, that should be handled by the subelements
        self._domain = domain

        # FIXME: This is just a temporary hack to 'support' tensor elements
        self._rank = 1

    def __add__(self, other):
        "Create mixed element"
        return MixedElement([self, other])

    def basis(self):
        "Return basis of finite element space"
        error("Basis cannot be accessed explicitly for a mixed element.")

    def component_element(self, component):
        "Return sub element and offset for given component."
        offset = 0
        for element in self.extract_elements():
            next_offset = offset + element.value_dimension(0)
            if next_offset > component:
                return (element, offset)
            offset = next_offset
        error("Unable to extract sub element for component %s of %s." % (str(component), str(self)))

    def dual_basis(self):
        """Return the representation of the dofs. We unnest the
        possibly nested dof types as for entity_dofs and shift the
        components according to the position of the basic elements."""

        # Calculate the shifts in components caused by the positioning
        # of the elements:
        dims = [element.value_dimension(0) for element in self._elements]
        shifts = [sum(dims[:i]) for i in range(len(dims))]
        n = sum(dims)

        # Shift the components for the separate dofs
        dofs = []
        for e in range(len(self._elements)):
            element = self._elements[e]
            shift = shifts[e]
            for d in element.dual_basis():
                dof = DofRepresentation(d)
                dof.shift_directions(n, shift)
                dofs += [dof]
        return dofs

    def entity_dofs(self):
        """Return the mapping from entities to dofs. Note that we
        unnest the possibly recursively nested entity_dofs here to
        generate just a list of entity dofs for basic elements."""
        return [entity_dofs for element in self._elements for entity_dofs in element.entity_dofs()]

    def extract_elements(self):
        "Extract list of all recursively nested elements."
        return _extract_elements(self)

    def num_sub_elements(self):
        "Return the number of sub elements"
        return len(self._elements)

    def space_dimension(self):
        "Return the dimension of the finite element function space"
        return sum([element.space_dimension() for element in self._elements])

    def sub_element(self, i):
        "Return sub element i"
        return self._elements[i]

    def tabulate(self, order, points):
        """Tabulate values on mixed element by appropriately reordering
        the tabulated values for the sub elements."""

        # Special case: only one element
        if len(self._elements) == 1:
            return elements[0].tabulate(order, points)

        # Iterate over sub elements and build mixed table from element tables
        mixed_table = []
        offset = 0
        for i in range(len(self._elements)):
            # Get current element and table
            element = self._elements[i]
            table = element.tabulate(order, points)
            # Iterate over the components corresponding to the current element
            if element.value_rank() == 0:
                component_table = _compute_component_table(table, offset, self.space_dimension())
                mixed_table.append(component_table)
            else:
                for i in range(element.value_dimension(0)):
                    component_table = _compute_component_table(table[i], offset, self.space_dimension())
                    mixed_table.append(component_table)
            # Add to offset, the number of the first basis function for the current element
            offset += element.space_dimension()

        return mixed_table

    def value_dimension(self, i):
        "Return the dimension of the value space for axis i"
        return sum([element.value_dimension(i) for element in self._elements])

    def value_rank(self):
        "Return the rank of the value space"
        return 1

def _compute_component_table(table, offset, space_dimension):

    print "TABLE"
    print table
    print "END TABLE"

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
            mixed_subtable = numpy.zeros((space_dimension, num_points), dtype = numpy.float)
            # Iterate over element basis functions and fill in non-zero values
            for i in range(len(element_subtable)):
                print "CHECK --------------"
                print "CHECK", element_subtable[i]
                print "CHECK", mixed_subtable[offset + i]
                mixed_subtable[offset + i] = element_subtable[i]
            # Add to dictionary
            component_table[dorder][dtuple] = mixed_subtable
    return component_table

def _compute_mixed_entity_dofs(elements):
    "Compute mixed entity dofs as a list of entity dof mappings"
    mixed_entity_dofs = []
    for element in elements:
        if isinstance(element.entity_dofs(), list):
            mixed_entity_dofs += element.entity_dofs()
        else:
            mixed_entity_dofs += [element.entity_dofs()]
    return mixed_entity_dofs

# def _extract_elements(element):
#     """This function extracts the basis elements recursively from vector elements and mixed elements.
#     Example, the following mixed element:

#     element1 = FiniteElement("Lagrange", "triangle", 1)
#     element2 = VectorElement("Lagrange", "triangle", 2)

#     element  = element2 + element1, has the structure:
#     mixed-element[mixed-element[Lagrange order 2, Lagrange order 2], Lagrange order 1]

#     This function returns the list of basis elements:
#     elements = [Lagrange order 2, Lagrange order 2, Lagrange order 1]"""

#     elements = []

#     # Import here to avoid cyclic dependency
#     #from finiteelement import FiniteElement

#     # If the element is not mixed (a basis element, add to list)
#     if isinstance(element, FiniteElement):
#         elements += [element]
#     # Else call this function again for each subelement
#     else:
#         for i in range(element.num_sub_elements()):
#             elements += _extract_elements(element.sub_element(i))

#     return elements
