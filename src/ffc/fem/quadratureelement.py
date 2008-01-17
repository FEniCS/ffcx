__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-12-10 -- 2007-01-16"
__copyright__ = "Copyright (C) 2007-2008 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# FFC fem modules
from finiteelement import *
from quadrature import *
from mapping import *

# We could remove this dummy class as explained later
class dummy_dual_basis:
    def __init__(self, points):
        self.pts = points

class QuadratureElement(FiniteElement):
    """Write description of QuadratureElement"""

    def __init__(self, shape, num_points_per_axis):
        "Create QuadratureElement"

        # Save incoming arguments
        self.__cell_shape = string_to_shape[shape]
        self.__num_axis_points = num_points_per_axis

        # Save element family
        self.__family = "Quadrature"

        # Set rank (is rank = 0 for this element?)
        self.__rank = 0

        # Save element degree (constant)
        self.__degree = 0

        # Set mapping to AFFINE (not important, I think, for this element)
        self.__mapping = Mapping.AFFINE

        # Create quadrature (only interested in points)
        points, weights = make_quadrature(self.__cell_shape, num_points_per_axis)

        # Save number of quadrature points
        self.__num_quad_points = len(points)

        # Create entity IDs, ripped from FIAT/DiscontinuousLagrange.py
        # Used by formdata.py to create the DofMap
        entity_ids = {}
        for d in range( shape_to_dim[ self.__cell_shape ] ):
            entity_ids[d] = {}
            for e in range(num_entities[self.__cell_shape][d]):
                entity_ids[d][e] = []
        entity_ids[ shape_to_dim[ self.__cell_shape ]] = {}
        entity_ids[ shape_to_dim[ self.__cell_shape ]][ 0 ] = range( self.__num_quad_points )
        self.__entity_dofs = [entity_ids]

        # FIXME: We could change calls to element.dual_basis().pts to
        # element.dof_coordinates() which would only support a limited number of
        # elements.
        # Initialise a dummy dual_basis.
        # Used in dofmap.py", line 158, in __compute_dof_coordinates
        self.__dummy_dual_basis = dummy_dual_basis(points)

    def family(self):
        "Return a string indentifying the finite element family"
        return self.__family

    def signature(self):
        "Return a string identifying the QuadratureElement"
        s = ""
        if self.cell_shape() > 1:
            s = " (%s)" % "x".join([str(self.__num_axis_points) for i in range(self.__cell_shape)])

        return "%s element with %d quadrature point(s)%s on a %s" % \
               (self.__family, self.__num_quad_points, s, shape_to_string[self.__cell_shape])

    def cell_shape(self):
        "Return the cell shape"
        return self.__cell_shape

    def space_dimension(self):
        "Return the total number of quadrature points"
        return self.__num_quad_points

    def value_rank(self):
        "Return the rank of the value space"
        return self.__rank

#     def value_dimension(self, i):
#         "Return the dimension of the value space for axis i"
#         if self.value_rank() == 0:
#             return 1
#         else:
#             return self.basis().tensor_dim()[i]

#     def num_sub_elements(self):
#         "Return the number of sub elements"
#         return 1

#     def sub_element(self, i):
#         "Return sub element i"
#         return self

    def degree(self):
        "Return degree of polynomial basis"
        return self.__degree

    def value_mapping(self, component):
        """Return the type of mapping associated with the i'th
        component of the element"""
        return self.__mapping

    def space_mapping(self, i):
        """Return the type of mapping associated with the i'th basis
        function of the element"""
        return self.__mapping

#     def value_offset(self, component):
#         """Given an absolute component (index), return the associated
#         subelement and offset of the component""" 
#         return (self, 0)

#     def space_offset(self, i):
#         """Given a basis function number i, return the associated
#         subelement and offset""" 
#         return (self, 0)
#     
#     def cell_dimension(self):
#         "Return dimension of shape"
#         return shape_to_dim[self.cell_shape()]

#     def facet_shape(self):
#         "Return shape of facet"
#         return shape_to_facet[self.cell_shape()]

#     def num_facets(self):
#         "Return number of facets for shape of element"
#         return shape_to_num_facets[self.cell_shape()]

    def entity_dofs(self):
        "Return the mapping from entities to dofs"
        return self.__entity_dofs

#     def basis(self):
#         "Return basis of finite element space"
#         return self.__fiat_element.function_space()

    # FIXME: Possibly remove this function
    def dual_basis(self):
        "Return dummy dual basis of finite element space"
        # This function only makes sense in a call like: element.dual_basis().pts
        return self.__dummy_dual_basis

    def tabulate(self, order, points):
        """Return the identity matrix of size (num_quad_points, num_quad_points),
        in a format that monomialintegration and monomialtabulation understands."""

        # Derivatives are not defined on a QuadratureElement
        if order:
            raise RuntimeError("Derivatives are not defined on a QuadratureElement")

        # Check if incoming points are equal to quadrature points. (If they are,
        # then monomialintegration or monomialtabulation is probably the caller.)
        # This 'direct' check is very dangerous, a small change in the values
        # from make_quadrature() will result in disaster!!
        if points == self.dual_basis().pts:
            values = numpy.identity(self.__num_quad_points, float)
            table = [{(0,)*self.__cell_shape: values}]
            return table
        else:
            raise RuntimeError("points must be equal to coordinates of quadrature points \n %s \n %s" %(points, self.dual_basis().pts))

#     def basis_elements(self):
#         "Returns a list of all basis elements"
#         return [self]

    def num_axis_points(self):
        "Return the number of quadrature points per axis as specified by user"
        return self.__num_axis_points

#     def __add__(self, other):
#         "Create mixed element"
#         return mixedelement.MixedElement([self, other])

    def __repr__(self):
        "Pretty print"
        return self.signature()

