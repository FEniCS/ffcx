# -*- coding: utf-8 -*-

# Copyright (C) 2009-2016 Kristian B. Oelgaard and Anders Logg
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
# Modified by Garth N. Wells, 2009.
# Modified by Marie Rognes, 2009-2013.
# Modified by Martin Sandve Alnæs, 2013
# Modified by Lizao Li, 2015, 2016

# Python modules
import six
import numpy
from numpy import array

# UFL and FIAT modules
import ufl
from ufl.utils.sorting import sorted_by_key
import FIAT
from FIAT.hdiv_trace import HDivTrace

# FFC modules
from ffc.log import debug, error
from ffc.mixedelement import MixedElement
from ffc.restrictedelement import RestrictedElement
from ffc.enrichedelement import EnrichedElement, SpaceOfReals

# Dictionary mapping from cellname to dimension
from ufl.cell import cellname2dim

# Element families supported by FFC
supported_families = ("Brezzi-Douglas-Marini",
                      "Brezzi-Douglas-Fortin-Marini",
                      "Crouzeix-Raviart",
                      "Discontinuous Lagrange",
                      "Discontinuous Raviart-Thomas",
                      "HDiv Trace",
                      "Lagrange",
                      "Lobatto",
                      "Nedelec 1st kind H(curl)",
                      "Nedelec 2nd kind H(curl)",
                      "Radau",
                      "Raviart-Thomas",
                      "Real",
                      "Bubble",
                      "Quadrature",
                      "Regge",
                      "Hellan-Herrmann-Johnson")

# Cache for computed elements
_cache = {}


def reference_cell(dim):
    if isinstance(dim, int):
        return FIAT.ufc_simplex(dim)
    else:
        return FIAT.ufc_simplex(cellname2dim[dim])


def reference_cell_vertices(dim):
    "Return dict of coordinates of reference cell vertices for this 'dim'."
    cell = reference_cell(dim)
    return cell.get_vertices()


def create_element(ufl_element):

    # Create element signature for caching (just use UFL element)
    element_signature = ufl_element

    # Check cache
    if element_signature in _cache:
        debug("Reusing element from cache")
        return _cache[element_signature]

    # Create regular FIAT finite element
    if isinstance(ufl_element, ufl.FiniteElement):
        element = _create_fiat_element(ufl_element)

    # Create mixed element (implemented by FFC)
    elif isinstance(ufl_element, ufl.MixedElement):
        elements = _extract_elements(ufl_element)
        element = MixedElement(elements)

    # Create element union (implemented by FFC)
    elif isinstance(ufl_element, ufl.EnrichedElement):
        elements = [create_element(e) for e in ufl_element._elements]
        element = EnrichedElement(elements)

    # Create restricted element(implemented by FFC)
    elif isinstance(ufl_element, ufl.RestrictedElement):
        element = _create_restricted_element(ufl_element)

    else:
        error("Cannot handle this element type: %s" % str(ufl_element))

    # Store in cache
    _cache[element_signature] = element

    return element


def _create_fiat_element(ufl_element):
    "Create FIAT element corresponding to given finite element."

    # Get element data
    family = ufl_element.family()
    cell = ufl_element.cell()
    cellname = cell.cellname()
    degree = ufl_element.degree()

    # Check that FFC supports this element
    if family not in supported_families:
        error("This element family (%s) is not supported by FFC." % family)

    # Handle the space of the constant
    if family == "Real":
        dg0_element = ufl.FiniteElement("DG", cell, 0)
        constant = _create_fiat_element(dg0_element)
        element = SpaceOfReals(constant)

    # FIXME: AL: Should this really be here?
    # Handle QuadratureElement
    elif family == "Quadrature":
        element = QuadratureElement(ufl_element)

    else:
        # Create FIAT cell
        fiat_cell = reference_cell(cellname)

        # Handle Bubble element as RestrictedElement of P_{k} to interior
        if family == "Bubble":
            V = FIAT.supported_elements["Lagrange"](fiat_cell, degree)
            tdim = cell.topological_dimension()
            return RestrictedElement(V, _indices(V, "interior", tdim), None)

        # Check if finite element family is supported by FIAT
        if family not in FIAT.supported_elements:
            error("Sorry, finite element of type \"%s\" are not supported by FIAT.", family)

        # Create FIAT finite element
        ElementClass = FIAT.supported_elements[family]
        if degree is None:
            element = ElementClass(fiat_cell)
        else:
            element = ElementClass(fiat_cell, degree)

    # Consistency check between UFL and FIAT elements.
    if element.value_shape() != ufl_element.reference_value_shape():
        error("Something went wrong in the construction of FIAT element from UFL element." +
              "Shapes are %s and %s." % (element.value_shape(), ufl_element.reference_value_shape()))

    return element


def create_quadrature(shape, degree, scheme="default"):
    """
    Generate quadrature rule (points, weights) for given shape
    that will integrate an polynomial of order 'degree' exactly.
    """
    if isinstance(shape, int) and shape == 0:
        return (numpy.zeros((1, 0)), numpy.ones((1,)))

    if shape in cellname2dim and cellname2dim[shape] == 0:
        return (numpy.zeros((1, 0)), numpy.ones((1,)))

    if scheme == "vertex":
        # The vertex scheme, i.e., averaging the function value in the vertices
        # and multiplying with the simplex volume, is only of order 1 and
        # inferior to other generic schemes in terms of error reduction.
        # Equation systems generated with the vertex scheme have some
        # properties that other schemes lack, e.g., the mass matrix is
        # a simple diagonal matrix. This may be prescribed in certain cases.
        if degree > 1:
            from warnings import warn
            warn(("Explicitly selected vertex quadrature (degree 1), "
                 + "but requested degree is %d.") % degree)
        if shape == "tetrahedron":
            return (array([[0.0, 0.0, 0.0],
                           [1.0, 0.0, 0.0],
                           [0.0, 1.0, 0.0],
                           [0.0, 0.0, 1.0]]),
                    array([1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0])
                    )
        elif shape == "triangle":
            return (array([[0.0, 0.0],
                           [1.0, 0.0],
                           [0.0, 1.0]]),
                    array([1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0])
                    )
        elif shape == "interval":
            # Trapezoidal rule.
            return (array([[0.0],
                           [1.0]]),
                    array([1.0 / 2.0, 1.0 / 2.0])
                    )

    quad_rule = FIAT.create_quadrature(reference_cell(shape), degree, scheme)
    points = numpy.asarray(quad_rule.get_points())
    weights = numpy.asarray(quad_rule.get_weights())
    return points, weights


def map_facet_points(points, facet):
    """
    Map points from the e (UFC) reference simplex of dimension d - 1
    to a given facet on the (UFC) reference simplex of dimension d.
    This may be used to transform points tabulated for example on the
    2D reference triangle to points on a given facet of the reference
    tetrahedron.
    """

    # Extract the geometric dimension of the points we want to map
    dim = len(points[0]) + 1

    # Special case, don't need to map coordinates on vertices
    if dim == 1:
        return [[(0.0,), (1.0,)][facet]]

    # Get the FIAT reference cell for this dimension
    fiat_cell = reference_cell(dim)

    # Extract vertex coordinates from cell and map of facet index to
    # indicent vertex indices
    coordinate_dofs = fiat_cell.get_vertices()
    facet_vertices = fiat_cell.get_topology()[dim - 1]

    # coordinate_dofs = \
    #    {1: ((0.,), (1.,)),
    #     2: ((0., 0.), (1., 0.), (0., 1.)),
    #     3: ((0., 0., 0.), (1., 0., 0.),(0., 1., 0.), (0., 0., 1))}

    # Facet vertices
    # facet_vertices = \
    #    {2: ((1, 2), (0, 2), (0, 1)),
    #     3: ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2))}

    # Compute coordinates and map the points
    coordinates = [coordinate_dofs[v] for v in facet_vertices[facet]]
    new_points = []
    for point in points:
        w = (1.0 - sum(point),) + tuple(point)
        x = tuple(sum([w[i] * array(coordinates[i]) for i in range(len(w))]))
        new_points += [x]

    return new_points


def _extract_elements(ufl_element, restriction_domain=None):
    "Recursively extract un-nested list of (component) elements."

    elements = []
    if isinstance(ufl_element, ufl.MixedElement):
        for sub_element in ufl_element.sub_elements():
            elements += _extract_elements(sub_element, restriction_domain)
        return elements

    # Handle restricted elements since they might be mixed elements too.
    if isinstance(ufl_element, ufl.RestrictedElement):
        base_element = ufl_element.sub_element()
        restriction_domain = ufl_element.restriction_domain()
        return _extract_elements(base_element, restriction_domain)

    if restriction_domain:
        ufl_element = ufl.RestrictedElement(ufl_element, restriction_domain)

    elements += [create_element(ufl_element)]

    return elements


def _create_restricted_element(ufl_element):
    "Create an FFC representation for an UFL RestrictedElement."

    if not isinstance(ufl_element, ufl.RestrictedElement):
        error("create_restricted_element expects an ufl.RestrictedElement")

    base_element = ufl_element.sub_element()
    restriction_domain = ufl_element.restriction_domain()

    # If simple element -> create RestrictedElement from fiat_element
    if isinstance(base_element, ufl.FiniteElement):
        element = _create_fiat_element(base_element)
        tdim = ufl_element.cell().topological_dimension()
        return RestrictedElement(element, _indices(element, restriction_domain, tdim), restriction_domain)

    # If restricted mixed element -> convert to mixed restricted element
    if isinstance(base_element, ufl.MixedElement):
        elements = _extract_elements(base_element, restriction_domain)
        return MixedElement(elements)

    error("Cannot create restricted element from %s" % str(ufl_element))


def _indices(element, restriction_domain, tdim):
    "Extract basis functions indices that correspond to restriction_domain."

    # FIXME: The restriction_domain argument in FFC/UFL needs to be re-thought and
    # cleaned-up.

    # If restriction_domain is "interior", pick basis functions associated with
    # cell.
    if restriction_domain == "interior":
        return element.entity_dofs()[tdim][0]

    # Pick basis functions associated with
    # the topological degree of the restriction_domain and of all lower
    # dimensions.
    if restriction_domain == "facet":
        rdim = tdim - 1
    elif restriction_domain == "face":
        rdim = 2
    elif restriction_domain == "edge":
        rdim = 1
    elif restriction_domain == "vertex":
        rdim = 0
    else:
        error("Restriction to domain: %s, is not supported." % repr(restriction_domain))

    entity_dofs = element.entity_dofs()
    indices = []
    for dim in range(rdim + 1):
        entities = entity_dofs[dim]
        for (entity, index) in sorted_by_key(entities):
            indices += index
    return indices


# Import FFC module with circular dependency
from ffc.quadratureelement import QuadratureElement
