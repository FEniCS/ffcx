# Copyright (C) 2009-2017 Kristian B. Oelgaard and Anders Logg
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import logging
import warnings

import numpy

import FIAT
import ufl
from FIAT.enriched import EnrichedElement
from FIAT.mixed import MixedElement
from FIAT.nodal_enriched import NodalEnrichedElement
from FIAT.quadrature_element import QuadratureElement
from FIAT.restricted import RestrictedElement
from FIAT.tensor_product import FlattenedDimensions

logger = logging.getLogger(__name__)

# Element families supported by FFCX
supported_families = ("Brezzi-Douglas-Marini", "Brezzi-Douglas-Fortin-Marini", "Crouzeix-Raviart",
                      "Discontinuous Lagrange", "Discontinuous Raviart-Thomas", "HDiv Trace",
                      "Lagrange", "Lobatto", "Nedelec 1st kind H(curl)", "Nedelec 2nd kind H(curl)",
                      "Radau", "Raviart-Thomas", "Real", "Bubble", "Quadrature", "Regge",
                      "Hellan-Herrmann-Johnson", "Q", "DQ", "TensorProductElement")

# Cache for computed elements
_cache = {}


class SpaceOfReals(object):
    """Constant over the entire domain, rather than just cellwise."""


def reference_cell(cellname):
    """Return FIAT reference cell."""
    return FIAT.ufc_cell(cellname)


def reference_cell_vertices(cellname):
    """Return dict of coordinates of reference cell vertices for this 'cellname'."""
    cell = reference_cell(cellname)
    return cell.get_vertices()


def create_element(ufl_element):

    # Create element signature for caching (just use UFL element)
    element_signature = ufl_element

    # Check cache
    if element_signature in _cache:
        logger.debug("Reusing element from cache")
        return _cache[element_signature]

    if isinstance(ufl_element, ufl.FiniteElement):
        element = _create_fiat_element(ufl_element)
    elif isinstance(ufl_element, ufl.MixedElement):
        elements = _extract_elements(ufl_element)
        element = MixedElement(elements)
    elif isinstance(ufl_element, ufl.EnrichedElement):
        elements = [create_element(e) for e in ufl_element._elements]
        element = EnrichedElement(*elements)
    elif isinstance(ufl_element, ufl.NodalEnrichedElement):
        elements = [create_element(e) for e in ufl_element._elements]
        element = NodalEnrichedElement(*elements)
    elif isinstance(ufl_element, ufl.RestrictedElement):
        element = _create_restricted_element(ufl_element)
        raise RuntimeError("Cannot handle this element type: {}".format(ufl_element))

    # Store in cache
    _cache[element_signature] = element

    return element


def _create_fiat_element(ufl_element):
    """Create FIAT element corresponding to given finite element."""

    # Get element data
    family = ufl_element.family()
    cell = ufl_element.cell()
    cellname = cell.cellname()
    degree = ufl_element.degree()

    # Check that FFCX supports this element
    if family not in supported_families:
        raise RuntimeError("This element family (%s) is not supported by FFC." % family)

    # Create FIAT cell
    fiat_cell = reference_cell(cellname)

    # Handle the space of the constant
    if family == "Real":
        element = _create_fiat_element(ufl.FiniteElement("DG", cell, 0))
        element.__class__ = type('SpaceOfReals', (type(element), SpaceOfReals), {})
        return element

    if cellname == "quadrilateral":
        # Handle quadrilateral case by reconstructing the element with
        # cell TensorProductCell (interval x interval)
        quadrilateral_tpc = ufl.TensorProductCell(ufl.Cell("interval"), ufl.Cell("interval"))
        return FlattenedDimensions(
            _create_fiat_element(ufl_element.reconstruct(cell=quadrilateral_tpc)))
    elif cellname == "hexahedron":
        # Handle hexahedron case by reconstructing the element with cell
        # TensorProductCell (quadrilateral x interval). This creates
        # TensorProductElement(TensorProductElement(interval, interval),
        # interval) Therefore dof entities consists of nested tuples,
        # example: ((0, 1), 1)
        hexahedron_tpc = ufl.TensorProductCell(ufl.Cell("quadrilateral"), ufl.Cell("interval"))
        return FlattenedDimensions(
            _create_fiat_element(ufl_element.reconstruct(cell=hexahedron_tpc)))

    # FIXME: AL: Should this really be here?
    # Handle QuadratureElement
    if family == "Quadrature":
        # Compute number of points per axis from the degree of the element
        scheme = ufl_element.quadrature_scheme()
        assert degree is not None
        assert scheme is not None

        # Create quadrature (only interested in points)
        # TODO: KBO: What should we do about quadrature functions that live on ds, dS?
        # Get cell and facet names.
        points, weights = create_quadrature(cellname, degree, scheme)

        # Make element
        element = QuadratureElement(fiat_cell, points)
    else:
        # Check if finite element family is supported by FIAT
        if family not in FIAT.supported_elements:
            raise RuntimeError("Sorry, finite element of type \"%s\" are not supported by FIAT.",
                               family)

        ElementClass = FIAT.supported_elements[family]

        # Create tensor product FIAT finite element
        if isinstance(ufl_element, ufl.TensorProductElement):
            A = create_element(ufl_element.sub_elements()[0])
            B = create_element(ufl_element.sub_elements()[1])
            element = ElementClass(A, B)

        # Create normal FIAT finite element
        else:
            if degree is None:
                element = ElementClass(fiat_cell)
            else:
                element = ElementClass(fiat_cell, degree)

    if element.value_shape() != ufl_element.reference_value_shape():
        # Consistency check between UFL and FIAT elements.
        raise RuntimeError("Something went wrong in the construction of FIAT element from UFL element."
                           + "Shapes are {} and {}.".format(element.value_shape(),
                                                            ufl_element.reference_value_shape()))

    return element


def create_quadrature(shape, degree, scheme="default"):
    """Generate quadrature rule.

    Quadrature rule(points, weights) for given shape
    that will integrate an polynomial of order 'degree' exactly.

    """
    if isinstance(shape, int) and shape == 0:
        return (numpy.zeros((1, 0)), numpy.ones((1, )))

    if shape in ufl.cell.cellname2dim and ufl.cell.cellname2dim[shape] == 0:
        return (numpy.zeros((1, 0)), numpy.ones((1, )))

    if scheme == "vertex":
        # The vertex scheme, i.e., averaging the function value in the
        # vertices and multiplying with the simplex volume, is only of
        # order 1 and inferior to other generic schemes in terms of
        # error reduction. Equation systems generated with the vertex
        # scheme have some properties that other schemes lack, e.g., the
        # mass matrix is a simple diagonal matrix. This may be
        # prescribed in certain cases.
        if degree > 1:
            warnings.warn(
                "Explicitly selected vertex quadrature (degree 1), but requested degree is {}.".
                format(degree))
        if shape == "tetrahedron":
            return (numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0,
                                                                                     1.0]]),
                    numpy.array([1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0]))
        elif shape == "triangle":
            return (numpy.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
                    numpy.array([1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0]))
        elif shape == "interval":
            # Trapezoidal rule.
            return (numpy.array([[0.0], [1.0]]), numpy.array([1.0 / 2.0, 1.0 / 2.0]))

    quad_rule = FIAT.create_quadrature(reference_cell(shape), degree, scheme)
    points = numpy.asarray(quad_rule.get_points())
    weights = numpy.asarray(quad_rule.get_weights())
    return points, weights


def map_facet_points(points, facet, cellname):
    """Map points from a facet to a cell.

    Map points from the e (UFC) reference simplex of dimension d - 1
    to a given facet on the (UFC) reference simplex of dimension d. This
    may be used to transform points tabulated for example on the 2D
    reference triangle to points on a given facet of the reference
    tetrahedron.

    """

    # Extract the geometric dimension of the points we want to map
    dim = len(points[0]) + 1

    # Special case, don't need to map coordinates on vertices
    if dim == 1:
        return [[(0.0, ), (1.0, )][facet]]

    # Get the FIAT reference cell
    fiat_cell = reference_cell(cellname)

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
        w = (1.0 - sum(point), ) + tuple(point)
        x = tuple(sum([w[i] * numpy.array(coordinates[i]) for i in range(len(w))]))
        new_points += [x]

    return new_points


def _extract_elements(ufl_element, restriction_domain=None):
    """Recursively extract un-nested list of (component) elements."""

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
    """Create an FFCX representation for an UFL RestrictedElement."""

    if not isinstance(ufl_element, ufl.RestrictedElement):
        raise RuntimeError("create_restricted_element expects an ufl.RestrictedElement")

    base_element = ufl_element.sub_element()
    restriction_domain = ufl_element.restriction_domain()

    # If simple element -> create RestrictedElement from fiat_element
    if isinstance(base_element, ufl.FiniteElement):
        element = _create_fiat_element(base_element)
        return RestrictedElement(element, restriction_domain=restriction_domain)

    # If restricted mixed element -> convert to mixed restricted element
    if isinstance(base_element, ufl.MixedElement):
        elements = _extract_elements(base_element, restriction_domain)
        return MixedElement(elements)

    raise RuntimeError("Cannot create restricted element from: {}".format(ufl_element))
