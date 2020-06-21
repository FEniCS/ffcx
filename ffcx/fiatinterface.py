# Copyright (C) 2009-2017 Kristian B. Oelgaard and Anders Logg
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import functools
import logging

import numpy

import FIAT
import ufl

logger = logging.getLogger("ffcx")

# Element families supported by FFCX
supported_families = ("Brezzi-Douglas-Marini", "Brezzi-Douglas-Fortin-Marini", "Crouzeix-Raviart",
                      "Discontinuous Lagrange", "Discontinuous Raviart-Thomas", "HDiv Trace",
                      "Lagrange", "Lobatto", "Nedelec 1st kind H(curl)", "Nedelec 2nd kind H(curl)",
                      "Radau", "Raviart-Thomas", "Real", "Bubble", "Quadrature", "Regge",
                      "Hellan-Herrmann-Johnson", "Q", "DQ", "TensorProductElement", "Gauss-Lobatto-Legendre")

# Cache for computed elements
_cache = {}

_tpc_quadrilateral = ufl.TensorProductCell(ufl.interval, ufl.interval)
_tpc_hexahedron = ufl.TensorProductCell(ufl.quadrilateral, ufl.interval)


class SpaceOfReals(object):
    """Constant over the entire domain, rather than just cellwise."""


def reference_cell(cellname):
    """Return FIAT reference cell."""
    return FIAT.ufc_cell(cellname)


def reference_cell_vertices(cellname):
    """Return dict of coordinates of reference cell vertices for this 'cellname'."""
    cell = reference_cell(cellname)
    return cell.get_vertices()


@functools.singledispatch
def _create_element(element):
    raise ValueError("Unsupported element")


@_create_element.register(ufl.FiniteElement)
def _create_finiteelement(element: ufl.FiniteElement) -> FIAT.FiniteElement:
    """Base FiniteElement"""
    if element.family() == "Real":
        e = create_element(ufl.FiniteElement("DG", element.cell(), 0))
        e.__class__ = type('SpaceOfReals', (type(e), SpaceOfReals), {})
        return e

    if element.family() == "Quadrature":
        # Compute number of points per axis from the degree of the element
        # scheme = element.quadrature_scheme()
        assert element.degree() is not None
        assert element.quadrature_scheme() is not None

        # Create quadrature (only interested in points)
        # TODO: KBO: What should we do about quadrature functions that live on ds, dS?
        # Get cell and facet names.
        points, weights = create_quadrature(element.cell().cellname(), element.degree(),
                                            element.quadrature_scheme())
        return FIAT.QuadratureElement(FIAT.ufc_cell(element.cell().cellname()), points)

    # Handle tensor-product structured elements
    if element.cell().cellname() == "quadrilateral":
        e = _create_element(element.reconstruct(cell=_tpc_quadrilateral))
        return FIAT.tensor_product.FlattenedDimensions(e)
    elif element.cell().cellname() == "hexahedron":
        e = _create_element(element.reconstruct(cell=_tpc_hexahedron))
        return FIAT.tensor_product.FlattenedDimensions(e)

    if element.family() not in FIAT.supported_elements:
        raise ValueError("Finite element of type \"{}\" is not supported by FIAT.".format(family))

    element_class = FIAT.supported_elements[element.family()]
    assert element.degree() is not None
    return element_class(FIAT.ufc_cell(element.cell().cellname()), element.degree())


@_create_element.register(ufl.MixedElement)
def _create_mixed_finiteelement(element: ufl.MixedElement) -> FIAT.MixedElement:
    elements = _extract_elements(element)
    return FIAT.MixedElement(elements)


@_create_element.register(ufl.EnrichedElement)
def _create_enriched_finiteelement(element: ufl.EnrichedElement) -> FIAT.EnrichedElement:
    elements = [create_element(e) for e in element._elements]
    return FIAT.EnrichedElement(*elements)


@_create_element.register(ufl.NodalEnrichedElement)
def _create_nodelenriched_finiteelement(element: ufl.NodalEnrichedElement) -> FIAT.NodalEnrichedElement:
    elements = [create_element(e) for e in element._elements]
    return FIAT.NodalEnrichedElement(*elements)


@_create_element.register(ufl.RestrictedElement)
def _create_restricted_finiteelement(element: ufl.RestrictedElement):
    # element = _create_restricted_element(element)
    raise RuntimeError("Cannot handle this element type: {}".format(element))


@_create_element.register(ufl.TensorProductElement)
def _create_tp_finiteelement(element) -> FIAT.TensorProductElement:
    e0, e1 = element.sub_elements()
    return FIAT.TensorProductElement(_create_element(e0), _create_element(e1))


def create_element(ufl_element: ufl.finiteelement) -> FIAT.FiniteElement:
    """Create a FIAT finite element for a given UFL element."""

    # Use UFL element as cache key
    element_signature = ufl_element
    if element_signature in _cache:
        return _cache[element_signature]

    # Create element and add to cache
    element = _create_element(ufl_element)
    _cache[element_signature] = element

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
        element = _create_element(base_element)
        return FIAT.RestrictedElement(element, restriction_domain=restriction_domain)

    # If restricted mixed element -> convert to mixed restricted element
    if isinstance(base_element, ufl.MixedElement):
        elements = _extract_elements(base_element, restriction_domain)
        return FIAT.MixedElement(elements)

    raise RuntimeError("Cannot create restricted element from: {}".format(ufl_element))
