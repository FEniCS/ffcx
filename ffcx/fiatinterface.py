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


def reference_cell_vertices(cellname):
    """Return dict of coordinates of reference cell vertices for this 'cellname'."""
    return FIAT.ufc_cell(cellname).get_vertices()


@functools.singledispatch
def _create_element(element):
    raise ValueError("Element type is not supported.")


@_create_element.register(ufl.FiniteElement)
def _create_finiteelement(element: ufl.FiniteElement) -> FIAT.FiniteElement:
    """Create FIAT element for UFL base type ufl.FiniteElement."""
    if element.family() == "Real":
        e = create_element(ufl.FiniteElement("DG", element.cell(), 0))
        e.__class__ = type('SpaceOfReals', (type(e), SpaceOfReals), {})
        return e

    if element.family() == "Quadrature":
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
        raise ValueError("Finite element of type \"{}\" is not supported by FIAT.".format(element.family()))

    # Handle Lagrange variants
    if element.family() == "Lagrange" and element.variant() == "spectral":
        assert element.cell().cellname() == "interval"
        element_class = FIAT.GaussLobattoLegendre
    else:
        element_class = FIAT.supported_elements[element.family()]

    assert element.degree() is not None
    return element_class(FIAT.ufc_cell(element.cell().cellname()), element.degree())


@_create_element.register(ufl.MixedElement)
def _create_mixed_finiteelement(element: ufl.MixedElement) -> FIAT.MixedElement:
    elements = []

    def rextract(els):
        for e in els:
            if isinstance(e, ufl.MixedElement):
                rextract(e.sub_elements())
            else:
                elements.append(e)

    rextract(element.sub_elements())
    return FIAT.MixedElement(map(_create_element, elements))


@_create_element.register(ufl.EnrichedElement)
def _create_enriched_finiteelement(element: ufl.EnrichedElement) -> FIAT.EnrichedElement:
    elements = [create_element(e) for e in element._elements]
    return FIAT.EnrichedElement(*elements)


@_create_element.register(ufl.NodalEnrichedElement)
def _create_nodalenriched_finiteelement(element: ufl.NodalEnrichedElement) -> FIAT.NodalEnrichedElement:
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
    print(ufl_element, type(ufl_element))
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

    quad_rule = FIAT.create_quadrature(FIAT.ufc_cell(shape), degree, scheme)
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
    fiat_cell = FIAT.ufc_cell(cellname)

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
