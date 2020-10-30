# Copyright (C) 2009-2020 Kristian B. Oelgaard, Anders Logg and Garth N. Wells
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import functools
import logging
import types

import numpy

import libtab
import FIAT
import ufl

logger = logging.getLogger("ffcx")

# Element families supported by FFCX
supported_families = ("Brezzi-Douglas-Marini", "Brezzi-Douglas-Fortin-Marini", "Crouzeix-Raviart",
                      "Discontinuous Lagrange", "Discontinuous Raviart-Thomas", "HDiv Trace",
                      "Lagrange", "Lobatto", "Nedelec 1st kind H(curl)", "Nedelec 2nd kind H(curl)",
                      "Radau", "Raviart-Thomas", "Real", "Bubble", "Quadrature", "Regge",
                      "Hellan-Herrmann-Johnson", "Q", "DQ", "TensorProductElement", "Gauss-Lobatto-Legendre",
                      "RTCF", "NCF",
                      "RTCE", "NCE")
# The following elements are not supported in FIAT yet, but will be supported here once they are
#     "BDMCF", "BDMCE"

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
    typename = type(element).__module__ + "." + type(element).__name__
    raise ValueError("Element type " + typename + " is not supported.")


@_create_element.register(ufl.FiniteElement)
def _create_finiteelement(element: ufl.FiniteElement):
    """Create FIAT element for UFL base type ufl.FiniteElement."""

    print("Family = ", element.family())
    print("Cell = ", element.cell().cellname())
    print("Degree = ", element.degree())
    libtab_element = libtab.create_element(element.family(),
                                           element.cell().cellname(),
                                           element.degree())

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
        return _flatten_quad(element)
    elif element.cell().cellname() == "hexahedron":
        return _flatten_hex(element)

    if element.family() not in FIAT.supported_elements:
        raise ValueError("Finite element of type \"{}\" is not supported by FIAT.".format(element.family()))

    # Handle Lagrange variants
    if element.family() == "Lagrange" and element.variant() == "spectral":
        assert element.cell().cellname() == "interval"
        element_class = FIAT.GaussLobattoLegendre
    else:
        element_class = FIAT.supported_elements[element.family()]

    assert element.degree() is not None
    print('call class')
    el = element_class(FIAT.ufc_cell(element.cell().cellname()), element.degree())
    print('got el', el)
    el.libtab_element = libtab_element
    return el


@_create_element.register(ufl.MixedElement)
def _create_mixed_finiteelement(element: ufl.MixedElement) -> FIAT.MixedElement:
    elements = []

    def rextract(els):
        for e in els:
            if isinstance(e, ufl.MixedElement) \
                    and not isinstance(e, ufl.VectorElement) \
                    and not isinstance(e, ufl.TensorElement):
                rextract(e.sub_elements())
            else:
                elements.append(e)

    rextract(element.sub_elements())
    print('ELement = ', elements)
    return FIAT.MixedElement(map(_create_element, elements))


@_create_element.register(ufl.VectorElement)
def _create_vector_finiteelement(element: ufl.VectorElement) -> FIAT.MixedElement:
    fiat_element = FIAT.MixedElement(map(_create_element, element.sub_elements()))

    def reorder_for_vector_element(item, block_size):
        """Reorder the elements in item from XXYYZZ ordering to XYZXYZ."""
        space_dim = len(item) // block_size
        return [item[i] for block in range(space_dim)
                for i in range(block, len(item), space_dim)]

    def calculate_entity_dofs_of_vector_element(entity_dofs, block_size):
        """Get the entity DOFs of a VectorElement with XYZXYZ ordering."""
        return {
            dim: {
                entity: [block_size * i + j for i in e_dofs for j in range(block_size)]
                for entity, e_dofs in dofs.items()
            } for dim, dofs in entity_dofs.items()
        }

    # Reorder from XXYYZZ to XYZXYZ
    block_size = fiat_element.num_sub_elements()
    fiat_element.mapping = types.MethodType(
        lambda self: [m for m in self._elements[0].mapping() for e in self._elements],
        fiat_element)
    fiat_element.dual.nodes = reorder_for_vector_element(fiat_element.dual.nodes, block_size)
    fiat_element.dual.entity_ids = calculate_entity_dofs_of_vector_element(
        fiat_element.elements()[0].dual.entity_ids, block_size)
    fiat_element.dual.entity_closure_ids = calculate_entity_dofs_of_vector_element(
        fiat_element.elements()[0].dual.entity_closure_ids, block_size)
    fiat_element.old_tabulate = fiat_element.tabulate

    def tabulate(self, order, points, entity=None):
        block_size = self.num_sub_elements()
        scalar_dofs = len(self.dual.nodes) // block_size
        return {
            i: numpy.array([item[j] for dim in range(scalar_dofs)
                            for j in range(dim, len(item), scalar_dofs)])
            for i, item in self.old_tabulate(order, points, entity=entity).items()
        }

    fiat_element.tabulate = types.MethodType(tabulate, fiat_element)

    print('element = ', fiat_element)

    return fiat_element


@_create_element.register(ufl.HDivElement)
def _create_hdiv_finiteelement(element: ufl.HDivElement) -> FIAT.TensorProductElement:
    tp = _create_element(element._element)
    return FIAT.Hdiv(tp)


@_create_element.register(ufl.HCurlElement)
def _create_hcurl_finiteelement(element: ufl.HCurlElement) -> FIAT.TensorProductElement:
    tp = _create_element(element._element)
    return FIAT.Hcurl(tp)


@_create_element.register(ufl.TensorElement)
def _create_tensor_finiteelement(element: ufl.TensorElement) -> FIAT.MixedElement:
    return _create_vector_finiteelement(element)


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
    element = _create_element(ufl_element)
    _cache[element_signature] = element

    return element


def create_quadrature(shape, degree: int, scheme: str = "default"):
    """Generate quadrature rule.

    Quadrature rule(points, weights) for given shape that will integrate
    an polynomial of order 'degree' exactly.

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

    Map points from the (UFC) reference simplex of dimension d - 1 to a
    given facet on the (UFC) reference simplex of dimension d. This may
    be used to transform points tabulated for example on the 2D
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


def _flatten_quad(element):
    e = _create_element(element.reconstruct(cell=_tpc_quadrilateral))
    flat = FIAT.tensor_product.FlattenedDimensions(e)

    # Overwrite undefined DOF types of Hdiv and Hcurl spaces with correct types
    if element.family() == "RTCF" or element.family() == "BDMCF":
        for dofs in flat.entity_dofs()[1].values():
            for d in dofs:
                flat.dual.nodes[d].functional_type = "PointNormalEval"

    # Overwrite undefined DOF types of Hdiv and Hcurl spaces with correct types
    if element.family() == "RTCE" or element.family() == "BDMCE":
        for dofs in flat.entity_dofs()[1].values():
            for d in dofs:
                flat.dual.nodes[d].functional_type = "PointEdgeTangent"

    return flat


def _flatten_hex(element):
    e = _create_element(element.reconstruct(cell=_tpc_hexahedron))
    flat = FIAT.tensor_product.FlattenedDimensions(e)

    # Overwrite undefined DOF types of Hdiv and Hcurl spaces with correct types
    if element.family() == "NCF":
        for dofs in flat.entity_dofs()[2].values():
            for d in dofs:
                flat.dual.nodes[d].functional_type = "PointNormalEval"

    if element.family() == "NCE":
        for dofs in flat.entity_dofs()[1].values():
            for d in dofs:
                flat.dual.nodes[d].functional_type = "PointEdgeTangent"
        for dofs in flat.entity_dofs()[2].values():
            for d in dofs:
                flat.dual.nodes[d].functional_type = "PointFaceTangent"

    return flat
