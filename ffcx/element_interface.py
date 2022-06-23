
# Copyright (C) 2021 Matthew W. Scroggs and Chris Richardson
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Finite element interface."""

from __future__ import annotations

import typing

import warnings

import basix
import numpy
import ufl
import basix.ufl_wrapper


def create_element(element: ufl.finiteelement.FiniteElementBase) -> basix.ufl_wrapper._BasixElementBase:
    """Create an FFCx element from a UFL element.

    Args:
        element: A UFL finite element

    Returns:
        A FFCx finite element
    """
    # TODO: EnrichedElement

    if element.family() == "Quadrature":
        return QuadratureElement(element)

    return basix.ufl_wrapper.convert_ufl_element(element)


def basix_index(indices: typing.Tuple[int]) -> int:
    """Get the Basix index of a derivative."""
    return basix.index(*indices)


def create_quadrature(cellname, degree, rule) -> typing.Tuple[numpy.typing.NDArray[numpy.float64],
                                                              numpy.typing.NDArray[numpy.float64]]:
    """Create a quadrature rule."""
    if cellname == "vertex":
        return (numpy.ones((1, 0), dtype=numpy.float64), numpy.ones(1, dtype=numpy.float64))

    quadrature = basix.make_quadrature(
        basix.quadrature.string_to_type(rule), basix.cell.string_to_type(cellname), degree)

    # The quadrature degree from UFL can be very high for some
    # integrals.  Print warning if number of quadrature points
    # exceeds 100.
    num_points = quadrature[1].size
    if num_points >= 100:
        warnings.warn(
            f"Number of integration points per cell is: {num_points}. Consider using 'quadrature_degree' "
            "to reduce number.")
    return quadrature


def reference_cell_vertices(cellname: str) -> numpy.typing.NDArray[numpy.float64]:
    """Get the vertices of a reference cell."""
    return basix.geometry(basix.cell.string_to_type(cellname))


def map_facet_points(points: numpy.typing.NDArray[numpy.float64], facet: int,
                     cellname: str) -> numpy.typing.NDArray[numpy.float64]:
    """Map points from a reference facet to a physical facet."""
    geom = basix.geometry(basix.cell.string_to_type(cellname))
    facet_vertices = [geom[i] for i in basix.topology(basix.cell.string_to_type(cellname))[-2][facet]]
    return numpy.asarray([facet_vertices[0] + sum((i - facet_vertices[0]) * j for i, j in zip(facet_vertices[1:], p))
                          for p in points], dtype=numpy.float64)


class QuadratureElement(basix.ufl_wrapper._BasixElementBase):
    """A quadrature element."""

    _points: basix.ufl_wrapper._nda_f64
    _element: basix.ufl_wrapper._BasixElementBase

    def __init__(self, element: basix.ufl_wrapper._BasixElementBase):
        """Initialise the element."""
        self._points, _ = basix.make_quadrature(element.cell_type, element.degree)
        self._element = element

    def tabulate(
        self, nderivs: int, points: basix.ufl_wrapper._nda_f64
    ) -> basix.ufl_wrapper._nda_f64:
        """Tabulate the basis functions of the element.

        Args:
            nderivs: Number of derivatives to tabulate.
            points: Points to tabulate at

        Returns:
            Tabulated basis functions
        """
        if nderivs > 0:
            raise ValueError("Cannot take derivatives of Quadrature element.")

        if points.shape != self._points.shape:
            raise ValueError("Mismatch of tabulation points and element points.")
        tables = numpy.asarray([numpy.eye(points.shape[0], points.shape[0])])
        return tables

    def get_component_element(self, flat_component: int) -> typing.Tuple[basix.ufl_wrapper._BasixElementBase, int, int]:
        """Get element that represents a component of the element, and the offset and stride of the component.

        Args:
            flat_component: The component

        Returns:
            component element, offset of the component, stride of the component
        """
        return self, 0, 1

    @property
    def ufcx_element_type(self) -> str:
        """Element type."""
        return "ufcx_quadrature_element"

    @property
    def dim(self) -> int:
        """Number of DOFs the element has."""
        return self._points.shape[0]

    @property
    def value_size(self) -> int:
        """Value size of the element.

        Equal to ``numpy.prod(value_shape)``.

        """
        return 1

    @property
    def value_shape(self) -> typing.Tuple[int, ...]:
        """Value shape of the element basis function.

        Note:
            For scalar elements, ``(1,)`` is returned. This is different
            from Basix where the value shape for scalar elements is
            ``(,)``.
        """
        return (1,)

    @property
    def num_entity_dofs(self) -> typing.List[typing.List[int]]:
        """Number of DOFs associated with each entity."""
        dofs = []
        tdim = self._element.cell().topological_dimension()

        if tdim >= 1:
            dofs += [[0] * self._element.cell().num_vertices()]

        if tdim >= 2:
            dofs += [[0] * self._element.cell().num_edges()]

        if tdim >= 3:
            dofs += [[0] * self._element.cell().num_facets()]

        dofs += [[self.dim]]
        return dofs

    @property
    def entity_dofs(self) -> typing.List[typing.List[typing.List[int]]]:
        """DOF numbers associated with each entity."""
        start_dof = 0
        entity_dofs = []
        for i in self.num_entity_dofs:
            dofs_list = []
            for j in i:
                dofs_list.append([start_dof + k for k in range(j)])
                start_dof += j
            entity_dofs.append(dofs_list)
        return entity_dofs

    @property
    def num_entity_closure_dofs(self) -> typing.List[typing.List[int]]:
        """Number of DOFs associated with the closure of each entity."""
        return self.num_entity_dofs

    @property
    def entity_closure_dofs(self) -> typing.List[typing.List[typing.List[int]]]:
        """DOF numbers associated with the closure of each entity."""
        return self.entity_dofs

    @property
    def num_global_support_dofs(self) -> int:
        """Get the number of global support DOFs."""
        return 0

    @property
    def reference_topology(self) -> typing.List[typing.List[typing.List[int]]]:
        """Topology of the reference element."""
        raise NotImplementedError()

    @property
    def reference_geometry(self) -> basix.ufl_wrapper._nda_f64:
        """Geometry of the reference element."""
        raise NotImplementedError()

    @property
    def family_name(self) -> str:
        """Family name of the element."""
        return self._element.family()

    @property
    def lagrange_variant(self) -> basix.LagrangeVariant:
        """Basix Lagrange variant used to initialise the element."""
        return None

    @property
    def dpc_variant(self) -> basix.DPCVariant:
        """Basix DPC variant used to initialise the element."""
        return None

    @property
    def element_family(self) -> basix.ElementFamily:
        """Basix element family used to initialise the element."""
        return None

    @property
    def cell_type(self) -> basix.CellType:
        """Basix cell type used to initialise the element."""
        return None

    @property
    def discontinuous(self) -> bool:
        """True if the discontinuous version of the element is used."""
        return False

    @property
    def interpolation_nderivs(self) -> int:
        """The number of derivatives needed when interpolating."""
        return 0
