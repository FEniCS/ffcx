# Copyright (C) 2021 Matthew W. Scroggs and Chris Richardson
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Finite element interface."""

from __future__ import annotations

import typing
import warnings
from functools import lru_cache


import basix
import basix.ufl
import ufl
import numpy as np
import numpy.typing as npt


def convert_element(element: ufl.finiteelement.FiniteElementBase) -> basix.ufl._ElementBase:
    """Convert and element to a FFCx element."""
    if isinstance(element, basix.ufl._ElementBase):
        return element
    else:
        return _cached_conversion(element)


@lru_cache()
def _cached_conversion(element: ufl.finiteelement.FiniteElementBase) -> basix.ufl._ElementBase:
    """Create an FFCx element from a UFL element.

    Args:
        element: A UFL finite element

    Returns:
        A Basix finite element
    """
    warnings.warn(
        "Use of elements created by UFL is deprecated. You should create elements directly using Basix.",
        DeprecationWarning)

    # Tackle compositional elements, e.g. VectorElement first, then elements
    # implemented by FFCx, then finally elements convertible by Basix.
    if hasattr(ufl, "VectorElement") and isinstance(element, ufl.VectorElement):
        return basix.ufl.blocked_element(
            _cached_conversion(element.sub_elements()[0]), shape=(element.num_sub_elements(), ))
    elif hasattr(ufl, "TensorElement") and isinstance(element, ufl.TensorElement):
        if len(element.symmetry()) == 0:
            return basix.ufl.blocked_element(_cached_conversion(element.sub_elements()[0]), shape=element._value_shape)
        else:
            assert element.symmetry()[(1, 0)] == (0, 1)
            return basix.ufl.blocked_element(_cached_conversion(element.sub_elements()[0]),
                                             element._value_shape, symmetry=True)
    elif hasattr(ufl, "MixedElement") and isinstance(element, ufl.MixedElement):
        return basix.ufl.mixed_element([_cached_conversion(e) for e in element.sub_elements()])
    elif hasattr(ufl, "EnrichedElement") and isinstance(element, ufl.EnrichedElement):
        return basix.ufl.enriched_element([_cached_conversion(e) for e in element._elements])
    elif element.family() == "Quadrature":
        return QuadratureElement(element.cell().cellname(), element.value_shape(), scheme=element.quadrature_scheme(),
                                 degree=element.degree())
    elif element.family() == "Real":
        return RealElement(element)
    else:
        return basix.ufl.convert_ufl_element(element)


def basix_index(indices: typing.Tuple[int]) -> int:
    """Get the Basix index of a derivative."""
    return basix.index(*indices)


def create_quadrature(cellname, degree, rule) -> typing.Tuple[npt.NDArray[np.float64],
                                                              npt.NDArray[np.float64]]:
    """Create a quadrature rule."""
    if cellname == "vertex":
        return (np.ones((1, 0), dtype=np.float64), np.ones(1, dtype=np.float64))
    else:
        return basix.make_quadrature(
            basix.cell.string_to_type(cellname), degree,
            rule=basix.quadrature.string_to_type(rule))


def reference_cell_vertices(cellname: str) -> npt.NDArray[np.float64]:
    """Get the vertices of a reference cell."""
    return basix.geometry(basix.cell.string_to_type(cellname))


def map_facet_points(points: npt.NDArray[np.float64], facet: int, cellname: str) -> npt.NDArray[np.float64]:
    """Map points from a reference facet to a physical facet."""
    geom = basix.geometry(basix.cell.string_to_type(cellname))
    facet_vertices = [geom[i] for i in basix.topology(basix.cell.string_to_type(cellname))[-2][facet]]
    return np.asarray([facet_vertices[0] + sum((i - facet_vertices[0]) * j for i, j in zip(facet_vertices[1:], p))
                       for p in points], dtype=np.float64)


class QuadratureElement(basix.ufl._ElementBase):
    """A quadrature element."""

    _points: npt.NDArray[np.float64]
    _weights: npt.NDArray[np.float64]
    _entity_counts: typing.List[int]
    _cellname: str

    def __init__(self, cellname: str, value_shape: typing.Tuple[int, ...], scheme: typing.Optional[str] = None,
                 degree: typing.Optional[int] = None, points: typing.Optional[npt.NDArray[np.float64]] = None,
                 weights: typing.Optional[npt.NDArray[np.float64]] = None, mapname: str = "identity"):
        """Initialise the element."""
        if scheme is not None:
            assert degree is not None
            assert points is None
            assert weights is None
            repr = f"QuadratureElement({cellname}, {scheme}, {degree})"
            self._points, self._weights = create_quadrature(cellname, degree, scheme)
        else:
            assert degree is None
            assert points is not None
            assert weights is not None
            self._points = points
            self._weights = weights
            repr = f"QuadratureElement({cellname}, {points}, {weights})"
            degree = len(points)

        self._cellname = cellname
        basix_cell = basix.cell.string_to_type(cellname)
        self._entity_counts = [len(i) for i in basix.topology(basix_cell)]

        super().__init__(repr, "quadrature element", cellname, value_shape, degree, mapname=mapname)

    def basix_sobolev_space(self):
        """Return the underlying Sobolev space."""
        return basix.sobolev_spaces.L2

    def __eq__(self, other) -> bool:
        """Check if two elements are equal."""
        return isinstance(other, QuadratureElement) and np.allclose(self._points, other._points) and \
            np.allclose(self._weights, other._weights)

    def __hash__(self) -> int:
        """Return a hash."""
        return super().__hash__()

    def tabulate(self, nderivs: int, points: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
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
        tables = np.asarray([np.eye(points.shape[0], points.shape[0])])
        return tables

    def get_component_element(self, flat_component: int) -> typing.Tuple[basix.ufl._ElementBase, int, int]:
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
    def num_entity_dofs(self) -> typing.List[typing.List[int]]:
        """Number of DOFs associated with each entity."""
        dofs = []
        for d in self._entity_counts[:-1]:
            dofs += [[0] * d]

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
    def reference_geometry(self) -> npt.NDArray[np.float64]:
        """Geometry of the reference element."""
        raise NotImplementedError()

    @property
    def family_name(self) -> str:
        """Family name of the element."""
        return "quadrature"

    @property
    def lagrange_variant(self) -> typing.Union[basix.LagrangeVariant, None]:
        """Basix Lagrange variant used to initialise the element."""
        return None

    @property
    def dpc_variant(self) -> typing.Union[basix.DPCVariant, None]:
        """Basix DPC variant used to initialise the element."""
        return None

    @property
    def element_family(self) -> typing.Union[basix.ElementFamily, None]:
        """Basix element family used to initialise the element."""
        return None

    @property
    def cell_type(self) -> basix.CellType:
        """Basix cell type used to initialise the element."""
        return basix.cell.string_to_type(self._cellname)

    @property
    def discontinuous(self) -> bool:
        """True if the discontinuous version of the element is used."""
        return False

    @property
    def map_type(self) -> basix.MapType:
        """The Basix map type."""
        return basix.MapType.identity


class RealElement(basix.ufl._ElementBase):
    """A real element."""

    _family_name: str
    _cellname: str
    _entity_counts: typing.List[int]

    def __init__(self, element: ufl.finiteelement.FiniteElementBase):
        """Initialise the element."""
        self._cellname = element.cell().cellname()
        self._family_name = element.family()
        tdim = element.cell().topological_dimension()

        self._entity_counts = []
        if tdim >= 1:
            self._entity_counts.append(element.cell().num_vertices())
        if tdim >= 2:
            self._entity_counts.append(element.cell().num_edges())
        if tdim >= 3:
            self._entity_counts.append(element.cell().num_facets())
        self._entity_counts.append(1)

        super().__init__(
            f"RealElement({element})", "real element", element.cell().cellname(), element.value_shape(),
            element.degree())

    def __eq__(self, other) -> bool:
        """Check if two elements are equal."""
        return isinstance(other, RealElement)

    def __hash__(self) -> int:
        """Return a hash."""
        return super().__hash__()

    def tabulate(self, nderivs: int, points: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """Tabulate the basis functions of the element.

        Args:
            nderivs: Number of derivatives to tabulate.
            points: Points to tabulate at

        Returns:
            Tabulated basis functions
        """
        out = np.zeros((nderivs + 1, len(points), 1))
        out[0, :] = 1.
        return out

    def get_component_element(self, flat_component: int) -> typing.Tuple[basix.ufl._ElementBase, int, int]:
        """Get element that represents a component of the element, and the offset and stride of the component.

        Args:
            flat_component: The component

        Returns:
            component element, offset of the component, stride of the component

        """
        assert flat_component < self.value_size
        return self, 0, 1

    @property
    def ufcx_element_type(self) -> str:
        """Element type."""
        return "ufcx_real_element"

    @property
    def dim(self) -> int:
        """Number of DOFs the element has."""
        return 0

    @property
    def num_entity_dofs(self) -> typing.List[typing.List[int]]:
        """Number of DOFs associated with each entity."""
        dofs = []
        for d in self._entity_counts[:-1]:
            dofs += [[0] * d]

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
        return 1

    @property
    def reference_topology(self) -> typing.List[typing.List[typing.List[int]]]:
        """Topology of the reference element."""
        raise NotImplementedError()

    @property
    def reference_geometry(self) -> npt.NDArray[np.float64]:
        """Geometry of the reference element."""
        raise NotImplementedError()

    @property
    def family_name(self) -> str:
        """Family name of the element."""
        return self._family_name

    @property
    def lagrange_variant(self) -> typing.Union[basix.LagrangeVariant, None]:
        """Basix Lagrange variant used to initialise the element."""
        return None

    @property
    def dpc_variant(self) -> typing.Union[basix.DPCVariant, None]:
        """Basix DPC variant used to initialise the element."""
        return None

    @property
    def element_family(self) -> typing.Union[basix.ElementFamily, None]:
        """Basix element family used to initialise the element."""
        return None

    @property
    def cell_type(self) -> basix.CellType:
        """Basix cell type used to initialise the element."""
        return basix.cell.string_to_type(self._cellname)

    @property
    def discontinuous(self) -> bool:
        """True if the discontinuous version of the element is used."""
        return False

    def basix_sobolev_space(self):
        """Return the underlying Sobolev space."""
        return basix.sobolev_spaces.Hinf

    @property
    def map_type(self) -> basix.MapType:
        """The Basix map type."""
        return basix.MapType.identity
