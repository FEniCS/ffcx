# Copyright (C) 2021 Matthew W. Scroggs and Chris Richardson
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Finite element interface."""

from __future__ import annotations

import typing

if typing.TYPE_CHECKING:
    import ufl.finiteelement.FiniteElementBase

import warnings
from abc import ABC, abstractmethod

import basix
import numpy
import ufl
import basix.ufl_wrapper
import functools


@functools.lru_cache()
def create_basix_element(family_type, cell_type, degree, variant_info, discontinuous):
    """Create a basix element."""
    return basix.create_element(family_type, cell_type, degree, *variant_info, discontinuous)


def create_element(element: ufl.finiteelement.FiniteElementBase) -> BaseElement:
    """Create an FFCx element from a UFL element.

    Args:
        element: A UFL finite element

    Returns:
        A FFCx finite element
    """
    # TODO: EnrichedElement

    if isinstance(element, basix.ufl_wrapper.BasixElement):
        return BasixElement(element.basix_element)

    if isinstance(element, ufl.VectorElement):
        return BlockedElement(create_element(element.sub_elements()[0]),
                              element.num_sub_elements())
    if isinstance(element, ufl.TensorElement):
        return BlockedElement(create_element(element.sub_elements()[0]),
                              element.num_sub_elements(), None)  # TODO: block shape

    if isinstance(element, ufl.MixedElement):
        return MixedElement([create_element(e) for e in element.sub_elements()])

    if element.family() == "Quadrature":
        return QuadratureElement(element)

    family_name = element.family()
    discontinuous = False
    if family_name.startswith("Discontinuous "):
        family_name = family_name[14:]
        discontinuous = True
    if family_name == "DP":
        family_name = "P"
        discontinuous = True
    if family_name == "DQ":
        family_name = "Q"
        discontinuous = True
    if family_name == "DPC":
        discontinuous = True

    family_type = basix.finite_element.string_to_family(family_name, element.cell().cellname())
    cell_type = basix.cell.string_to_type(element.cell().cellname())

    variant_info = []
    if family_type == basix.ElementFamily.P and element.variant() == "equispaced":
        # This is used for elements defining cells
        variant_info = [basix.LagrangeVariant.equispaced]
    else:
        if element.variant() is not None:
            raise ValueError("UFL variants are not supported by FFCx. Please wrap a Basix element directly.")

        EF = basix.ElementFamily
        if family_type == EF.P:
            variant_info = [basix.LagrangeVariant.gll_warped]
        elif family_type in [EF.RT, EF.N1E]:
            variant_info = [basix.LagrangeVariant.legendre]
        elif family_type in [EF.serendipity, EF.BDM, EF.N2E]:
            variant_info = [basix.LagrangeVariant.legendre, basix.DPCVariant.legendre]
        elif family_type == EF.DPC:
            variant_info = [basix.DPCVariant.diagonal_gll]

    return BasixElement(create_basix_element(
        family_type, cell_type, element.degree(), tuple(variant_info), discontinuous))


def basix_index(*args):
    """Get the Basix index of a derivative."""
    return basix.index(*args)


def create_quadrature(cellname, degree, rule):
    """Create a quadrature rule."""
    if cellname == "vertex":
        return [[]], [1]

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


def reference_cell_vertices(cellname):
    """Get the vertices of a reference cell."""
    return basix.geometry(basix.cell.string_to_type(cellname))


def map_facet_points(points, facet, cellname):
    """Map points from a reference facet to a physical facet."""
    geom = basix.geometry(basix.cell.string_to_type(cellname))
    facet_vertices = [geom[i] for i in basix.topology(basix.cell.string_to_type(cellname))[-2][facet]]
    return [facet_vertices[0] + sum((i - facet_vertices[0]) * j for i, j in zip(facet_vertices[1:], p))
            for p in points]


class BaseElement(ABC):
    """An abstract element class."""

    @abstractmethod
    def tabulate(self, nderivs: int, points: numpy.ndarray):
        """Tabulate the basis functions of the element.

        Args:
            nderivs: Number of derivatives to tabulate.
            points: Points to tabulate at

        Returns:
            Tabulated basis functions
        """
        pass

    @abstractmethod
    def get_component_element(self, flat_component: int) -> tuple[BaseElement, int, int]:
        """Get element that represents a component of the element, and the offset and stride of the component.

        For example, for a MixedElement, this will return the
        sub-element that represents the given component, the offset of
        that sub-element, and a stride of 1. For a BlockedElement, this
        will return the sub-element, an offset equal to the component
        number, and a stride equal to the block size. For vector-valued
        element (eg H(curl) and H(div) elements), this returns a
        ComponentElement (and as offset of 0 and a stride of 1). When
        tabulate is called on the ComponentElement, only the part of the
        table for the given component is returned.

        Args:
            flat_component: The component

        Returns:
            component element, offset of the component, stride of the component
        """
        pass

    @property
    def element_type(self):
        """Element type."""
        raise NotImplementedError

    @property
    @abstractmethod
    def dim(self) -> int:
        """Number of DOFs the element has."""
        pass

    @property
    @abstractmethod
    def value_size(self) -> int:
        """Value size of the element.

        Equal to ``numpy.prod(value_shape)``.

        """
        pass

    @property
    @abstractmethod
    def value_shape(self) -> typing.Tuple[int, ...]:
        """Value shape of the element basis function.

        Note:
            For scalar elements, ``(1,)`` is returned. This is different
            from Basix where the value shape for scalar elements is
            ``(,)``.

        """
        pass

    @property
    @abstractmethod
    def num_entity_dofs(self):
        """Number of DOFs associated with each entity."""
        pass

    @property
    @abstractmethod
    def entity_dofs(self):
        """DOF numbers associated with each entity."""
        pass

    @property
    @abstractmethod
    def num_entity_closure_dofs(self):
        """Number of DOFs associated with the closure of each entity."""
        pass

    @property
    @abstractmethod
    def entity_closure_dofs(self):
        """DOF numbers associated with the closure of each entity."""
        pass

    @property
    @abstractmethod
    def num_global_support_dofs(self):
        pass

    @property
    @abstractmethod
    def reference_topology(self):
        """Topology of the reference element."""
        pass

    @property
    @abstractmethod
    def reference_geometry(self):
        """Geometry of the reference element."""
        pass

    @property
    @abstractmethod
    def family_name(self) -> str:
        """Family name of the element."""
        pass

    @property
    @abstractmethod
    def element_family(self):
        """Basix element family used to initialise the element."""
        pass

    @property
    @abstractmethod
    def lagrange_variant(self):
        """Basix Lagrange variant used to initialise the element."""
        pass

    @property
    @abstractmethod
    def dpc_variant(self):
        """Basix DPC variant used to initialise the element."""
        pass

    @property
    @abstractmethod
    def cell_type(self):
        """Basix cell type used to initialise the element."""
        pass

    @property
    @abstractmethod
    def discontinuous(self) -> bool:
        """True if the discontinuous version of the element is used."""
        pass

    @property
    def is_custom_element(self) -> bool:
        """True if the element is a custom Basix element."""
        return False


class BasixElement(BaseElement):
    """An element defined by Basix."""

    def __init__(self, element):
        self.element = element

    def tabulate(self, nderivs, points):
        tab = self.element.tabulate(nderivs, points)
        return tab.transpose((0, 1, 3, 2)).reshape((tab.shape[0], tab.shape[1], -1))

    def get_component_element(self, flat_component):
        assert flat_component < self.value_size
        return ComponentElement(self, flat_component), 0, 1

    @property
    def element_type(self) -> str:
        """Element type."""
        if self.is_custom_element:
            return "ufcx_basix_custom_element"
        else:
            return "ufcx_basix_element"

    @property
    def dim(self):
        return self.element.dim

    @property
    def value_size(self):
        return self.element.value_size

    @property
    def value_shape(self):
        """Get the value shape of the element."""
        if len(self.element.value_shape) == 0:
            return (1,)
        else:
            return self.element.value_shape

    @property
    def num_entity_dofs(self):
        return self.element.num_entity_dofs

    @property
    def entity_dofs(self):
        return self.element.entity_dofs

    @property
    def num_entity_closure_dofs(self):
        return self.element.num_entity_closure_dofs

    @property
    def entity_closure_dofs(self):
        return self.element.entity_closure_dofs

    @property
    def num_global_support_dofs(self):
        # TODO
        return 0

    @property
    def family_name(self):
        return self.element.family.name

    @property
    def reference_topology(self):
        return basix.topology(self.element.cell_type)

    @property
    def reference_geometry(self):
        return basix.geometry(self.element.cell_type)

    @property
    def element_family(self):
        return self.element.family

    @property
    def lagrange_variant(self):
        return self.element.lagrange_variant

    @property
    def dpc_variant(self):
        return self.element.dpc_variant

    @property
    def cell_type(self):
        return self.element.cell_type

    @property
    def discontinuous(self):
        return self.element.discontinuous

    @property
    def is_custom_element(self) -> bool:
        """True if the element is a custom Basix element."""
        return self.element.family == basix.ElementFamily.custom


class ComponentElement(BaseElement):
    """An element representing one component of a BasixElement."""

    def __init__(self, element, component):
        self.element = element
        self.component = component

    def tabulate(self, nderivs, points):
        tables = self.element.tabulate(nderivs, points)
        output = []
        for tbl in tables:
            shape = (tbl.shape[0],) + tuple(self.element.value_shape) + (-1,)
            tbl = tbl.reshape(shape)
            if len(self.element.value_shape) == 1:
                output.append(tbl[:, self.component, :])
            elif len(self.element.value_shape) == 2:
                # TODO: Something different may need doing here if
                # tensor is symmetric
                vs0 = self.element.value_shape[0]
                output.append(tbl[:, self.component // vs0, self.component % vs0, :])
            else:
                raise NotImplementedError
        return output

    def get_component_element(self, flat_component):
        if flat_component == 0:
            return self, 0, 1
        raise NotImplementedError

    @property
    def dim(self):
        raise NotImplementedError

    @property
    def value_size(self):
        raise NotImplementedError

    @property
    def value_shape(self):
        raise NotImplementedError

    @property
    def num_entity_dofs(self):
        raise NotImplementedError

    @property
    def entity_dofs(self):
        raise NotImplementedError

    @property
    def num_entity_closure_dofs(self):
        raise NotImplementedError

    @property
    def entity_closure_dofs(self):
        raise NotImplementedError

    @property
    def num_global_support_dofs(self):
        raise NotImplementedError

    @property
    def family_name(self) -> str:
        raise NotImplementedError

    @property
    def reference_topology(self):
        raise NotImplementedError

    @property
    def reference_geometry(self):
        raise NotImplementedError

    @property
    def element_family(self):
        return self.element.element_family

    @property
    def lagrange_variant(self):
        return self.element.lagrange_variant

    @property
    def dpc_variant(self):
        return self.element.dpc_variant

    @property
    def cell_type(self):
        return self.element.cell_type

    @property
    def discontinuous(self):
        return self.element.discontinuous


class MixedElement(BaseElement):
    """A mixed element that combines two or more elements."""

    def __init__(self, sub_elements):
        assert len(sub_elements) > 0
        self.sub_elements = sub_elements

    def tabulate(self, nderivs, points):
        tables = []
        results = [e.tabulate(nderivs, points) for e in self.sub_elements]
        for deriv_tables in zip(*results):
            new_table = numpy.zeros((len(points), self.value_size * self.dim))
            start = 0
            for e, t in zip(self.sub_elements, deriv_tables):
                for i in range(0, e.dim, e.value_size):
                    new_table[:, start: start + e.value_size] = t[:, i: i + e.value_size]
                    start += self.value_size
            tables.append(new_table)
        return tables

    def get_component_element(self, flat_component):
        sub_dims = [0] + [e.dim for e in self.sub_elements]
        sub_cmps = [0] + [e.value_size for e in self.sub_elements]

        irange = numpy.cumsum(sub_dims)
        crange = numpy.cumsum(sub_cmps)

        # Find index of sub element which corresponds to the current
        # flat component
        component_element_index = numpy.where(
            crange <= flat_component)[0].shape[0] - 1

        sub_e = self.sub_elements[component_element_index]

        e, offset, stride = sub_e.get_component_element(flat_component - crange[component_element_index])
        # TODO: is this offset correct?
        return e, irange[component_element_index] + offset, stride

    @property
    def element_type(self) -> str:
        """Get the element type."""
        return "ufcx_mixed_element"

    @property
    def dim(self):
        return sum(e.dim for e in self.sub_elements)

    @property
    def value_size(self):
        return sum(e.value_size for e in self.sub_elements)

    @property
    def value_shape(self):
        return (sum(e.value_size for e in self.sub_elements), )

    @property
    def num_entity_dofs(self):
        data = [e.num_entity_dofs for e in self.sub_elements]
        return [[sum(d[tdim][entity_n] for d in data) for entity_n, _ in enumerate(entities)]
                for tdim, entities in enumerate(data[0])]

    @property
    def entity_dofs(self):
        dofs = [[[] for i in entities] for entities in self.sub_elements[0].entity_dofs]
        start_dof = 0
        for e in self.sub_elements:
            for tdim, entities in enumerate(e.entity_dofs):
                for entity_n, entity_dofs in enumerate(entities):
                    dofs[tdim][entity_n] += [start_dof + i for i in entity_dofs]
            start_dof += e.dim
        return dofs

    @property
    def num_entity_closure_dofs(self):
        data = [e.num_entity_closure_dofs for e in self.sub_elements]
        return [[sum(d[tdim][entity_n] for d in data) for entity_n, _ in enumerate(entities)]
                for tdim, entities in enumerate(data[0])]

    @property
    def entity_closure_dofs(self):
        dofs = [[[] for i in entities] for entities in self.sub_elements[0].entity_closure_dofs]
        start_dof = 0
        for e in self.sub_elements:
            for tdim, entities in enumerate(e.entity_closure_dofs):
                for entity_n, entity_dofs in enumerate(entities):
                    dofs[tdim][entity_n] += [start_dof + i for i in entity_dofs]
            start_dof += e.dim
        return dofs

    @property
    def num_global_support_dofs(self):
        return sum(e.num_global_support_dofs for e in self.sub_elements)

    @property
    def family_name(self):
        return "mixed element"

    @property
    def reference_topology(self):
        return self.sub_elements[0].reference_topology

    @property
    def reference_geometry(self):
        return self.sub_elements[0].reference_geometry

    @property
    def lagrange_variant(self):
        return None

    @property
    def dpc_variant(self):
        return None

    @property
    def element_family(self):
        return None

    @property
    def cell_type(self):
        return None

    @property
    def discontinuous(self):
        return False


class BlockedElement(BaseElement):
    """An element with a block size that contains multiple copies of a sub element."""

    def __init__(self, sub_element, block_size, block_shape=None):
        assert block_size > 0
        if sub_element.value_size != 1:
            raise ValueError("Blocked elements (VectorElement and TensorElement) of "
                             "non-scalar elements are not supported. Try using MixedElement "
                             "instead.")

        self.sub_element = sub_element
        self.block_size = block_size
        if block_shape is None:
            self.block_shape = (block_size, )
        else:
            self.block_shape = block_shape

    def tabulate(self, nderivs, points):
        assert len(self.block_shape) == 1  # TODO: block shape
        assert self.value_size == self.block_size  # TODO: remove this assumption
        output = []
        for table in self.sub_element.tabulate(nderivs, points):
            new_table = numpy.zeros((table.shape[0], table.shape[1] * self.block_size**2))
            for block in range(self.block_size):
                col = block * (self.block_size + 1)
                new_table[:, col: col + table.shape[1] * self.block_size**2: self.block_size**2] = table
            output.append(new_table)
        return output

    def get_component_element(self, flat_component):
        return self.sub_element, flat_component, self.block_size

    @property
    def element_type(self):
        """Element type."""
        return self.sub_element.element_type

    @property
    def dim(self):
        return self.sub_element.dim * self.block_size

    @property
    def value_size(self):
        return self.block_size * self.sub_element.value_size

    @property
    def value_shape(self):
        return (self.value_size, )

    @property
    def num_entity_dofs(self):
        return [[j * self.block_size for j in i] for i in self.sub_element.num_entity_dofs]

    @property
    def entity_dofs(self):
        # TODO: should this return this, or should it take blocks into
        # account?
        return [[[k * self.block_size + b for k in j for b in range(self.block_size)]
                 for j in i] for i in self.sub_element.entity_dofs]

    @property
    def num_entity_closure_dofs(self):
        return [[j * self.block_size for j in i] for i in self.sub_element.num_entity_closure_dofs]

    @property
    def entity_closure_dofs(self):
        # TODO: should this return this, or should it take blocks into
        # account?
        return [[[k * self.block_size + b for k in j for b in range(self.block_size)]
                 for j in i] for i in self.sub_element.entity_closure_dofs]

    @property
    def num_global_support_dofs(self):
        return self.sub_element.num_global_support_dofs * self.block_size

    @property
    def family_name(self):
        return self.sub_element.family_name

    @property
    def reference_topology(self):
        return self.sub_element.reference_topology

    @property
    def reference_geometry(self):
        return self.sub_element.reference_geometry

    @property
    def lagrange_variant(self):
        return self.sub_element.lagrange_variant

    @property
    def dpc_variant(self):
        return self.sub_element.dpc_variant

    @property
    def element_family(self):
        return self.sub_element.element_family

    @property
    def cell_type(self):
        return self.sub_element.cell_type

    @property
    def discontinuous(self):
        return self.sub_element.discontinuous


class QuadratureElement(BaseElement):
    """A quadrature element."""

    def __init__(self, ufl_element):
        self._points, _ = create_quadrature(ufl_element.cell().cellname(),
                                            ufl_element.degree(), ufl_element.quadrature_scheme())
        self._ufl_element = ufl_element

    def tabulate(self, nderivs, points):
        if nderivs > 0:
            raise ValueError("Cannot take derivatives of Quadrature element.")

        if points.shape != self._points.shape:
            raise ValueError("Mismatch of tabulation points and element points.")
        tables = [numpy.eye(points.shape[0], points.shape[0])]
        return tables

    def get_component_element(self, flat_component):
        return self, 0, 1

    @property
    def element_type(self) -> str:
        """Element type."""
        return "ufcx_quadrature_element"

    @property
    def dim(self):
        return self._points.shape[0]

    @property
    def value_size(self):
        return 1

    @property
    def value_shape(self):
        return (1,)

    @property
    def num_entity_dofs(self):
        dofs = []
        tdim = self._ufl_element.cell().topological_dimension()

        if tdim >= 1:
            dofs += [[0] * self._ufl_element.cell().num_vertices()]

        if tdim >= 2:
            dofs += [[0] * self._ufl_element.cell().num_edges()]

        if tdim >= 3:
            dofs += [[0] * self._ufl_element.cell().num_facets()]

        dofs += [[self.dim]]
        return dofs

    @property
    def entity_dofs(self):
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
    def num_entity_closure_dofs(self):
        return self.num_entity_dofs

    @property
    def entity_closure_dofs(self):
        return self.entity_dofs

    @property
    def num_global_support_dofs(self):
        return 0

    @property
    def reference_topology(self):
        raise NotImplementedError

    @property
    def reference_geometry(self):
        raise NotImplementedError

    @property
    def family_name(self):
        return self._ufl_element.family()

    @property
    def lagrange_variant(self):
        return None

    @property
    def dpc_variant(self):
        return None

    @property
    def element_family(self):
        return None

    @property
    def cell_type(self) -> None:
        return None

    @property
    def discontinuous(self):
        return False
