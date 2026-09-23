# Copyright (C) 2021 Matthew Scroggs
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Geometry."""

from collections.abc import Iterable

import basix
import basix.ufl
import numpy as np
import ufl

import ffcx.codegeneration.lnodes as L
from ffcx.definitions import entity_types
from ffcx.ir.integral import IntermediateIntegrandIR

# Per entity_type: the table a mixed-dimensional submesh's coordinate
# dofs are gathered through, and the codimension it must have. Keyed by
# "vertex" rather than "peak", which is not an entity_type in FFCx.
# Codimension 2 is the cap: `ffcx/ir/elementtables.py` rejects anything
# higher, so a 3D parent never reaches codegen.
CLOSURE_DOFS_TABLES: dict[str, tuple[str, int]] = {
    "facet": ("facet_closure_dofs", 1),
    "ridge": ("ridge_closure_dofs", 2),
    "vertex": ("peak_closure_dofs", 2),
}

# Orientations ("quadrature_permutation" values) a sub-entity can have:
# rotations x reflections. A property of the entity, not of its parent cell.
_ENTITY_NPERM = {
    basix.CellType.point: 1,
    basix.CellType.interval: 2,
    basix.CellType.triangle: 6,
    basix.CellType.quadrilateral: 8,
}


_CLOSURE_TABLE_NAMES = {name for name, _ in CLOSURE_DOFS_TABLES.values()}


def write_table(tablename, cellname, coordinate_element=None):
    """Write a table.

    `coordinate_element` is required for the closure-dofs tables and
    ignored by every other table kind.
    """
    if tablename in _CLOSURE_TABLE_NAMES and coordinate_element is None:
        raise ValueError(f"Geometry table {tablename!r} requires a coordinate element.")
    if tablename == "facet_edge_vertices":
        return facet_edge_vertices(tablename, cellname)
    if tablename == "cell_facet_jacobian":
        return cell_facet_jacobian(tablename, cellname)
    if tablename == "cell_ridge_jacobian":
        return cell_ridge_jacobian(tablename, cellname)
    if tablename == "facet_closure_dofs":
        return facet_closure_dofs(tablename, cellname, coordinate_element)
    if tablename == "ridge_closure_dofs":
        return ridge_closure_dofs(tablename, cellname, coordinate_element)
    if tablename == "peak_closure_dofs":
        return peak_closure_dofs(tablename, cellname, coordinate_element)
    if tablename == "reference_cell_volume":
        return reference_cell_volume(tablename, cellname)
    if tablename == "reference_facet_volume":
        return reference_facet_volume(tablename, cellname)
    if tablename == "reference_cell_edge_vectors":
        return reference_cell_edge_vectors(tablename, cellname)
    if tablename == "reference_facet_edge_vectors":
        return reference_facet_edge_vectors(tablename, cellname)
    if tablename == "reference_normals":
        return reference_normals(tablename, cellname)
    if tablename == "facet_orientation":
        return facet_orientation(tablename, cellname)
    raise ValueError(f"Unknown geometry table name: {tablename}")


def facet_edge_vertices(tablename, cellname):
    """Write facet edge vertices."""
    celltype = getattr(basix.CellType, cellname)
    topology = basix.topology(celltype)
    triangle_edges = basix.topology(basix.CellType.triangle)[1]
    quadrilateral_edges = basix.topology(basix.CellType.quadrilateral)[1]

    if len(topology) != 4:
        raise ValueError("Can only get facet edges for 3D cells.")

    edge_vertices = []
    for facet in topology[-2]:
        if len(facet) == 3:
            edge_vertices += [[[facet[i] for i in edge] for edge in triangle_edges]]
        elif len(facet) == 4:
            edge_vertices += [[[facet[i] for i in edge] for edge in quadrilateral_edges]]
        else:
            raise ValueError("Only triangular and quadrilateral faces supported.")

    out = np.array(edge_vertices, dtype=int)
    symbol = L.Symbol(f"{cellname}_{tablename}", dtype=L.DataType.INT)
    return L.ArrayDecl(symbol, values=out, const=True)


def cell_facet_jacobian(tablename, cellname):
    """Write a reference facet jacobian."""
    celltype = getattr(basix.CellType, cellname)
    out = basix.cell.facet_jacobians(celltype)
    symbol = L.Symbol(f"{cellname}_{tablename}", dtype=L.DataType.REAL)
    return L.ArrayDecl(symbol, values=out, const=True)


def _scalar_basix_element(coordinate_element):
    """Return the scalar (unblocked) basix FiniteElement of a coordinate element."""
    sub_elements = coordinate_element.sub_elements
    scalar_element = sub_elements[0] if sub_elements else coordinate_element
    return scalar_element.basix_element


def _closure_dofs_table(tablename, cellname, coordinate_element, entity_dim):
    """Write a per-orientation table of a cell's own sub-entity closure dofs.

    One row per possible sub-entity orientation ("quadrature_permutation"
    value), one column per sub-entity of dimension `entity_dim`, listing
    the parent cell's own scalar coordinate dof indices in that
    sub-entity's closure.

    A co-dimensional entity's coordinate dofs are always a subset of its
    parent cell's. The permutation resolves the mismatch between the
    submesh's canonically oriented dofmap and the parent cell's local
    sub-entity closure.
    """
    celltype = getattr(basix.CellType, cellname)
    entity_celltypes = set(basix.cell.subentity_types(celltype)[entity_dim])
    if len(entity_celltypes) != 1:
        raise NotImplementedError(
            f"Mixed-dimensional integrals with a submesh domain are not supported for "
            f"cell type {cellname!r} (entity dimension {entity_dim}): its sub-entities are "
            f"not all the same cell type."
        )
    (entity_celltype,) = entity_celltypes
    nperm = _ENTITY_NPERM[entity_celltype]

    be = _scalar_basix_element(coordinate_element)
    if not be.dof_transformations_are_permutations:
        raise NotImplementedError(
            "Mixed-dimensional coordinate-dofs gathering requires a coordinate element whose "
            "dof transformations are permutations."
        )
    topology = basix.topology(celltype)
    num_entities = len(topology[entity_dim])
    closure_dofs = be.entity_closure_dofs[entity_dim]
    ndofs = len(closure_dofs[0])

    out = np.zeros((nperm, num_entities, ndofs), dtype=int)
    for e in range(num_entities):
        base = np.array(closure_dofs[e], dtype=np.int32)
        for p in range(nperm):
            permuted = base.copy()
            be.permute_subentity_closure_inv(permuted, p, entity_celltype, None)
            out[p, e, :] = permuted

    symbol = L.Symbol(f"{cellname}_{tablename}", dtype=L.DataType.INT)
    return L.ArrayDecl(symbol, values=out, const=True)


def _cell_tdim(cellname):
    return len(basix.topology(getattr(basix.CellType, cellname))) - 1


def facet_closure_dofs(tablename, cellname, coordinate_element):
    """Write a facet-closure-dofs table (see `_closure_dofs_table`)."""
    return _closure_dofs_table(tablename, cellname, coordinate_element, _cell_tdim(cellname) - 1)


def ridge_closure_dofs(tablename, cellname, coordinate_element):
    """Write a ridge-closure-dofs table (see `_closure_dofs_table`)."""
    return _closure_dofs_table(tablename, cellname, coordinate_element, _cell_tdim(cellname) - 2)


def peak_closure_dofs(tablename, cellname, coordinate_element):
    """Write a peak-closure-dofs table (see `_closure_dofs_table`).

    A vertex has one orientation, so this table has a single row and
    `quadrature_permutation` is never consulted to index it.
    """
    return _closure_dofs_table(tablename, cellname, coordinate_element, 0)


def cell_ridge_jacobian(tablename, cellname):
    """Write a reference ridge jacobian."""
    celltype = getattr(basix.CellType, cellname)
    out = basix.cell.edge_jacobians(celltype)
    symbol = L.Symbol(f"{cellname}_{tablename}", dtype=L.DataType.REAL)
    return L.ArrayDecl(symbol, values=out, const=True)


def reference_cell_volume(tablename, cellname):
    """Write a reference cell volume."""
    celltype = getattr(basix.CellType, cellname)
    out = basix.cell.volume(celltype)
    symbol = L.Symbol(f"{cellname}_{tablename}", dtype=L.DataType.REAL)
    return L.VariableDecl(symbol, out)


def reference_facet_volume(tablename, cellname):
    """Write a reference facet volume."""
    celltype = getattr(basix.CellType, cellname)
    volumes = basix.cell.facet_reference_volumes(celltype)
    for i in volumes[1:]:
        if not np.isclose(i, volumes[0]):
            raise ValueError("Reference facet volume not supported for this cell type.")
    symbol = L.Symbol(f"{cellname}_{tablename}", dtype=L.DataType.REAL)
    return L.VariableDecl(symbol, volumes[0])


def reference_cell_edge_vectors(tablename, cellname):
    """Write reference edge vectors."""
    celltype = getattr(basix.CellType, cellname)
    topology = basix.topology(celltype)
    geometry = basix.geometry(celltype)
    edge_vectors = [geometry[j] - geometry[i] for i, j in topology[1]]
    out = np.array(edge_vectors)
    symbol = L.Symbol(f"{cellname}_{tablename}", dtype=L.DataType.REAL)
    return L.ArrayDecl(symbol, values=out, const=True)


def reference_facet_edge_vectors(tablename, cellname):
    """Write facet reference edge vectors."""
    celltype = getattr(basix.CellType, cellname)
    topology = basix.topology(celltype)
    geometry = basix.geometry(celltype)
    triangle_edges = basix.topology(basix.CellType.triangle)[1]
    quadrilateral_edges = basix.topology(basix.CellType.quadrilateral)[1]

    if len(topology) != 4:
        raise ValueError("Can only get facet edges for 3D cells.")

    edge_vectors = []
    for facet in topology[-2]:
        if len(facet) == 3:
            edge_vectors += [geometry[facet[j]] - geometry[facet[i]] for i, j in triangle_edges]
        elif len(facet) == 4:
            edge_vectors += [
                geometry[facet[j]] - geometry[facet[i]] for i, j in quadrilateral_edges
            ]
        else:
            raise ValueError("Only triangular and quadrilateral faces supported.")

    out = np.array(edge_vectors)
    symbol = L.Symbol(f"{cellname}_{tablename}", dtype=L.DataType.REAL)
    return L.ArrayDecl(symbol, values=out, const=True)


def reference_normals(tablename, cellname):
    """Write reference facet normals."""
    celltype = getattr(basix.CellType, cellname)
    out = basix.cell.facet_outward_normals(celltype)
    symbol = L.Symbol(f"{cellname}_{tablename}", dtype=L.DataType.REAL)
    return L.ArrayDecl(symbol, values=out, const=True)


def facet_orientation(tablename, cellname):
    """Write facet orientations."""
    celltype = getattr(basix.CellType, cellname)
    out = basix.cell.facet_orientations(celltype)
    symbol = L.Symbol(f"{cellname}_{tablename}", dtype=L.DataType.REAL)
    return L.ArrayDecl(symbol, values=np.asarray(out), const=True)


def closure_dofs_tables(entity_type, integrands, parent_element):
    """Write the closure-dofs tables a mixed-dimensional integrand needs.

    Returns the table declarations for gathering any lower-dimensional
    submesh's coordinate dofs out of the integration domain's own
    `coordinate_dofs` buffer, or an empty list when the integrand has no
    such submesh geometry.
    """
    if parent_element is None:
        return []
    table_kind, expected_codim = CLOSURE_DOFS_TABLES.get(entity_type, (None, None))
    if table_kind is None:
        return []

    parent_tdim = parent_element.cell.topological_dimension
    needed = False
    for integrand in integrands:
        for attr in integrand["factorization"].nodes.values():
            mt = attr.get("mt")
            if mt is None:
                continue
            if type(mt.terminal) not in (ufl.geometry.SpatialCoordinate, ufl.geometry.Jacobian):
                continue
            domain = ufl.domain.extract_unique_domain(mt.terminal)
            # Match on the exact codimension the table is built for, so
            # this agrees with the check in
            # `definitions._define_coordinate_dofs_lincomb` rather than
            # emitting a table that side would reject.
            if parent_tdim - domain.topological_dimension == expected_codim:
                needed = True

    if not needed:
        return []
    return [write_table(table_kind, parent_element.cell.cellname, parent_element)]


# Geometry defined on a sub-entity, and the entity a kernel must be over to
# use it. `access.py` indexes these by the kernel's local entity index,
# which is NULL in a cell kernel and the wrong kind of index in any other.
_ENTITY_GEOMETRY: tuple[tuple[type, str], ...] = (
    (ufl.geometry.GeometricFacetQuantity, "facet"),
    (ufl.geometry.GeometricRidgeQuantity, "ridge"),
)

# The geometry quantity a kernel can reference, and the table it is read
# from. These are the names `write_table` knows and the ones `access.py`
# emits symbols for; form and expression kernels share both.
_GEOMETRY_TABLES: dict[type, str] = {
    ufl.geometry.FacetEdgeVectors: "facet_edge_vertices",
    ufl.geometry.CellFacetJacobian: "cell_facet_jacobian",
    ufl.geometry.CellRidgeJacobian: "cell_ridge_jacobian",
    ufl.geometry.ReferenceCellVolume: "reference_cell_volume",
    ufl.geometry.ReferenceFacetVolume: "reference_facet_volume",
    ufl.geometry.ReferenceCellEdgeVectors: "reference_cell_edge_vectors",
    ufl.geometry.ReferenceFacetEdgeVectors: "reference_facet_edge_vectors",
    ufl.geometry.ReferenceNormal: "reference_normals",
    ufl.geometry.FacetOrientation: "facet_orientation",
}


def static_tables(
    entity_type: entity_types,
    integrands: Iterable[IntermediateIntegrandIR],
    parent_element: basix.ufl._ElementBase | None,
) -> list[L.ArrayDecl]:
    """Write the static tables of geometry data a kernel needs.

    The single entry point for both form and expression kernels: they read
    the same geometry quantities through the same tables, so keeping one
    copy of this mapping is what stops the two drifting apart.

    Args:
        entity_type: The entity the kernel is over.
        integrands: The kernel's integrands, whose factorization graphs are
            scanned for the geometry quantities actually referenced.
        parent_element: Coordinate element of the integration domain, used
            for the closure-dofs tables. May be None.

    Returns:
        Table declarations, in a deterministic order.
    """
    integrands = list(integrands)

    # One entry per quantity, created up front so tables are emitted in
    # `_GEOMETRY_TABLES` order rather than the order the graph is walked.
    # Quantities the integrands never reference keep an empty set.
    cellnames: dict[type, set[str]] = {terminal: set() for terminal in _GEOMETRY_TABLES}
    for integrand in integrands:
        for attr in integrand["factorization"].nodes.values():
            mt = attr.get("mt")
            if mt is None:
                continue
            terminal = type(mt.terminal)
            # UFL checks facet quantities in forms, but not ridge
            # quantities, and nothing in expressions, which have no measure
            for base, required in _ENTITY_GEOMETRY:
                if issubclass(terminal, base) and entity_type != required:
                    raise RuntimeError(
                        f"{terminal.__name__} is only defined on a {required}, "
                        f"but this kernel is over a {entity_type}."
                    )
            if terminal in _GEOMETRY_TABLES:
                ud = ufl.domain.extract_unique_domain(mt.terminal)
                assert ud is not None and isinstance(ud, ufl.Mesh)
                cellnames[terminal].add(ud.ufl_cell().cellname)

    tables: list[L.ArrayDecl] = []
    for terminal, names in cellnames.items():
        tablename = _GEOMETRY_TABLES[terminal]
        for cellname in sorted(names):
            tables.append(write_table(tablename, cellname))

    tables.extend(closure_dofs_tables(entity_type, integrands, parent_element))
    return tables
