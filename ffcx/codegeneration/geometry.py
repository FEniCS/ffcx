# Copyright (C) 2021 Matthew Scroggs
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Geometry."""

import basix
import numpy as np
import ufl

import ffcx.codegeneration.lnodes as L

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

# Number of possible facet/ridge orientations ("quadrature_permutation"
# values) FFCx/DOLFINx currently support permuting for -- matches the
# variant counts `ffcx/ir/elementtables.py`'s `build_optimized_tables`
# already permutes mixed-dimensional-submesh element tables over for the
# same cell types (interior_facet / mixed-dim codim-0 case).
_FACET_NPERM = {
    "triangle": 2,  # interval facets
    "quadrilateral": 2,  # interval facets
    "tetrahedron": 6,  # triangle facets
    "hexahedron": 8,  # quadrilateral facets
}
_RIDGE_NPERM = {
    # Every 3D cell's ridges (codim-2 entities) are edges.
    "tetrahedron": 2,
    "hexahedron": 2,
    "prism": 2,
    "pyramid": 2,
}


def write_table(tablename, cellname, coordinate_element=None):
    """Write a table."""
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
    if sub_elements:
        (scalar_element,) = set(sub_elements)
    else:
        scalar_element = coordinate_element
    return scalar_element.basix_element


def _closure_dofs_table(tablename, cellname, coordinate_element, entity_dim, nperm_by_cell):
    """Write a per-orientation table of a cell's own sub-entity closure dofs.

    One row per possible sub-entity orientation ("quadrature_permutation"
    value), one column per sub-entity of dimension `entity_dim`, listing
    the parent cell's own scalar coordinate dof indices in that
    sub-entity's closure.

    A co-dimensional entity's coordinate dofs are always a subset of its
    parent cell. The permutation handles the fact that the submesh's own
    dofmap is always canonically (globally) oriented, but the parent
    cell's local sub-entity closure is not. We use the
    `quadrature_permutation` to resolve this mismatch.
    """
    celltype = getattr(basix.CellType, cellname)
    nperm = nperm_by_cell.get(cellname)
    if nperm is None:
        raise NotImplementedError(
            f"Mixed-dimensional integrals with a submesh domain are not supported for "
            f"cell type {cellname!r} (entity dimension {entity_dim})."
        )
    be = _scalar_basix_element(coordinate_element)
    if not be.dof_transformations_are_permutations:
        raise NotImplementedError(
            "Mixed-dimensional coordinate-dofs gathering requires a coordinate element whose "
            "dof transformations are permutations."
        )
    topology = basix.topology(celltype)
    (entity_celltype,) = set(basix.cell.subentity_types(celltype)[entity_dim])
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


def facet_closure_dofs(tablename, cellname, coordinate_element):
    """Write a facet-closure-dofs table (see `_closure_dofs_table`)."""
    celltype = getattr(basix.CellType, cellname)
    tdim = len(basix.topology(celltype)) - 1
    return _closure_dofs_table(tablename, cellname, coordinate_element, tdim - 1, _FACET_NPERM)


def ridge_closure_dofs(tablename, cellname, coordinate_element):
    """Write a ridge-closure-dofs table (see `_closure_dofs_table`)."""
    celltype = getattr(basix.CellType, cellname)
    tdim = len(basix.topology(celltype)) - 1
    if tdim - 2 == 0:
        # A 2D cell's ridges are its vertices, with no orientation to
        # permute over. Keeps the `ridge_closure_dofs` name, since the
        # symbol is derived from entity_type, not from entity dimension.
        return _vertex_closure_dofs(tablename, cellname, coordinate_element)
    return _closure_dofs_table(tablename, cellname, coordinate_element, tdim - 2, _RIDGE_NPERM)


def _vertex_closure_dofs(tablename, cellname, coordinate_element):
    """Write a single-row table of each of a cell's vertices' own closure dofs.

    A vertex has no orientation ambiguity, unlike `facet_closure_dofs`
    and `ridge_closure_dofs` on a 3D cell. Therefore this table has
    a single row and `quadrature_permutation` is never consulted to index it.
    """
    celltype = getattr(basix.CellType, cellname)
    be = _scalar_basix_element(coordinate_element)
    num_vertices = len(basix.topology(celltype)[0])
    closure_dofs = be.entity_closure_dofs[0]

    if any(len(closure_dofs[v]) != 1 for v in range(num_vertices)):
        raise NotImplementedError(
            "Mixed-dimensional coordinate-dofs gathering is only implemented for coordinate "
            "elements whose vertex closure is a single dof."
        )

    out = np.array([[closure_dofs[v] for v in range(num_vertices)]], dtype=int)
    symbol = L.Symbol(f"{cellname}_{tablename}", dtype=L.DataType.INT)
    return L.ArrayDecl(symbol, values=out, const=True)


def peak_closure_dofs(tablename, cellname, coordinate_element):
    """Write a peak-closure-dofs table (peak = a cell's own vertex)."""
    return _vertex_closure_dofs(tablename, cellname, coordinate_element)


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
