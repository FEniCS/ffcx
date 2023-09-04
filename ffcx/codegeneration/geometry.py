# Copyright (C) 2021 Matthew Scroggs
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numpy as np
import ffcx.codegeneration.lnodes as L
import basix
import ufl


def generate_geometry_tables(integrands):
    """Generate static tables of geometry data."""
    ufl_geometry = {
        ufl.geometry.FacetEdgeVectors: facet_edge_vertices,
        ufl.geometry.CellFacetJacobian: reference_facet_jacobian,
        ufl.geometry.ReferenceCellVolume: reference_cell_volume,
        ufl.geometry.ReferenceFacetVolume: reference_facet_volume,
        ufl.geometry.ReferenceCellEdgeVectors: reference_edge_vectors,
        ufl.geometry.ReferenceFacetEdgeVectors: facet_reference_edge_vectors,
        ufl.geometry.ReferenceNormal: reference_facet_normals,
        ufl.geometry.FacetOrientation: facet_orientation,
    }
    cells = {t: set() for t in ufl_geometry.keys()}

    for integrand in integrands:
        for attr in integrand["factorization"].nodes.values():
            mt = attr.get("mt", False)
            if mt:
                t = type(mt.terminal)
                if t in ufl_geometry:
                    cells[t].add(ufl.domain.extract_unique_domain(mt.terminal).ufl_cell().cellname())

    parts = []
    for i, cell_list in cells.items():
        for c in cell_list:
            fn = ufl_geometry[i]
            parts.append(fn(c))

    return parts


def facet_edge_vertices(cellname):
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
    arr_symbol = L.Symbol(f"{cellname}_facet_edge_vertices", dtype=L.DataType.REAL)
    return L.ArrayDecl(arr_symbol, values=out, const=True)


def reference_facet_jacobian(cellname):
    celltype = getattr(basix.CellType, cellname)
    out = basix.cell.facet_jacobians(celltype)
    arr_symbol = L.Symbol(f"{cellname}_reference_facet_jacobian", dtype=L.DataType.REAL)
    return L.ArrayDecl(arr_symbol, values=out, const=True)


def reference_cell_volume(cellname):
    celltype = getattr(basix.CellType, cellname)
    out = basix.cell.volume(celltype)
    symbol = L.Symbol(f"{cellname}_reference_cell_volume", dtype=L.DataType.REAL)
    return L.VariableDecl(symbol, out)


def reference_facet_volume(cellname):
    celltype = getattr(basix.CellType, cellname)
    volumes = basix.cell.facet_reference_volumes(celltype)
    for i in volumes[1:]:
        if not np.isclose(i, volumes[0]):
            raise ValueError("Reference facet volume not supported for this cell type.")
    symbol = L.Symbol(f"{cellname}_reference_facet_volume", L.DataType.REAL)
    return L.VariableDecl(symbol, volumes[0])


def reference_edge_vectors(cellname):
    celltype = getattr(basix.CellType, cellname)
    topology = basix.topology(celltype)
    geometry = basix.geometry(celltype)
    edge_vectors = [geometry[j] - geometry[i] for i, j in topology[1]]
    out = np.array(edge_vectors[cellname])
    arr_symbol = L.Symbol(f"{cellname}_reference_edge_vectors", dtype=L.DataType.REAL)
    return L.ArrayDecl(arr_symbol, values=out, const=True)


def facet_reference_edge_vectors(cellname):
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
            edge_vectors += [geometry[facet[j]] - geometry[facet[i]] for i, j in quadrilateral_edges]
        else:
            raise ValueError("Only triangular and quadrilateral faces supported.")

    out = np.array(edge_vectors)
    arr_symbol = L.Symbol(f"{cellname}_facet_reference_edge_vectors", dtype=L.DataType.REAL)
    return L.ArrayDecl(arr_symbol, values=out, const=True)


def reference_facet_normals(cellname):
    celltype = getattr(basix.CellType, cellname)
    out = basix.cell.facet_outward_normals(celltype)
    arr_symbol = L.Symbol(f"{cellname}_reference_facet_normals", dtype=L.DataType.REAL)
    return L.ArrayDecl(arr_symbol, values=out, const=True)


def facet_orientation(cellname):
    celltype = getattr(basix.CellType, cellname)
    out = np.array(basix.cell.facet_orientations(celltype))
    arr_symbol = L.Symbol(f"{cellname}_facet_orientation", dtype=L.DataType.REAL)
    return L.ArrayDecl(arr_symbol, values=out, const=True)
