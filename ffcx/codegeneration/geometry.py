# Copyright (C) 2021 Matthew Scroggs
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numpy as np

import basix


def write_table(L, tablename, cellname):
    if tablename == "facet_edge_vertices":
        return facet_edge_vertices(L, tablename, cellname)
    if tablename == "reference_facet_jacobian":
        return reference_facet_jacobian(L, tablename, cellname)
    if tablename == "reference_cell_volume":
        return reference_cell_volume(L, tablename, cellname)
    if tablename == "reference_facet_volume":
        return reference_facet_volume(L, tablename, cellname)
    if tablename == "reference_edge_vectors":
        return reference_edge_vectors(L, tablename, cellname)
    if tablename == "facet_reference_edge_vectors":
        return facet_reference_edge_vectors(L, tablename, cellname)
    if tablename == "reference_facet_normals":
        return reference_facet_normals(L, tablename, cellname)
    if tablename == "facet_orientation":
        return facet_orientation(L, tablename, cellname)
    raise ValueError(f"Unknown geometry table name: {tablename}")


def facet_edge_vertices(L, tablename, cellname):
    celltype = getattr(basix.CellType, cellname)
    topology = basix.topology(celltype)
    triangle_edges = basix.topology(basix.CellType.triangle)[1]
    quadrilateral_edges = basix.topology(basix.CellType.quadrilateral)[1]

    if len(topology) != 4:
        raise ValueError("Can only get facet edges for 3D cells.")

    edge_vertices = []
    for facet in topology[-2]:
        if len(facet) == 3:
            edge_vertices += [[facet[i] for i in edge] for edge in triangle_edges]
        elif len(facet) == 4:
            edge_vertices += [[facet[i] for i in edge] for edge in quadrilateral_edges]
        else:
            raise ValueError("Only triangular and quadrilateral faces supported.")

    out = np.array(edge_vertices, dtype=int)
    arr_symbol = L.Symbol(f"{cellname}_{tablename}", dtype=L.DataType.REAL)
    return L.ArrayDecl(arr_symbol, values=out, const=True)


def reference_facet_jacobian(L, tablename, cellname):
    celltype = getattr(basix.CellType, cellname)
    out = basix.cell.facet_jacobians(celltype)
    arr_symbol = L.Symbol(f"{cellname}_{tablename}", dtype=L.DataType.REAL)
    return L.ArrayDecl(arr_symbol, values=out, const=True)


def reference_cell_volume(L, tablename, cellname):
    celltype = getattr(basix.CellType, cellname)
    out = basix.cell.volume(celltype)
    symbol = L.Symbol(f"{cellname}_{tablename}", dtype=L.DataType.REAL)
    return L.VariableDecl(symbol, out)


def reference_facet_volume(L, tablename, cellname):
    celltype = getattr(basix.CellType, cellname)
    volumes = basix.cell.facet_reference_volumes(celltype)
    for i in volumes[1:]:
        if not np.isclose(i, volumes[0]):
            raise ValueError("Reference facet volume not supported for this cell type.")
    symbol = L.Symbol(f"{cellname}_{tablename}", L.DataType.REAL)
    return L.VariableDecl(symbol, volumes[0])


def reference_edge_vectors(L, tablename, cellname):
    celltype = getattr(basix.CellType, cellname)
    topology = basix.topology(celltype)
    geometry = basix.geometry(celltype)
    edge_vectors = [geometry[j] - geometry[i] for i, j in topology[1]]
    out = np.array(edge_vectors[cellname])
    arr_symbol = L.Symbol(f"{cellname}_{tablename}", dtype=L.DataType.REAL)
    return L.ArrayDecl(arr_symbol, values=out, const=True)


def facet_reference_edge_vectors(L, tablename, cellname):
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
    arr_symbol = L.Symbol(f"{cellname}_{tablename}", dtype=L.DataType.REAL)
    return L.ArrayDecl(arr_symbol, values=out, const=True)


def reference_facet_normals(L, tablename, cellname):
    celltype = getattr(basix.CellType, cellname)
    out = basix.cell.facet_outward_normals(celltype)
    arr_symbol = L.Symbol(f"{cellname}_{tablename}", dtype=L.DataType.REAL)
    return L.ArrayDecl(arr_symbol, values=out, const=True)


def facet_orientation(L, tablename, cellname):
    celltype = getattr(basix.CellType, cellname)
    out = np.array(basix.cell.facet_orientations(celltype))
    arr_symbol = L.Symbol(f"{cellname}_{tablename}", dtype=L.DataType.REAL)
    return L.ArrayDecl(arr_symbol, values=out, const=True)
