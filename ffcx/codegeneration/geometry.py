# Copyright (C) 2021 Matthew Scroggs
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import basix
import numpy


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

    out = numpy.array(edge_vertices, dtype=int)
    return L.ArrayDecl("static const unsigned int", f"{cellname}_{tablename}", out.shape, out)


def reference_facet_jacobian(L, tablename, cellname):
    celltype = getattr(basix.CellType, cellname)
    topology = basix.topology(celltype)
    geometry = basix.geometry(celltype)

    tdim = len(topology) - 1

    if tdim not in [2, 3]:
        raise ValueError("Can only get facet jacobians for 2D and 3D cells.")

    facet_jacobian = []
    if tdim == 2:
        for facet in topology[-2]:
            edge = geometry[facet[1]] - geometry[facet[0]]
            facet_jacobian.append(list(zip(edge)))

    else:
        for facet in topology[-2]:
            edge0 = geometry[facet[1]] - geometry[facet[0]]
            edge1 = geometry[facet[2]] - geometry[facet[0]]
            facet_jacobian.append(list(zip(edge0, edge1)))

    out = numpy.array(facet_jacobian)
    return L.ArrayDecl("static const double", f"{cellname}_{tablename}", out.shape, out)


def reference_cell_volume(L, tablename, cellname):
    cell_volume = {
        "interval": 1.0, "triangle": 0.5, "tetrahedron": 1 / 6,
        "quadrilateral": 1.0, "hexahedron": 1.0
    }
    out = cell_volume[cellname]
    return L.VariableDecl("static const double", f"{cellname}_{tablename}", out)


def reference_facet_volume(L, tablename, cellname):
    facet_volume = {
        "triangle": 1.0, "tetrahedron": 0.5,
        "quadrilateral": 1.0, "hexahedron": 1.0
    }
    out = facet_volume[cellname]
    return L.VariableDecl("static const double", f"{cellname}_{tablename}", out)


def reference_edge_vectors(L, tablename, cellname):
    celltype = getattr(basix.CellType, cellname)
    topology = basix.topology(celltype)
    geometry = basix.geometry(celltype)

    edge_vectors = [geometry[j] - geometry[i] for i, j in topology[1]]

    out = numpy.array(edge_vectors[cellname])
    return L.ArrayDecl("static const double", f"{cellname}_{tablename}", out.shape, out)


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

    out = numpy.array(edge_vectors)
    return L.ArrayDecl("static const double", f"{cellname}_{tablename}", out.shape, out)


def _make_normals(cellname):
    celltype = getattr(basix.CellType, cellname)
    topology = basix.topology(celltype)
    geometry = basix.geometry(celltype)
    tdim = len(topology) - 1

    midpoint = sum(geometry) / len(geometry)

    if tdim == 1:
        # Interval
        return [numpy.array([1.]), numpy.array([1.])], [-1, 1]

    normals = []
    orientations = []
    for facet in topology[-2]:
        if tdim == 2:
            # Facets are edges
            edge = geometry[facet[1]] - geometry[facet[0]]
            n = numpy.array([-edge[1], edge[0]])
        elif tdim == 3:
            # Facets are faces
            edge0 = geometry[facet[1]] - geometry[facet[0]]
            edge1 = geometry[facet[2]] - geometry[facet[0]]
            n = numpy.cross(edge0, edge1)
        else:
            raise ValueError("Normals not supported for this celltype.")
        n /= numpy.linalg.norm(n)
        normals.append(n)
        if numpy.dot(n, geometry[facet[0]] - midpoint) > 0:
            orientations.append(1)
        else:
            orientations.append(-1)
    return normals, orientations


def reference_facet_normals(L, tablename, cellname):
    normals, orientations = _make_normals(cellname)

    out = numpy.array([i * j for i, j in zip(normals, orientations)])
    return L.ArrayDecl("static const double", f"{cellname}_{tablename}", out.shape, out)


def facet_orientation(L, tablename, cellname):
    out = numpy.array(_make_normals(cellname)[1])
    return L.ArrayDecl("static const double", f"{cellname}_{tablename}", out.shape, out)
