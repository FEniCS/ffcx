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

    if tdim < 2:
        raise ValueError("Can only get facet jacobians for 2D and 3D cells.")

    facet_jacobian = []
    if tdim == 2:
        for facet in topology[-2]:
            edge = geometry[facet[1]] - geometry[facet[0]]
            facet_jacobian.append(list(zip(edge)))

    else:
        assert tdim == 3
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


def reference_facet_normals(L, tablename, cellname):
    facet_normals = {
        "interval": [[-1.0], [1.0]],
        "triangle": [[0.7071067811865476, 0.7071067811865476], [-1.0, 0.0], [0.0, -1.0]],
        "tetrahedron": [[0.5773502691896258, 0.5773502691896258, 0.5773502691896258],
                        [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
        "quadrilateral": [[0.0, -1.0], [-1.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        "hexahedron": [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0],
                       [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    }
    out = numpy.array(facet_normals[cellname])
    return L.ArrayDecl("static const double", f"{cellname}_{tablename}", out.shape, out)


def facet_orientation(L, tablename, cellname):
    facet_orientation = {
        "interval": [-1., 1.],
        "triangle": [1., -1., 1.],
        "tetrahedron": [1., -1., 1., -1.],
        "quadrilateral": [-1., 1., -1., 1.],
        "hexahedron": [-1., 1., -1., 1., -1., 1.]
    }
    out = numpy.array(facet_orientation[cellname])
    return L.ArrayDecl("static const double", f"{cellname}_{tablename}", out.shape, out)
