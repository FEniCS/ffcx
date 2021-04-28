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
    edge_vertices = {
        "tetrahedron": [[[2, 3], [1, 3], [1, 2]], [[2, 3], [0, 3], [0, 2]],
                        [[1, 3], [0, 3], [0, 1]], [[1, 2], [0, 2], [0, 1]]],
        "hexahedron": [[[0, 1], [0, 2], [1, 3], [2, 3]], [[0, 1], [0, 4], [1, 5], [4, 5]],
                       [[0, 2], [0, 4], [2, 6], [4, 6]], [[1, 3], [1, 5], [3, 7], [5, 7]],
                       [[2, 3], [2, 6], [3, 6], [3, 7]], [[4, 5], [4, 6], [5, 7], [6, 7]]]
    }
    out = numpy.array(edge_vertices[cellname], dtype=int)
    return L.ArrayDecl("static const unsigned int", f"{cellname}_{tablename}", out.shape, out)


def reference_facet_jacobian(L, tablename, cellname):
    facet_jacobian = {
        "triangle": [[[-1.0], [1.0]], [[0.0], [1.0]], [[1.0], [0.0]]],
        "tetrahedron": [[[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]],
                        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
                        [[1.0, 0.0], [0.0, 0.0], [0.0, 1.0]],
                        [[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]],
        "quadrilateral": [[[1.0], [0.0]], [[0.0], [1.0]],
                          [[0.0], [1.0]], [[1.0], [0.0]]],
        "hexahedron": [[[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]],
                       [[1.0, 0.0], [0.0, 0.0], [0.0, 1.0]],
                       [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
                       [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
                       [[1.0, 0.0], [0.0, 0.0], [0.0, 1.0]],
                       [[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]]
    }
    out = numpy.array(facet_jacobian[cellname])
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
        "interval": 1.0, "triangle": 1.0, "tetrahedron": 0.5,
        "quadrilateral": 1.0, "hexahedron": 1.0
    }
    out = facet_volume[cellname]
    return L.VariableDecl("static const double", f"{cellname}_{tablename}", out)


def reference_edge_vectors(L, tablename, cellname):
    edge_vectors = {
        "triangle": [[-1.0, 1.0], [0.0, 1.0], [1.0, 0.0]],
        "tetrahedron": [[0.0, -1.0, 1.0], [-1.0, 0.0, 1.0], [-1.0, 1.0, 0.0],
                        [0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]],
        "quadrilateral": [[1.0, 0.0], [0.0, 1.0], [0.0, 1.0], [1.0, 0.0]],
        "hexahedron": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0],
                       [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, 1.0],
                       [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]
    }
    out = numpy.array(edge_vectors[cellname])
    return L.ArrayDecl("static const double", f"{cellname}_{tablename}", out.shape, out)


def facet_reference_edge_vectors(L, tablename, cellname):
    edge_vectors = {
        "tetrahedron": [[[0.0, -1.0, 1.0], [-1.0, 0.0, 1.0], [-1.0, 1.0, 0.0]],
                        [[0.0, -1.0, 1.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]],
                        [[-1.0, 0.0, 1.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
                        [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]],
        "hexahedron": [[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]],
                       [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
                       [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]],
                       [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]],
                       [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
                       [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]]
    }
    out = numpy.array(edge_vectors[cellname])
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
