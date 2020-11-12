import numpy
import ufl
import libtab


libtab_cells = {
    "interval": libtab.CellType.interval,
    "triangle": libtab.CellType.triangle,
    "tetrahedron": libtab.CellType.tetrahedron,
    "quadrilateral": libtab.CellType.quadrilateral,
    "hexahedron": libtab.CellType.hexahedron
}


def create_libtab_element(ufl_element):
    if isinstance(ufl_element, ufl.VectorElement):
        return BlockedElement(create_libtab_element(ufl_element.sub_elements()[0]),
                              ufl_element.num_sub_elements())
    if isinstance(ufl_element, ufl.TensorElement):
        return BlockedElement(create_libtab_element(ufl_element.sub_elements()[0]),
                              ufl_element.num_sub_elements(), None)  # TODO: block shape

    if isinstance(ufl_element, ufl.MixedElement):
        return MixedElement([create_libtab_element(e) for e in ufl_element.sub_elements()])

    return libtab.create_element(
        ufl_element.family(), ufl_element.cell().cellname(), ufl_element.degree())


def create_quadrature(cellname, degree, rule):
    # TODO
    if cellname == "interval":
        return libtab.make_quadrature([[0], [1]], degree)
    if cellname == "triangle":
        return libtab.make_quadrature([[0, 0], [1, 0], [0, 1]], degree)
    if cellname == "tetrahedron":
        return libtab.make_quadrature([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], degree)
    return None, None


class LibtabBaseElement:
    def tabulate(self, nderivs, points):
        raise NotImplementedError

    @property
    def base_permutations(self):
        raise NotImplementedError

    @property
    def ndofs(self):
        raise NotImplementedError

    @property
    def value_size(self):
        raise NotImplementedError

    @property
    def value_shape(self):
        raise NotImplementedError


class MixedElement(LibtabBaseElement):
    def __init__(self, sub_elements):
        assert len(sub_elements) > 0
        self.sub_elements = sub_elements

    def tabulate(self, nderivs, points):
        results = [e.tabulate(nderivs, points) for e in self.sub_elements]
        return [numpy.hstack([r[i] for r in results]) for i, _ in enumerate(results[0])]

    @property
    def base_permutations(self):
        perms = [[] for i in self.sub_elements[0].base_permutations]
        for e in self.sub_elements:
            for i, b in enumerate(e.base_permutations):
                perms[i].append(b)

        output = []
        for p in perms:
            print(p)
            new_perm = numpy.zeros((sum(i.shape[0] for i in p), sum(i.shape[1] for i in p)))
            row_start = 0
            col_start = 0
            for i in p:
                new_perm[row_start: row_start + i.shape[0], col_start: col_start + i.shape[1]] = i
                row_start += i.shape[0]
                col_start += i.shape[1]
            output.append(new_perm)
        return output

    @property
    def ndofs(self):
        return sum(e.ndofs for e in self.sub_elements)

    @property
    def value_size(self):
        return sum(e.value_size for e in self.sub_elements)

    @property
    def value_shape(self):
        shape = tuple()
        for e in self.sub_elements:
            shape += e.value_shape
        return shape



class BlockedElement(LibtabBaseElement):
    def __init__(self, sub_element, block_size, block_shape=None):
        assert block_size > 0
        self.sub_element = sub_element
        self.block_size = block_size
        if block_shape is None:
            self.block_shape = (block_size, )
        else:
            self.block_shape = block_shape

    def tabulate(self, nderivs, points):
        assert len(self.block_shape) == 1  # TODO: block shape

        output = []
        for table in self.sub_element.tabulate(nderivs, points):
            new_table = numpy.zeros((self.block_size * table.shape[0], table.shape[1]))
            for i, row in enumerate(table):
                for j, item in enumerate(row):
                    new_table[i, j * self.block_size: (j + 1) * self.block_size] = item
            output.append(new_table)
        return output

    @property
    def base_permutations(self):
        assert len(self.block_shape) == 1  # TODO: block shape

        output = []
        for perm in self.sub_element.base_permutations:
            new_perm = numpy.zeros((perm.shape[0] * self.block_size, perm.shape[1] * self.block_size))
            for i, row in enumerate(perm):
                for j, entry in enumerate(row):
                    new_perm[i * self.block_size: (i + 1) * self.block_size,
                             j * self.block_size: (j + 1) * self.block_size] = entry
            output.append(new_perm)
        return output

    @property
    def ndofs(self):
        return self.sub_element.ndofs * self.block_size

    @property
    def sub_elements(self):
        return [self.sub_element] * self.block_size

    @property
    def value_size(self):
        return self.block_size * self.sub_element.value_size

    @property
    def value_shape(self):
        return (self.value_size, )


def reference_cell_vertices(cellname):
    return libtab.geometry(libtab_cells[cellname])


def map_facet_points(points, facet, cellname):
    geom = libtab.geometry(libtab_cells[cellname])
    facet_vertices = [geom[0] for i in libtab.topology(libtab_cells[cellname])[-2][facet]]

    return [facet_vertices[0] + sum((i - facet_vertices[0]) * j for i, j in zip(facet_vertices[1:], p))
            for p in points]
