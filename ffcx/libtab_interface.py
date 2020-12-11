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

# This dictionary can be used to map ufl element names to libtab element names.
# Currently all the names agree but this will not necessarily remian true.
ufl_to_libtab_names = {
    "Lagrange": "Lagrange"
}


def create_libtab_element(ufl_element):
    # TODO: EnrichedElement
    # TODO: Short/alternative names for elements

    if isinstance(ufl_element, ufl.VectorElement):
        return BlockedElement(create_libtab_element(ufl_element.sub_elements()[0]),
                              ufl_element.num_sub_elements())
    if isinstance(ufl_element, ufl.TensorElement):
        return BlockedElement(create_libtab_element(ufl_element.sub_elements()[0]),
                              ufl_element.num_sub_elements(), None)  # TODO: block shape

    if isinstance(ufl_element, ufl.MixedElement):
        return MixedElement([create_libtab_element(e) for e in ufl_element.sub_elements()])

    if ufl_element.family() in ufl_to_libtab_names:
        return LibtabElement(libtab.create_element(
            ufl_to_libtab_names[ufl_element.family()], ufl_element.cell().cellname(), ufl_element.degree()))

    return LibtabElement(libtab.create_element(
        ufl_element.family(), ufl_element.cell().cellname(), ufl_element.degree()))


def libtab_index(*args):
    return libtab.index(*args)


def create_quadrature(cellname, degree, rule):
    if cellname == "vertex":
        return [[]], [1]
    return libtab.make_quadrature(libtab_cells[cellname], degree)


def reference_cell_vertices(cellname):
    return libtab.geometry(libtab_cells[cellname])


def map_facet_points(points, facet, cellname):
    geom = libtab.geometry(libtab_cells[cellname])
    facet_vertices = [geom[i] for i in libtab.topology(libtab_cells[cellname])[-2][facet]]

    return [facet_vertices[0] + sum((i - facet_vertices[0]) * j for i, j in zip(facet_vertices[1:], p))
            for p in points]


class LibtabBaseElement:
    def tabulate(self, nderivs, points):
        raise NotImplementedError

    @property
    def base_permutations(self):
        raise NotImplementedError

    @property
    def interpolation_matrix(self):
        raise NotImplementedError

    @property
    def points(self):
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
    def entity_dofs(self):
        raise NotImplementedError

    @property
    def entity_dof_numbers(self):
        raise NotImplementedError

    @property
    def coeffs(self):
        raise NotImplementedError

    @property
    def num_global_support_dofs(self):
        raise NotImplementedError

    @property
    def family_name(self):
        raise NotImplementedError

    @property
    def reference_topology(self):
        raise NotImplementedError

    @property
    def reference_geometry(self):
        raise NotImplementedError

    @property
    def dof_mappings(self):
        raise NotImplementedError


class LibtabElement(LibtabBaseElement):
    def __init__(self, element):
        self.element = element

    def tabulate(self, nderivs, points):
        return self.element.tabulate(nderivs, points)

    @property
    def base_permutations(self):
        return self.element.base_permutations

    @property
    def interpolation_matrix(self):
        return self.element.interpolation_matrix

    @property
    def points(self):
        return self.element.points

    @property
    def dim(self):
        return self.element.dim

    @property
    def value_size(self):
        return self.element.value_size

    @property
    def value_shape(self):
        return self.element.value_shape

    @property
    def entity_dofs(self):
        return self.element.entity_dofs

    @property
    def entity_dof_numbers(self):
        # TODO: move this to libtab, then remove this wrapper class
        start_dof = 0
        entity_dofs = []
        for i in self.entity_dofs:
            dofs_list = []
            for j in i:
                dofs_list.append([start_dof + k for k in range(j)])
                start_dof += j
            entity_dofs.append(dofs_list)
        return entity_dofs

    @property
    def coeffs(self):
        return self.element.coeffs

    @property
    def num_global_support_dofs(self):
        # TODO
        return 0

    @property
    def family_name(self):
        return self.element.family_name

    @property
    def reference_topology(self):
        return libtab.topology(self.element.cell_type)

    @property
    def reference_geometry(self):
        return libtab.geometry(self.element.cell_type)

    @property
    def dof_mappings(self):
        return [self.element.mapping_name for i in range(self.dim)]


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
    def interpolation_matrix(self):
        matrix = numpy.zeros((self.dim, len(self.points) * self.value_size))
        start_row = 0
        start_col = 0
        for e in self.sub_elements:
            m = e.interpolation_matrix
            matrix[start_row: start_row + m.shape[0], start_col: start_col + m.shape[1]] = m
            start_row += m.shape[0]
            start_col += m.shape[1]
        return matrix

    @property
    def points(self):
        return numpy.vstack([e.points for e in self.sub_elements])

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
    def entity_dofs(self):
        data = [e.entity_dofs for e in self.sub_elements]
        return [[sum(d[tdim][entity_n] for d in data) for entity_n, _ in enumerate(entities)]
                for tdim, entities in enumerate(data[0])]

    @property
    def entity_dof_numbers(self):
        dofs = [[[] for i in entities] for entities in self.sub_elements[0].entity_dof_numbers]
        start_dof = 0
        for e in self.sub_elements:
            for tdim, entities in enumerate(e.entity_dof_numbers):
                for entity_n, entity_dofs in enumerate(entities):
                    dofs[tdim][entity_n] += [start_dof + i for i in entity_dofs]
            start_dof += e.dim
        return dofs

    @property
    def coeffs(self):
        coeff_matrix = numpy.zeros((
            sum(e.dim for e in self.sub_elements),
            max(e.coeffs.shape[1] for e in self.sub_elements)))
        start_dof = 0
        for e in self.sub_elements:
            coeff_matrix[start_dof:e.dim, :e.coeffs.shape[1]] = e.coeffs
            start_dof += e.dim
        return coeff_matrix

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
    def dof_mappings(self):
        out = []
        for e in self.sub_elements:
            out += e.dof_mappings
        return out


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
        assert self.value_size == self.block_size  # TODO: remove this assumption

        output = []
        for table in self.sub_element.tabulate(nderivs, points):
            new_table = numpy.zeros((table.shape[0], table.shape[1] * self.block_size * self.block_size))
            for i, row in enumerate(table):
                for j, item in enumerate(row):
                    new_table[i, j * self.block_size: (j + 1) * self.block_size] = [
                        item if k == j else 0 for k in range(self.block_size)]
            output.append(new_table)
        return output

    @property
    def base_permutations(self):
        assert len(self.block_shape) == 1  # TODO: block shape

        output = []
        for perm in self.sub_element.base_permutations:
            new_perm = numpy.zeros((perm.shape[0] * self.block_size, perm.shape[1] * self.block_size))
            for i in range(self.block_size):
                new_perm[i::self.block_size, i::self.block_size] = perm
            output.append(new_perm)
        return output

    @property
    def interpolation_matrix(self):
        sub_mat = self.sub_element.interpolation_matrix
        assert self.value_size == self.block_size  # TODO: remove this assumption
        mat = numpy.zeros((sub_mat.shape[0] * self.block_size, sub_mat.shape[1] * self.value_size))
        for i, row in enumerate(sub_mat):
            for j, entry in enumerate(row):
                mat[i * self.block_size: (i + 1) * self.block_size,
                    j::sub_mat.shape[1]] = entry * numpy.identity(self.block_size)
        return mat

    @property
    def points(self):
        return self.sub_element.points

    @property
    def dim(self):
        return self.sub_element.dim * self.block_size

    @property
    def sub_elements(self):
        return [self.sub_element] * self.block_size

    @property
    def value_size(self):
        return self.block_size * self.sub_element.value_size

    @property
    def value_shape(self):
        return (self.value_size, )

    @property
    def entity_dofs(self):
        return self.sub_element.entity_dofs

    @property
    def entity_dof_numbers(self):
        # TODO: should this return this, or should it take blocks into account?
        return self.sub_element.entity_dof_numbers

    @property
    def coeffs(self):
        # TODO: should this return this, or should it take blocks into account?
        return self.sub_element.coeffs

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
    def dof_mappings(self):
        return self.sub_element.dof_mappings * self.block_size
