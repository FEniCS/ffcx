import numpy
import ufl
import basix


basix_cells = {
    "interval": basix.CellType.interval,
    "triangle": basix.CellType.triangle,
    "tetrahedron": basix.CellType.tetrahedron,
    "quadrilateral": basix.CellType.quadrilateral,
    "hexahedron": basix.CellType.hexahedron
}

# This dictionary can be used to map ufl element names to basix element names.
# Currently all the names agree but this will not necessarily remian true.
ufl_to_basix_names = {
    "Q": "Lagrange",
    "DQ": "Discontinuous Lagrange",
    "RTCE": "Nedelec 1st kind H(curl)",
    "NCE": "Nedelec 1st kind H(curl)",
    "RTCF": "Raviart-Thomas",
    "NCF": "Raviart-Thomas",
}


def create_element(ufl_element):
    # TODO: EnrichedElement
    # TODO: Short/alternative names for elements

    if isinstance(ufl_element, ufl.VectorElement):
        return BlockedElement(create_element(ufl_element.sub_elements()[0]),
                              ufl_element.num_sub_elements())
    if isinstance(ufl_element, ufl.TensorElement):
        return BlockedElement(create_element(ufl_element.sub_elements()[0]),
                              ufl_element.num_sub_elements(), None)  # TODO: block shape

    if isinstance(ufl_element, ufl.MixedElement):
        return MixedElement([create_element(e) for e in ufl_element.sub_elements()])

    if ufl_element.family() in ufl_to_basix_names:
        return BasixElement(basix.create_element(
            ufl_to_basix_names[ufl_element.family()], ufl_element.cell().cellname(), ufl_element.degree()))

    if ufl_element.family() == "Quadrature":
        return QuadratureElement(ufl_element)

    return BasixElement(basix.create_element(
        ufl_element.family(), ufl_element.cell().cellname(), ufl_element.degree()))


def basix_index(*args):
    return basix.index(*args)


def create_quadrature(cellname, degree, rule):
    if cellname == "vertex":
        return [[]], [1]
    return basix.make_quadrature(rule, basix_cells[cellname], degree)


def reference_cell_vertices(cellname):
    return basix.geometry(basix_cells[cellname])


def map_facet_points(points, facet, cellname):
    geom = basix.geometry(basix_cells[cellname])
    facet_vertices = [geom[i] for i in basix.topology(basix_cells[cellname])[-2][facet]]

    return [facet_vertices[0] + sum((i - facet_vertices[0]) * j for i, j in zip(facet_vertices[1:], p))
            for p in points]


class BaseElement:
    def tabulate(self, nderivs, points):
        raise NotImplementedError

    @property
    def base_transformations(self):
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
    def is_blocked(self):
        raise NotImplementedError


class BasixElement(BaseElement):
    def __init__(self, element):
        self.element = element

    def tabulate(self, nderivs, points):
        return self.element.tabulate(nderivs, points)

    @property
    def base_transformations(self):
        return self.element.base_transformations()

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
        # TODO: move this to basix, then remove this wrapper class
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
        return basix.family_to_str(self.element.family)

    @property
    def reference_topology(self):
        return basix.topology(self.element.cell_type)

    @property
    def reference_geometry(self):
        return basix.geometry(self.element.cell_type)

    @property
    def is_blocked(self):
        return False


class MixedElement(BaseElement):
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

    @property
    def base_transformations(self):
        for e in self.sub_elements[1:]:
            assert len(e.base_transformations) == len(self.sub_elements[0].base_transformations)
        transformations = [[] for i in self.sub_elements[0].base_transformations]
        for e in self.sub_elements:
            for i, b in enumerate(e.base_transformations):
                transformations[i].append(b)

        output = []
        for p in transformations:
            new_transformation = numpy.zeros((sum(i.shape[0] for i in p), sum(i.shape[1] for i in p)))
            row_start = 0
            col_start = 0
            for i in p:
                new_transformation[row_start: row_start + i.shape[0], col_start: col_start + i.shape[1]] = i
                row_start += i.shape[0]
                col_start += i.shape[1]
            output.append(new_transformation)
        return output

    @property
    def interpolation_matrix(self):
        try:
            matrix = numpy.zeros((self.dim, len(self.points) * self.value_size))
            start_row = 0
            start_col = 0
            for e in self.sub_elements:
                m = e.interpolation_matrix
                matrix[start_row: start_row + m.shape[0], start_col: start_col + m.shape[1]] = m
                start_row += m.shape[0]
                start_col += m.shape[1]
            return matrix
        except ValueError:
            return numpy.zeros((0, 0))

    @property
    def points(self):
        try:
            return numpy.vstack([e.points for e in self.sub_elements])
        except ValueError:
            return numpy.zeros(0)

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
        # Return empty matrix to disable interpolation into mixed elements
        return numpy.zeros(0, 0)

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
    def is_blocked(self):
        return False


class BlockedElement(BaseElement):
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
            new_table = numpy.zeros((table.shape[0], table.shape[1] * self.block_size**2))
            for block in range(self.block_size):
                col = block * (self.block_size + 1)
                new_table[:, col: col + table.shape[1] * self.block_size**2: self.block_size**2] = table
            output.append(new_table)
        return output

    @property
    def is_blocked(self):
        return True

    @property
    def base_transformations(self):
        assert len(self.block_shape) == 1  # TODO: block shape

        output = []
        for transformation in self.sub_element.base_transformations:
            new_transformation = numpy.zeros((transformation.shape[0] * self.block_size,
                                              transformation.shape[1] * self.block_size))
            for i in range(self.block_size):
                new_transformation[i::self.block_size, i::self.block_size] = transformation
            output.append(new_transformation)
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
        return [[j * self.block_size for j in i] for i in self.sub_element.entity_dofs]

    @property
    def entity_dof_numbers(self):
        # TODO: should this return this, or should it take blocks into account?
        return [[[k * self.block_size + b for k in j for b in range(self.block_size)]
                 for j in i] for i in self.sub_element.entity_dof_numbers]

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

    @property
    def num_reference_components(self):
        return self.sub_element.num_reference_components


class QuadratureElement(BaseElement):
    def __init__(self, ufl_element):
        self._points, _ = create_quadrature(ufl_element.cell().cellname(),
                                            ufl_element.degree(), ufl_element.quadrature_scheme())
        self._ufl_element = ufl_element

    def tabulate(self, nderivs, points):
        if nderivs > 0:
            raise ValueError("Cannot take derivatives of Quadrature element.")

        if points.shape != self.points.shape:
            raise ValueError("Mismatch of tabulation points and element points.")
        tables = [numpy.eye(points.shape[0], points.shape[0])]
        return tables

    @property
    def base_transformations(self):
        perm_count = 0
        for i in range(1, self._ufl_element.cell().topological_dimension()):
            if i == 1:
                perm_count += self._ufl_element.cell().num_edges()
            if i == 2:
                perm_count += self._ufl_element.cell().num_facets() * 2

        return [numpy.identity(self.dim) for i in range(perm_count)]

    @property
    def interpolation_matrix(self):
        return numpy.eye(self.points.shape[0], self.points.shape[1])

    @property
    def points(self):
        return self._points

    @property
    def dim(self):
        return self.points.shape[0]

    @property
    def value_size(self):
        return 1

    @property
    def value_shape(self):
        return [1]

    @property
    def entity_dofs(self):
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
    def entity_dof_numbers(self):
        # TODO: move this to basix, then remove this wrapper class
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
        raise NotImplementedError

    @property
    def num_global_support_dofs(self):
        return 0

    @property
    def family_name(self):
        return self._ufl_element.family()

    @property
    def dof_mappings(self):
        assert self._ufl_element.mapping() == "identity"
        return ["affine"] * self.dim

    @property
    def num_reference_components(self):
        return {"affine": self.value_size}

    @property
    def is_blocked(self):
        return False
