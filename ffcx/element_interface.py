import numpy
import ufl
import basix


basix_cells = {
    "interval": basix.CellType.interval,
    "triangle": basix.CellType.triangle,
    "tetrahedron": basix.CellType.tetrahedron,
    "prism": basix.CellType.prism,
    "pyramid": basix.CellType.pyramid,
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
    """Create an element from a UFL element."""
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
    """Get the Basix index of a derivative."""
    return basix.index(*args)


def create_quadrature(cellname, degree, rule):
    """Create a quadrature rule."""
    if cellname == "vertex":
        return [[]], [1]
    return basix.make_quadrature(rule, basix_cells[cellname], degree)


def reference_cell_vertices(cellname):
    """Get the vertices of a reference cell."""
    return basix.geometry(basix_cells[cellname])


def map_facet_points(points, facet, cellname):
    """Map points from a reference facet to a physical facet."""
    geom = basix.geometry(basix_cells[cellname])
    facet_vertices = [geom[i] for i in basix.topology(basix_cells[cellname])[-2][facet]]

    return [facet_vertices[0] + sum((i - facet_vertices[0]) * j for i, j in zip(facet_vertices[1:], p))
            for p in points]


class BaseElement:
    """An abstract element class."""

    def tabulate(self, nderivs, points):
        """Tabulate the basis functions of the element.

        Parameters
        ----------
        nderivs : int
            The number of derivatives to tabulate.
        points : np.array
            The points to tabulate at
        """
        raise NotImplementedError

    def get_component_element(self, flat_component):
        """Get an element that represents a component of the element, and the offset and stride of the component.

        For example, for a MixedElement, this will return the sub-element that represents the given component,
        the offset of that sub-element, and a stride of 1. For a BlockedElement, this will return the sub-element,
        an offset equal to the component number, and a stride equal to the block size. For vector-valued element
        (eg H(curl) and H(div) elements), this returns a ComponentElement (and as offset of 0 and a stride of 1).
        When tabulate is called on the ComponentElement, only the part of the table for the given component is
        returned.

        Parameters
        ----------
        flat_component : int
            The component

        Returns
        -------
        BaseElement
            The component element
        int
            The offset of the component
        int
            The stride of the component
        """
        raise NotImplementedError

    @property
    def element_type(self):
        """Get the element type."""
        raise NotImplementedError

    @property
    def dim(self):
        """Get the number of DOFs the element has."""
        raise NotImplementedError

    @property
    def value_size(self):
        """Get the value size of the element."""
        raise NotImplementedError

    @property
    def value_shape(self):
        """Get the value shape of the element."""
        raise NotImplementedError

    @property
    def num_entity_dofs(self):
        """Get the number of DOFs associated with each entity."""
        raise NotImplementedError

    @property
    def entity_dofs(self):
        """Get the DOF numbers associated with each entity."""
        raise NotImplementedError

    @property
    def num_entity_closure_dofs(self):
        """Get the number of DOFs associated with the closure of each entity."""
        raise NotImplementedError

    @property
    def entity_closure_dofs(self):
        """Get the DOF numbers associated with the closure of each entity."""
        raise NotImplementedError

    @property
    def num_global_support_dofs(self):
        """Get the number of globally supported DOFs."""
        raise NotImplementedError

    @property
    def family_name(self):
        """Get the family name of the element."""
        raise NotImplementedError

    @property
    def reference_topology(self):
        """Get the topology of the reference element."""
        raise NotImplementedError

    @property
    def reference_geometry(self):
        """Get the geometry of the reference element."""
        raise NotImplementedError


class BasixElement(BaseElement):
    """An element defined by Basix."""

    def __init__(self, element):
        self.element = element

    def tabulate(self, nderivs, points):
        """Tabulate the basis functions of the element.

        Parameters
        ----------
        nderivs : int
            The number of derivatives to tabulate.
        points : np.array
            The points to tabulate at
        """
        return self.element.tabulate(nderivs, points)

    def get_component_element(self, flat_component):
        """Get an element that represents a component of the element, and the offset and stride of the component.

        For vector-valued elements (eg H(curl) and H(div) elements), this returns a ComponentElement (and an
        offset of 0 and a stride of 1). When tabulate is called on the ComponentElement, only the part of the
        table for the given component is returned.

        Parameters
        ----------
        flat_component : int
            The component

        Returns
        -------
        BaseElement
            The component element
        int
            The offset of the component
        int
            The stride of the component
        """
        assert flat_component < self.value_size
        return ComponentElement(self, flat_component), 0, 1

    @property
    def element_type(self):
        """Get the element type."""
        return "ufc_basix_element"

    @property
    def dim(self):
        """Get the number of DOFs the element has."""
        return self.element.dim

    @property
    def value_size(self):
        """Get the value size of the element."""
        return self.element.value_size

    @property
    def value_shape(self):
        """Get the value shape of the element."""
        return self.element.value_shape

    @property
    def num_entity_dofs(self):
        """Get the number of DOFs associated with each entity."""
        return self.element.num_entity_dofs

    @property
    def entity_dofs(self):
        """Get the DOF numbers associated with each entity."""
        return self.element.entity_dofs

    @property
    def num_entity_closure_dofs(self):
        """Get the number of DOFs associated with the closure of each entity."""
        return self.element.num_entity_closure_dofs

    @property
    def entity_closure_dofs(self):
        """Get the DOF numbers associated with the closure of each entity."""
        return self.element.entity_closure_dofs

    @property
    def num_global_support_dofs(self):
        """Get the number of globally supported DOFs."""
        # TODO
        return 0

    @property
    def family_name(self):
        """Get the family name of the element."""
        return basix.family_to_str(self.element.family)

    @property
    def reference_topology(self):
        """Get the topology of the reference element."""
        return basix.topology(self.element.cell_type)

    @property
    def reference_geometry(self):
        """Get the geometry of the reference element."""
        return basix.geometry(self.element.cell_type)


class ComponentElement(BaseElement):
    """An element representing one component of a BasixElement."""

    def __init__(self, element, component):
        self.element = element
        self.component = component

    def tabulate(self, nderivs, points):
        """Tabulate the basis functions of the element.

        Parameters
        ----------
        nderivs : int
            The number of derivatives to tabulate.
        points : np.array
            The points to tabulate at
        """
        tables = self.element.tabulate(nderivs, points)
        output = []
        for tbl in tables:
            shape = (tbl.shape[0],) + tuple(self.element.value_shape) + (-1,)
            tbl = tbl.reshape(shape)
            if len(self.element.value_shape) == 1:
                output.append(tbl[:, self.component, :])
            elif len(self.element.value_shape) == 2:
                # TODO: Something different may need doing here if tensor is symmetric
                vs0 = self.element.value_shape[0]
                output.append(tbl[:, self.component // vs0, self.component % vs0, :])
            else:
                raise NotImplementedError
        return output

    def get_component_element(self, flat_component):
        """Get an element that represents a component of the element, and the offset and stride of the component.

        Parameters
        ----------
        flat_component : int
            The component

        Returns
        -------
        BaseElement
            The component element
        int
            The offset of the component
        int
            The stride of the component
        """
        if flat_component == 0:
            return self, 0, 1
        raise NotImplementedError


class MixedElement(BaseElement):
    """A mixed element that combines two or more elements."""

    def __init__(self, sub_elements):
        assert len(sub_elements) > 0
        self.sub_elements = sub_elements

    def tabulate(self, nderivs, points):
        """Tabulate the basis functions of the element.

        Parameters
        ----------
        nderivs : int
            The number of derivatives to tabulate.
        points : np.array
            The points to tabulate at
        """
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

    def get_component_element(self, flat_component):
        """Get an element that represents a component of the element, and the offset and stride of the component.

        For a MixedElement, this will return the sub-element that represents the given component,
        the offset of that sub-element, and a stride of 1.

        Parameters
        ----------
        flat_component : int
            The component

        Returns
        -------
        BaseElement
            The component element
        int
            The offset of the component
        int
            The stride of the component
        """
        sub_dims = [0] + [e.dim for e in self.sub_elements]
        sub_cmps = [0] + [e.value_size for e in self.sub_elements]

        irange = numpy.cumsum(sub_dims)
        crange = numpy.cumsum(sub_cmps)

        # Find index of sub element which corresponds to the current flat component
        component_element_index = numpy.where(
            crange <= flat_component)[0].shape[0] - 1

        sub_e = self.sub_elements[component_element_index]

        e, offset, stride = sub_e.get_component_element(flat_component - crange[component_element_index])
        # TODO: is this offset correct?
        return e, irange[component_element_index] + offset, stride

    @property
    def element_type(self):
        """Get the element type."""
        return "ufc_mixed_element"

    @property
    def dim(self):
        """Get the number of DOFs the element has."""
        return sum(e.dim for e in self.sub_elements)

    @property
    def value_size(self):
        """Get the value size of the element."""
        return sum(e.value_size for e in self.sub_elements)

    @property
    def value_shape(self):
        """Get the value shape of the element."""
        return (sum(e.value_size for e in self.sub_elements), )

    @property
    def num_entity_dofs(self):
        """Get the number of DOFs associated with each entity."""
        data = [e.num_entity_dofs for e in self.sub_elements]
        return [[sum(d[tdim][entity_n] for d in data) for entity_n, _ in enumerate(entities)]
                for tdim, entities in enumerate(data[0])]

    @property
    def entity_dofs(self):
        """Get the DOF numbers associated with each entity."""
        dofs = [[[] for i in entities] for entities in self.sub_elements[0].entity_dofs]
        start_dof = 0
        for e in self.sub_elements:
            for tdim, entities in enumerate(e.entity_dofs):
                for entity_n, entity_dofs in enumerate(entities):
                    dofs[tdim][entity_n] += [start_dof + i for i in entity_dofs]
            start_dof += e.dim
        return dofs

    @property
    def num_entity_closure_dofs(self):
        """Get the number of DOFs associated with the closure of each entity."""
        data = [e.num_entity_closure_dofs for e in self.sub_elements]
        return [[sum(d[tdim][entity_n] for d in data) for entity_n, _ in enumerate(entities)]
                for tdim, entities in enumerate(data[0])]

    @property
    def entity_closure_dofs(self):
        """Get the DOF numbers associated with the closure of each entity."""
        dofs = [[[] for i in entities] for entities in self.sub_elements[0].entity_closure_dofs]
        start_dof = 0
        for e in self.sub_elements:
            for tdim, entities in enumerate(e.entity_closure_dofs):
                for entity_n, entity_dofs in enumerate(entities):
                    dofs[tdim][entity_n] += [start_dof + i for i in entity_dofs]
            start_dof += e.dim
        return dofs

    @property
    def num_global_support_dofs(self):
        """Get the number of globally supported DOFs."""
        return sum(e.num_global_support_dofs for e in self.sub_elements)

    @property
    def family_name(self):
        """Get the family name of the element."""
        return "mixed element"

    @property
    def reference_topology(self):
        """Get the topology of the reference element."""
        return self.sub_elements[0].reference_topology

    @property
    def reference_geometry(self):
        """Get the geometry of the reference element."""
        return self.sub_elements[0].reference_geometry


class BlockedElement(BaseElement):
    """An element with a block size that contains multiple copies of a sub element."""

    def __init__(self, sub_element, block_size, block_shape=None):
        assert block_size > 0
        self.sub_element = sub_element
        self.block_size = block_size
        if block_shape is None:
            self.block_shape = (block_size, )
        else:
            self.block_shape = block_shape

    def tabulate(self, nderivs, points):
        """Tabulate the basis functions of the element.

        Parameters
        ----------
        nderivs : int
            The number of derivatives to tabulate.
        points : np.array
            The points to tabulate at
        """
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

    def get_component_element(self, flat_component):
        """Get an element that represents a component of the element, and the offset and stride of the component.

        For a BlockedElement, this will return the sub-element, an offset equal to the component number, and a
        stride equal to the block size.

        Parameters
        ----------
        flat_component : int
            The component

        Returns
        -------
        BaseElement
            The component element
        int
            The offset of the component
        int
            The stride of the component
        """
        return self.sub_element, flat_component, self.block_size

    @property
    def element_type(self):
        """Get the element type."""
        return "ufc_blocked_element"

    @property
    def dim(self):
        """Get the number of DOFs the element has."""
        return self.sub_element.dim * self.block_size

    @property
    def value_size(self):
        """Get the value size of the element."""
        return self.block_size * self.sub_element.value_size

    @property
    def value_shape(self):
        """Get the value shape of the element."""
        return (self.value_size, )

    @property
    def num_entity_dofs(self):
        """Get the number of DOFs associated with each entity."""
        return [[j * self.block_size for j in i] for i in self.sub_element.num_entity_dofs]

    @property
    def entity_dofs(self):
        """Get the DOF numbers associated with each entity."""
        # TODO: should this return this, or should it take blocks into account?
        return [[[k * self.block_size + b for k in j for b in range(self.block_size)]
                 for j in i] for i in self.sub_element.entity_dofs]

    @property
    def num_entity_closure_dofs(self):
        """Get the number of DOFs associated with the closure of each entity."""
        return [[j * self.block_size for j in i] for i in self.sub_element.num_entity_closure_dofs]

    @property
    def entity_closure_dofs(self):
        """Get the DOF numbers associated with the closure of each entity."""
        # TODO: should this return this, or should it take blocks into account?
        return [[[k * self.block_size + b for k in j for b in range(self.block_size)]
                 for j in i] for i in self.sub_element.entity_closure_dofs]

    @property
    def num_global_support_dofs(self):
        """Get the number of globally supported DOFs."""
        return self.sub_element.num_global_support_dofs * self.block_size

    @property
    def family_name(self):
        """Get the family name of the element."""
        return self.sub_element.family_name

    @property
    def reference_topology(self):
        """Get the topology of the reference element."""
        return self.sub_element.reference_topology

    @property
    def reference_geometry(self):
        """Get the geometry of the reference element."""
        return self.sub_element.reference_geometry


class QuadratureElement(BaseElement):
    """A quadrature element."""

    def __init__(self, ufl_element):
        self._points, _ = create_quadrature(ufl_element.cell().cellname(),
                                            ufl_element.degree(), ufl_element.quadrature_scheme())
        self._ufl_element = ufl_element

    def tabulate(self, nderivs, points):
        """Tabulate the basis functions of the element.

        Parameters
        ----------
        nderivs : int
            The number of derivatives to tabulate.
        points : np.array
            The points to tabulate at
        """
        if nderivs > 0:
            raise ValueError("Cannot take derivatives of Quadrature element.")

        if points.shape != self._points.shape:
            raise ValueError("Mismatch of tabulation points and element points.")
        tables = [numpy.eye(points.shape[0], points.shape[0])]
        return tables

    def get_component_element(self, flat_component):
        """Get an element that represents a component of the element, and the offset and stride of the component.

        Parameters
        ----------
        flat_component : int
            The component

        Returns
        -------
        BaseElement
            The component element
        int
            The offset of the component
        int
            The stride of the component
        """
        return self, 0, 1

    @property
    def element_type(self):
        """Get the element type."""
        return "ufc_quadrature_element"

    @property
    def dim(self):
        """Get the number of DOFs the element has."""
        return self._points.shape[0]

    @property
    def value_size(self):
        """Get the value size of the element."""
        return 1

    @property
    def value_shape(self):
        """Get the value shape of the element."""
        return [1]

    @property
    def num_entity_dofs(self):
        """Get the number of DOFs associated with each entity."""
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
    def entity_dofs(self):
        """Get the DOF numbers associated with each entity."""
        start_dof = 0
        entity_dofs = []
        for i in self.num_entity_dofs:
            dofs_list = []
            for j in i:
                dofs_list.append([start_dof + k for k in range(j)])
                start_dof += j
            entity_dofs.append(dofs_list)
        return entity_dofs

    @property
    def num_entity_closure_dofs(self):
        """Get the number of DOFs associated with the closure of each entity."""
        return self.num_entity_dofs

    @property
    def entity_closure_dofs(self):
        """Get the DOF numbers associated with the closure of each entity."""
        return self.entity_dofs

    @property
    def num_global_support_dofs(self):
        """Get the number of globally supported DOFs."""
        return 0

    @property
    def family_name(self):
        """Get the family name of the element."""
        return self._ufl_element.family()
