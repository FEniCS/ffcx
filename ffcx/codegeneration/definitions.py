# Copyright (C) 2011-2023 Martin Sandve Alnæs, Igor A. Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""FFCx/UFCx specific variable definitions."""

import logging

import ufl

import ffcx.codegeneration.lnodes as L
from ffcx.codegeneration import geometry
from ffcx.definitions import entity_types
from ffcx.ir.analysis.modified_terminals import ModifiedTerminal
from ffcx.ir.elementtables import UniqueTableReferenceT
from ffcx.ir.representationutils import QuadratureRule

logger = logging.getLogger("ffcx")


# Cap on unrolling, so a curved high-order simplex doesn't blow up compile
# time. Tensor-product geometry always keeps the loop path regardless.
_MAX_UNROLLED_COORDINATE_DOFS = 16


def _should_unroll_coordinate_dofs(num_dofs: int, has_tensor_factorisation: bool) -> bool:
    """Return whether a coordinate-dof linear combination should be unrolled."""
    return not has_tensor_factorisation and num_dofs <= _MAX_UNROLLED_COORDINATE_DOFS


def create_quadrature_index(quadrature_rule, quadrature_index_symbol):
    """Create a multi index for the quadrature loop."""
    ranges = [0]
    name = quadrature_index_symbol.name
    indices = [L.Symbol(name, dtype=L.DataType.INT)]
    if quadrature_rule:
        ranges = [quadrature_rule.weights.size]
        if quadrature_rule.has_tensor_factors:
            dim = len(quadrature_rule.tensor_factors)
            ranges = [factor[1].size for factor in quadrature_rule.tensor_factors]
            indices = [L.Symbol(name + f"{i}", dtype=L.DataType.INT) for i in range(dim)]

    return L.MultiIndex(indices, ranges)


def create_dof_index(tabledata, dof_index_symbol):
    """Create a multi index for the coefficient dofs."""
    name = dof_index_symbol.name
    if tabledata.has_tensor_factorisation:
        dim = len(tabledata.tensor_factors)
        ranges = [factor.values.shape[-1] for factor in tabledata.tensor_factors]
        indices = [L.Symbol(f"{name}{i}", dtype=L.DataType.INT) for i in range(dim)]
    else:
        ranges = [tabledata.values.shape[-1]]
        indices = [L.Symbol(name, dtype=L.DataType.INT)]

    return L.MultiIndex(indices, ranges)


class FFCXBackendDefinitions:
    """FFCx specific code definitions."""

    entity_type: entity_types

    def __init__(self, entity_type: entity_types, integral_type: str, access, options):
        """Initialise."""
        # Store ir and options
        self.integral_type = integral_type
        self.entity_type = entity_type
        self.access = access
        self.options = options

        # called, depending on the first argument type.
        self.handler_lookup = {
            ufl.coefficient.Coefficient: self.coefficient,
            ufl.geometry.Jacobian: self._define_coordinate_dofs_lincomb,
            ufl.geometry.SpatialCoordinate: self.spatial_coordinate,
            ufl.constant.Constant: self.pass_through,
            ufl.geometry.CellVertices: self.pass_through,
            ufl.geometry.FacetEdgeVectors: self.pass_through,
            ufl.geometry.CellEdgeVectors: self.pass_through,
            ufl.geometry.CellFacetJacobian: self.pass_through,
            ufl.geometry.CellRidgeJacobian: self.pass_through,
            ufl.geometry.ReferenceCellVolume: self.pass_through,
            ufl.geometry.ReferenceFacetVolume: self.pass_through,
            ufl.geometry.ReferenceCellEdgeVectors: self.pass_through,
            ufl.geometry.ReferenceFacetEdgeVectors: self.pass_through,
            ufl.geometry.ReferenceNormal: self.pass_through,
            ufl.geometry.CellOrientation: self.pass_through,
            ufl.geometry.FacetOrientation: self.pass_through,
        }

    @property
    def symbols(self):
        """Return formatter."""
        return self.access.symbols

    def get(
        self,
        mt: ModifiedTerminal,
        tabledata: UniqueTableReferenceT,
        quadrature_rule: QuadratureRule,
        access: L.Symbol,
    ) -> L.Section | list:
        """Return definition code for a terminal."""
        # Call appropriate handler, depending on the type of terminal
        terminal = mt.terminal
        ttype = type(terminal)

        # Look for parent class of ttype or direct handler
        while ttype not in self.handler_lookup and ttype.__bases__:
            ttype = ttype.__bases__[0]

        # Get the handler from the lookup, or None if not found
        handler = self.handler_lookup.get(ttype)  # type: ignore

        if handler is None:
            raise NotImplementedError(f"No handler for terminal type: {ttype}")

        # Call the handler
        return handler(mt, tabledata, quadrature_rule, access)  # type: ignore

    def coefficient(
        self,
        mt: ModifiedTerminal,
        tabledata: UniqueTableReferenceT,
        quadrature_rule: QuadratureRule,
        access: L.Symbol,
    ) -> L.Section | list:
        """Return definition code for coefficients."""
        # For applying tensor product to coefficients, we need to know
        # if the coefficient has a tensor factorisation and if the
        # quadrature rule has a tensor factorisation. If both are true,
        # we can apply the tensor product to the coefficient.

        iq_symbol = self.symbols.quadrature_loop_index
        ic_symbol = self.symbols.coefficient_dof_sum_index

        iq = create_quadrature_index(quadrature_rule, iq_symbol)
        ic = create_dof_index(tabledata, ic_symbol)

        # Get properties of tables
        ttype = tabledata.ttype
        num_dofs = tabledata.values.shape[3]
        bs = tabledata.block_size
        begin = tabledata.offset
        assert bs is not None
        assert begin is not None
        end = begin + bs * (num_dofs - 1) + 1

        if ttype == "zeros":
            logger.debug("Not expecting zero coefficients to get this far.")
            return []

        # For a constant coefficient we reference the dofs directly, so
        # no definition needed
        if ttype == "ones" and end - begin == 1:
            return []

        assert begin < end

        # Get access to element table
        FE, tables = self.access.table_access(tabledata, self.entity_type, mt.restriction, iq, ic)
        dof_access: L.ArrayAccess = self.symbols.coefficient_dof_access(
            mt.terminal, (ic.global_index) * bs + begin
        )

        declaration: list[L.Declaration] = [L.VariableDecl(access, 0.0)]
        body = [L.AssignAdd(access, dof_access * FE)]
        code = [L.create_nested_for_loops([ic], body)]

        name = type(mt.terminal).__name__
        input = [dof_access.array, *tables]
        output = [access]
        annotations = [L.Annotation.fuse]

        # assert input and output are Symbol objects
        assert all(isinstance(i, L.Symbol) for i in input)
        assert all(isinstance(o, L.Symbol) for o in output)

        return L.Section(name, code, declaration, input, output, annotations)

    def _define_coordinate_dofs_lincomb(
        self,
        mt: ModifiedTerminal,
        tabledata: UniqueTableReferenceT,
        quadrature_rule: QuadratureRule,
        access: L.Symbol,
    ) -> L.Section | list:
        """Define x or J as a linear combination of coordinate dofs with given table data."""
        # Get properties of domain
        domain = ufl.domain.extract_unique_domain(mt.terminal)
        assert isinstance(domain, ufl.Mesh)
        coordinate_element = domain.ufl_coordinate_element()
        num_scalar_dofs = coordinate_element.sub_elements[0].dim

        num_dofs = tabledata.values.shape[3]
        begin = tabledata.offset
        assert begin is not None

        assert num_scalar_dofs == num_dofs

        # Find table name
        ttype = tabledata.ttype

        assert ttype != "zeros"

        # Get access to element table
        ic_symbol = self.symbols.coefficient_dof_sum_index
        iq_symbol = self.symbols.quadrature_loop_index
        iq = create_quadrature_index(quadrature_rule, iq_symbol)

        dof_access = L.Symbol("coordinate_dofs", dtype=L.DataType.REAL)

        parent_element = self.access.integration_domain_coordinate_element

        # coordinate dofs is always 3d
        dim = 3
        offset = 0
        if mt.restriction == "-":
            # `coordinate_dofs` holds the two cells of the *integration
            # domain* back to back. If terminal is on a submesh,
            # we need the offset to get the second parent cell
            restriction_dofs = (
                num_scalar_dofs if parent_element is None else parent_element.sub_elements[0].dim
            )
            offset = restriction_dofs * dim

        # A submesh's coordinate dofs are a subset of the integration
        # domain's, so gather them through a closure-dofs table rather
        # than indexing `coordinate_dofs` directly.
        closure_table = None
        if (
            parent_element is not None
            and (codim := parent_element.cell.topological_dimension - domain.topological_dimension)
            != 0
        ):
            table_kind, expected_codim = geometry.CLOSURE_DOFS_TABLES.get(
                self.entity_type, (None, None)
            )
            if table_kind is None or codim != expected_codim:
                raise NotImplementedError(
                    "Cannot gather a mixed-dimensional submesh's coordinate dofs: coefficient "
                    f"domain has codimension {codim} relative to the integration domain, but "
                    f"the integral's entity type is {self.entity_type!r}."
                )
            parent_cellname = parent_element.cell.cellname
            closure_table = L.Symbol(f"{parent_cellname}_{table_kind}", dtype=L.DataType.INT)
            entity = self.symbols.entity(self.entity_type, mt.restriction)
            # A vertex table has one row, so the runtime permutation
            # would be meaningless and possibly out of range.
            perm = (
                L.LiteralInt(0)
                if domain.topological_dimension == 0
                else self.symbols.entity_permutation(mt.restriction)
            )
            closure_index = closure_table[perm][entity]

        # Map a submesh-local scalar dof index to its coordinate_dofs index.
        if closure_table is None:

            def _dof(local_index):
                return local_index
        else:

            def _dof(local_index):
                return closure_index[local_index]

        code: list[L.LNode]
        if ttype == "ones":
            # Point meshes have DG-0 basis functions, so the table
            # for this element has been dropped (as it is all ones).
            # A coordinate basis is a partition of unity, so an all-ones
            # table can only come from a single-dof element.
            assert num_dofs == 1
            declaration = [L.VariableDecl(access, dof_access[_dof(0) * dim + begin + offset])]
            code = []
            input = [dof_access]
        elif not _should_unroll_coordinate_dofs(num_dofs, tabledata.has_tensor_factorisation):
            # Many coordinate dofs: keep the runtime loop instead of one
            # literal term per dof.
            ic = create_dof_index(tabledata, ic_symbol)
            FE, tables = self.access.table_access(
                tabledata, self.entity_type, mt.restriction, iq, ic
            )
            code = []
            declaration = [L.VariableDecl(access, 0.0)]
            body = [
                L.AssignAdd(access, dof_access[_dof(ic.global_index) * dim + begin + offset] * FE)
            ]
            code = [L.create_nested_for_loops([ic], body)]
            input = [dof_access, *tables]
        else:
            # Few dofs: emit a single literal-indexed sum instead of a
            # runtime loop. This avoids GCC's vectoriser mis-vectorising a
            # tiny fixed-trip-count reduction into an expensive
            # permute-then-horizontal-sum sequence, and lets the table's
            # exact 0.0/+-1.0 low-order values constant-fold away.
            terms = []
            coord_tables: list[L.Symbol] = []
            for k in range(num_dofs):
                ic_k = L.MultiIndex([L.LiteralInt(k)], [num_dofs])
                FE_k, tables_k = self.access.table_access(
                    tabledata, self.entity_type, mt.restriction, iq, ic_k
                )
                for t in tables_k:
                    if t not in coord_tables:
                        coord_tables.append(t)
                terms.append(dof_access[_dof(k) * dim + begin + offset] * FE_k)
            declaration = [L.VariableDecl(access, L.Sum(terms))]
            code = []
            input = [dof_access, *coord_tables]

        if closure_table is not None:
            input.append(closure_table)

        name = type(mt.terminal).__name__
        output = [access]
        annotations = [L.Annotation.fuse]

        # assert input and output are Symbol objects
        assert all(isinstance(i, L.Symbol) for i in input)
        assert all(isinstance(o, L.Symbol) for o in output)

        return L.Section(name, code, declaration, input, output, annotations)

    def spatial_coordinate(
        self,
        mt: ModifiedTerminal,
        tabledata: UniqueTableReferenceT,
        quadrature_rule: QuadratureRule,
        access: L.Symbol,
    ) -> L.Section | list:
        """Return definition code for the physical spatial coordinates.

        If physical coordinates are given:
          No definition needed.

        If reference coordinates are given:
          x = sum_k xdof_k xphi_k(X)

        If reference facet coordinates are given:
          x = sum_k xdof_k xphi_k(Xf)
        """
        if self.integral_type in ufl.custom_integral_types:
            # FIXME: Jacobian may need adjustment for custom_integral_types
            if mt.local_derivatives:
                logger.error("FIXME: Jacobian in custom integrals is not implemented.")
            return []
        else:
            return self._define_coordinate_dofs_lincomb(mt, tabledata, quadrature_rule, access)

    def jacobian(
        self,
        mt: ModifiedTerminal,
        tabledata: UniqueTableReferenceT,
        quadrature_rule: QuadratureRule,
        access: L.Symbol,
    ) -> L.Section | list:
        """Return definition code for the Jacobian of x(X)."""
        return self._define_coordinate_dofs_lincomb(mt, tabledata, quadrature_rule, access)

    def pass_through(
        self,
        mt: ModifiedTerminal,
        tabledata: UniqueTableReferenceT,
        quadrature_rule: QuadratureRule,
        access: L.Symbol,
    ) -> L.Section | list:
        """Return definition code for pass through terminals."""
        return []
