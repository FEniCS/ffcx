# Copyright (C) 2011-2023 Martin Sandve Alnæs, Igor A. Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""FFCx/UFC specific variable definitions."""

import logging

import ffcx.codegeneration.lnodes as L
from ffcx.codegeneration.backend import FFCXBackend
from ffcx.ir.elementtables import UniqueTableReferenceT
from ffcx.ir.representationutils import QuadratureRule
from ffcx.ir.analysis.modified_terminals import ModifiedTerminal
from typing import List, Union
import ufl


logger = logging.getLogger("ffcx")


def create_quadrature_index(quadrature_rule: QuadratureRule, quadrature_index_symbol):
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


def get(backend: FFCXBackend, mt: ModifiedTerminal, tabledata: UniqueTableReferenceT,
        quadrature_rule: QuadratureRule, access: L.Symbol) -> Union[L.Section, List]:
    """Return definition code for a terminal."""
    # Call appropriate handler, depending on the type of terminal
    terminal = mt.terminal
    ttype = type(terminal)

    # Look for parent class of ttype or direct handler
    while ttype not in handler_lookup and ttype.__bases__:  # noqa: E721
        ttype = ttype.__bases__[0]

    # Get the handler from the lookup, or None if not found
    handler = handler_lookup.get(ttype)

    if handler is None:
        raise NotImplementedError(f"No handler for terminal type: {ttype}")

    # Call the handler
    return handler(backend, mt, tabledata, quadrature_rule, access)


def coefficient(backend, mt: ModifiedTerminal, tabledata: UniqueTableReferenceT,
                quadrature_rule: QuadratureRule, symbol: L.Symbol) -> Union[L.Section, List]:
    """Return definition code for coefficients."""
    # For applying tensor product to coefficients, we need to know if the coefficient
    # has a tensor factorisation and if the quadrature rule has a tensor factorisation.
    # If both are true, we can apply the tensor product to the coefficient.

    iq_symbol = backend.symbols.quadrature_loop_index
    ic_symbol = backend.symbols.coefficient_dof_sum_index

    iq = create_quadrature_index(quadrature_rule, iq_symbol)
    ic = create_dof_index(tabledata, ic_symbol)
    symbol.shape = tuple(iq.sizes)

    # Get properties of tables
    ttype = tabledata.ttype
    num_dofs = tabledata.values.shape[3]
    bs = tabledata.block_size
    begin = tabledata.offset
    end = begin + bs * (num_dofs - 1) + 1

    if ttype == "zeros":
        logging.debug("Not expecting zero coefficients to get this far.")
        return []

    # For a constant coefficient we reference the dofs directly, so no definition needed
    if ttype == "ones" and end - begin == 1:
        return []

    assert begin < end

    entitytype = backend.ir.entitytype

    # Get access to element table
    FE, tables = backend.access.table_access(tabledata, entitytype, mt.restriction, iq, ic)
    dof_access: L.ArrayAccess = backend.symbols.coefficient_dof_access(mt.terminal, (ic.global_index) * bs + begin)

    declaration: List[L.Declaration] = [L.declaration(symbol, [0.0])]
    access = symbol[iq] if symbol.size() else symbol

    body = [L.AssignAdd(access, dof_access * FE)]
    code = [L.create_nested_for_loops([ic], body)]

    name = type(mt.terminal).__name__
    input = [dof_access.array, *tables]
    output = [symbol]
    annotations = [L.Annotation.fuse] if not backend.ir.sum_factorization else []

    # assert input and output are Symbol objects
    assert all(isinstance(i, L.Symbol) for i in input)
    assert all(isinstance(o, L.Symbol) for o in output)

    return L.Section(name, code, declaration, input, output, annotations)


def _define_coordinate_dofs_lincomb(backend, mt: ModifiedTerminal, tabledata: UniqueTableReferenceT,
                                    quadrature_rule: QuadratureRule, symbol: L.Symbol) -> Union[L.Section, List]:
    """Define x or J as a linear combination of coordinate dofs with given table data."""
    # Get properties of domain
    domain = ufl.domain.extract_unique_domain(mt.terminal)
    coordinate_element = domain.ufl_coordinate_element()
    num_scalar_dofs = coordinate_element._sub_element.dim

    num_dofs = tabledata.values.shape[3]
    begin = tabledata.offset

    assert num_scalar_dofs == num_dofs

    # Find table name
    ttype = tabledata.ttype

    assert ttype != "zeros"
    assert ttype != "ones"

    # Get access to element table
    ic_symbol = backend.symbols.coefficient_dof_sum_index
    iq_symbol = backend.symbols.quadrature_loop_index
    ic = create_dof_index(tabledata, ic_symbol)
    iq = create_quadrature_index(quadrature_rule, iq_symbol)
    FE, tables = backend.access.table_access(tabledata, backend.ir.entitytype, mt.restriction, iq, ic)
    symbol.shape = tuple(iq.sizes)

    dof_access = L.Symbol("coordinate_dofs", dtype=L.DataType.REAL)

    # coordinate dofs is always 3d
    dim = 3
    offset = 0
    if mt.restriction == "-":
        offset = num_scalar_dofs * dim

    code = []
    declaration = [L.declaration(symbol, [0.0])]
    access = symbol[iq] if symbol.size() else symbol

    body = [L.AssignAdd(access, dof_access[ic.global_index * dim + begin + offset] * FE)]
    code = [L.create_nested_for_loops([ic], body)]

    name = type(mt.terminal).__name__
    output = [symbol]
    input = [dof_access, *tables]
    annotations = [L.Annotation.fuse] if not backend.ir.sum_factorization else []

    # assert input and output are Symbol objects
    assert all(isinstance(i, L.Symbol) for i in input)
    assert all(isinstance(o, L.Symbol) for o in output)

    return L.Section(name, code, declaration, input, output, annotations)


def spatial_coordinate(backend, mt: ModifiedTerminal, tabledata: UniqueTableReferenceT,
                       quadrature_rule: QuadratureRule, access: L.Symbol) -> Union[L.Section, List]:
    """Return definition code for the physical spatial coordinates.

    If physical coordinates are given:
        No definition needed.

    If reference coordinates are given:
        x = sum_k xdof_k xphi_k(X)

    If reference facet coordinates are given:
        x = sum_k xdof_k xphi_k(Xf)
    """
    integral_type = backend.ir.integral_type
    if integral_type in ufl.custom_integral_types:
        # FIXME: Jacobian may need adjustment for custom_integral_types
        if mt.local_derivatives:
            logging.exception("FIXME: Jacobian in custom integrals is not implemented.")
        return []
    else:
        return _define_coordinate_dofs_lincomb(backend, mt, tabledata, quadrature_rule, access)


def jacobian(backend, mt: ModifiedTerminal, tabledata: UniqueTableReferenceT,
             quadrature_rule: QuadratureRule, access: L.Symbol) -> Union[L.Section, List]:
    """Return definition code for the Jacobian of x(X)."""
    return _define_coordinate_dofs_lincomb(backend, mt, tabledata, quadrature_rule, access)


def pass_through(backend, mt: ModifiedTerminal, tabledata: UniqueTableReferenceT,
                 quadrature_rule: QuadratureRule, access: L.Symbol) -> Union[L.Section, List]:
    """Return definition code for pass through terminals."""
    return []


# called, depending on the first argument type.
handler_lookup = {ufl.coefficient.Coefficient: coefficient,
                  ufl.geometry.Jacobian: _define_coordinate_dofs_lincomb,
                  ufl.geometry.SpatialCoordinate: spatial_coordinate,
                  ufl.constant.Constant: pass_through,
                  ufl.geometry.CellVertices: pass_through,
                  ufl.geometry.FacetEdgeVectors: pass_through,
                  ufl.geometry.CellEdgeVectors: pass_through,
                  ufl.geometry.CellFacetJacobian: pass_through,
                  ufl.geometry.ReferenceCellVolume: pass_through,
                  ufl.geometry.ReferenceFacetVolume: pass_through,
                  ufl.geometry.ReferenceCellEdgeVectors: pass_through,
                  ufl.geometry.ReferenceFacetEdgeVectors: pass_through,
                  ufl.geometry.ReferenceNormal: pass_through,
                  ufl.geometry.CellOrientation: pass_through,
                  ufl.geometry.FacetOrientation: pass_through}
