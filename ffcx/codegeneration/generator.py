# Copyright (C) 2024 Igor Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

from typing import Any, Dict, List, Set
from numbers import Integral
import ffcx.codegeneration.lnodes as L
# from ffcx.ir.analysis.graph import ExpressionGraph
from ffcx.codegeneration import geometry
from ffcx.ir.elementtables import piecewise_ttypes
from ffcx.codegeneration import definitions
from ffcx.codegeneration import access
import ufl


def extract_dtype(v, vops: List[Any]):
    """Extract dtype from ufl expression v and its operands."""
    dtypes = []
    for op in vops:
        if hasattr(op, "dtype"):
            dtypes.append(op.dtype)
        elif hasattr(op, "symbol"):
            dtypes.append(op.symbol.dtype)
        elif isinstance(op, Integral):
            dtypes.append(L.DataType.INT)
        else:
            raise RuntimeError(f"Not expecting this type of operand {type(op)}")
    is_cond = isinstance(v, ufl.classes.Condition)
    return L.DataType.BOOL if is_cond else L.merge_dtypes(dtypes)


def generate_quadrature_tables(ir, backend):
    """Generate static tables of quadrature points and weights."""
    parts = []

    # No quadrature tables for custom (given argument) or point
    # (evaluation in single vertex)
    skip = ufl.measure.point_integral_types
    if ir.integral_type in skip:
        return parts

    # Loop over quadrature rules
    for quadrature_rule, integrand in ir.integrand.items():
        # Generate quadrature weights array
        wsym = backend.symbols.weights_table(quadrature_rule)
        parts += [L.ArrayDecl(wsym, values=quadrature_rule.weights, const=True)]

    # Add leading comment if there are any tables
    parts = L.commented_code_list(parts, "Quadrature rules")
    return parts


def generate_element_tables(ir, backend):
    """Generate static tables with precomputed element basisfunction values in quadrature points."""
    parts = []
    tables = ir.unique_tables
    table_names = sorted(tables)

    for name in table_names:
        table = tables[name]
        table_symbol = L.Symbol(name, dtype=L.DataType.REAL, shape=table.shape)
        backend.symbols.element_tables[name] = table_symbol
        parts += [L.ArrayDecl(table_symbol, values=table, const=True)]

    # Add leading comment if there are any tables
    parts = L.commented_code_list(parts, [
        "Precomputed values of basis functions and precomputations",
        "FE* dimensions: [permutation][entities][points][dofs]"])
    return parts


def generate_geometry_tables(ir, backend):
    """Generate static tables of geometry data."""
    ufl_geometry = {
        ufl.geometry.FacetEdgeVectors: "facet_edge_vertices",
        ufl.geometry.CellFacetJacobian: "reference_facet_jacobian",
        ufl.geometry.ReferenceCellVolume: "reference_cell_volume",
        ufl.geometry.ReferenceFacetVolume: "reference_facet_volume",
        ufl.geometry.ReferenceCellEdgeVectors: "reference_edge_vectors",
        ufl.geometry.ReferenceFacetEdgeVectors: "facet_reference_edge_vectors",
        ufl.geometry.ReferenceNormal: "reference_facet_normals",
        ufl.geometry.FacetOrientation: "facet_orientation"
    }
    cells: Dict[Any, Set[Any]] = {t: set() for t in ufl_geometry.keys()}  # type: ignore

    for integrand in ir.integrand.values():
        for attr in integrand["factorization"].nodes.values():
            mt = attr.get("mt")
            if mt is not None:
                t = type(mt.terminal)
                if t in ufl_geometry:
                    cells[t].add(ufl.domain.extract_unique_domain(mt.terminal).ufl_cell().cellname())

    parts = []
    for i, cell_list in cells.items():
        for c in cell_list:
            parts.append(geometry.write_table(ufl_geometry[i], c))

    return parts


def generate_piecewise_partition(ir, backend, quadrature_rule):
    # Get annotated graph of factorisation
    F = ir.integrand[quadrature_rule]["factorization"]
    nodes = F.get_mode("piecewise")
    name = f"sp_{quadrature_rule.id()}"
    return generate_partition(backend, name, nodes, None)


def generate_varying_partition(ir, backend, quadrature_rule):
    # Get annotated graph of factorisation
    F = ir.integrand[quadrature_rule]["factorization"]
    nodes = F.get_mode("varying")
    name = f"sv_{quadrature_rule.id()}"
    return generate_partition(backend, name, nodes, quadrature_rule)


def generate_partition(backend, name, nodes, quadrature_rule):
    defs: List[L.LNode] = []
    intermediates: List[L.LNode] = []

    iq = backend.quadrature_indices[quadrature_rule]

    for i, attr in nodes.items():
        expr = attr['expression']
        if backend.get_var(quadrature_rule, expr):
            continue
        elif expr._ufl_is_literal_:
            lhs = L.ufl_to_lnodes(expr)
        elif (mt := attr.get('mt')):
            tabledata = attr.get('tr')
            lhs = access.get(backend, mt, tabledata, quadrature_rule)
            code = definitions.get(backend, mt, tabledata, quadrature_rule, lhs)
            defs.append(code)
        else:
            vops = [backend.get_var(quadrature_rule, op) for op in expr.ufl_operands]
            vops = [op[iq] if (hasattr(op, "size") and op.size() > 1) else op for op in vops]
            dtype = extract_dtype(expr, vops)
            rhs = L.ufl_to_lnodes(expr, *vops)
            lhs = L.Symbol(f"{name}_{len(intermediates)}", dtype=dtype)
            intermediates.append(L.VariableDecl(lhs, rhs))

        # Add to scope
        backend.set_var(quadrature_rule, expr, lhs)

    return defs, intermediates


def get_arg_factors(backend, blockdata, block_rank, quadrature_rule, iq, indices):

    ir = backend.ir
    arg_factors = []
    tables = []
    for i in range(block_rank):
        mad = blockdata.ma_data[i]
        td = mad.tabledata
        scope = ir.integrand[quadrature_rule]["modified_arguments"]
        mt = scope[mad.ma_index]
        arg_tables = []

        assert td.ttype != "zeros"

        if td.ttype == "ones":
            arg_factor = 1
        else:
            # Assuming B sparsity follows element table sparsity
            arg_factor, arg_tables = access.table_access(backend,
                                                         td, ir.entitytype, mt.restriction, iq, indices[i])

        tables += arg_tables
        arg_factors.append(arg_factor)

    return arg_factors, tables
