# -*- coding: utf-8 -*-
# Copyright (C) 2013-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Main algorithm for building the uflacs intermediate representation."""

import collections
import itertools
import logging

import numpy

import ufl
from ffc.ir.uflacs.analysis.factorization import compute_argument_factorization
from ffc.ir.uflacs.analysis.graph import build_scalar_graph
from ffc.ir.uflacs.analysis.modified_terminals import (analyse_modified_terminal,
                                                       is_modified_terminal)
from ffc.ir.uflacs.analysis.visualise import visualise
from ffc.ir.uflacs.elementtables import (build_optimized_tables,
                                         clamp_table_small_numbers,
                                         piecewise_ttypes)
from ufl.algorithms.balancing import balance_modifiers
from ufl.checks import is_cellwise_constant
from ufl.classes import CellCoordinate, FacetCoordinate, QuadratureWeight
from ufl.measure import (custom_integral_types, facet_integral_types,
                         point_integral_types)

logger = logging.getLogger(__name__)

ma_data_t = collections.namedtuple("ma_data_t", ["ma_index", "tabledata"])

block_data_t = collections.namedtuple("block_data_t",
                                      ["block_mode",
                                       # "safe" | "full" | "preintegrated" | "premultiplied"
                                       "ttypes",  # list of table types for each block rank
                                       "factor_index",  # int: index of factor in vertex array
                                       "factor_is_piecewise",
                                       # bool: factor is found in piecewise vertex array
                                       # instead of quadloop specific vertex array
                                       "unames",  # list of unique FE table names for each block rank
                                       "restrictions",  # restriction "+" | "-" | None for each block rank
                                       "transposed",  # block is the transpose of another
                                       "is_uniform",  # used in "preintegrated" and "premultiplied"
                                       "name",  # used in "preintegrated" and "premultiplied"
                                       "ma_data",  # used in "full", "safe" and "partial"
                                       "piecewise_ma_index"  # used in "partial"
                                       ])


def multiply_block_interior_facets(point_index, unames, ttypes, unique_tables,
                                   unique_table_num_dofs):
    rank = len(unames)
    tables = [unique_tables.get(name) for name in unames]
    num_dofs = tuple(unique_table_num_dofs[name] for name in unames)

    num_entities = max([1] + [tbl.shape[0] for tbl in tables if tbl is not None])
    ptable = numpy.zeros((num_entities, ) * rank + num_dofs)
    for facets in itertools.product(*[range(num_entities)] * rank):
        vectors = []
        for i, tbl in enumerate(tables):
            if tbl is None:
                assert ttypes[i] == "ones"
                vectors.append(numpy.ones((num_dofs[i], )))
            else:
                # Some tables are compacted along entities or points
                e = 0 if tbl.shape[0] == 1 else facets[i]
                q = 0 if tbl.shape[1] == 1 else point_index
                vectors.append(tbl[e, q, :])
        if rank > 1:
            assert rank == 2
            ptable[facets[0], facets[1], ...] = numpy.outer(*vectors)
        elif rank == 1:
            ptable[facets[0], :] = vectors[0]
        else:
            raise RuntimeError("Nothing to multiply!")

    return ptable


def multiply_block(point_index, unames, ttypes, unique_tables, unique_table_num_dofs):
    rank = len(unames)
    tables = [unique_tables.get(name) for name in unames]
    num_dofs = tuple(unique_table_num_dofs[name] for name in unames)

    num_entities = max([1] + [tbl.shape[0] for tbl in tables if tbl is not None])
    ptable = numpy.zeros((num_entities, ) + num_dofs)
    for entity in range(num_entities):
        vectors = []
        for i, tbl in enumerate(tables):
            if tbl is None:
                assert ttypes[i] == "ones"
                vectors.append(numpy.ones((num_dofs[i], )))
            else:
                # Some tables are compacted along entities or points
                e = 0 if tbl.shape[0] == 1 else entity
                q = 0 if tbl.shape[1] == 1 else point_index
                vectors.append(tbl[e, q, :])
        if rank > 1:
            ptable[entity, ...] = numpy.outer(*vectors)
        elif rank == 1:
            ptable[entity, :] = vectors[0]
        else:
            raise RuntimeError("Nothing to multiply!")

    return ptable


def integrate_block(weights, unames, ttypes, unique_tables, unique_table_num_dofs):
    tables = [unique_tables.get(name) for name in unames]
    num_dofs = tuple(unique_table_num_dofs[name] for name in unames)

    num_entities = max([1] + [tbl.shape[0] for tbl in tables if tbl is not None])
    ptable = numpy.zeros((num_entities, ) + num_dofs)
    for iq, w in enumerate(weights):
        ptable[...] += w * multiply_block(iq, unames, ttypes, unique_tables, unique_table_num_dofs)

    return ptable


def integrate_block_interior_facets(weights, unames, ttypes, unique_tables, unique_table_num_dofs):
    rank = len(unames)
    tables = [unique_tables.get(name) for name in unames]
    num_dofs = tuple(unique_table_num_dofs[name] for name in unames)

    num_entities = max([1] + [tbl.shape[0] for tbl in tables if tbl is not None])
    ptable = numpy.zeros((num_entities, ) * rank + num_dofs)
    for iq, w in enumerate(weights):
        mtable = multiply_block_interior_facets(iq, unames, ttypes, unique_tables,
                                                unique_table_num_dofs)
        ptable[...] += w * mtable

    return ptable


def uflacs_default_parameters(optimize):
    """Default parameters for tuning of uflacs code generation.

    These are considered experimental and may change without deprecation
    mechanism at any time.
    """
    p = {
        # Relative precision to use when comparing finite element table
        # values for table reuse
        "table_rtol": 1e-6,

        # Absolute precision to use when comparing finite element table
        # values for table reuse and dropping of table zeros
        "table_atol": 1e-9,

        # Point chunk size for custom integrals
        "chunk_size": 8,

        # Optimization parameters used in representation building
        # TODO: The names of these parameters can be a bit misleading
        "enable_preintegration": False,
        "enable_premultiplication": False,
        "enable_sum_factorization": False,
        "enable_block_transpose_reuse": False,
        "enable_table_zero_compression": False,

        # Code generation parameters
        "vectorize": False,
        "alignas": 0,
        "padlen": 1,
        "use_symbol_array": True,
        "tensor_init_mode": "upfront",  # interleaved | direct | upfront
    }
    if optimize:
        # Override defaults if optimization is turned on
        p.update({
            # Optimization parameters used in representation building
            # TODO: The names of these parameters can be a bit misleading
            "enable_preintegration": True,
            "enable_premultiplication": False,
            "enable_sum_factorization": True,
            "enable_block_transpose_reuse": True,
            "enable_table_zero_compression": True,

            # Code generation parameters
            "vectorize": False,
            "alignas": 32,
            "padlen": 1,
            "use_symbol_array": True,
            "tensor_init_mode": "interleaved",  # interleaved | direct | upfront
        })
    return p


def parse_uflacs_optimization_parameters(parameters, integral_type):
    """Following model from quadrature representation, extracting
    uflacs specific parameters from the global parameters dict."""

    # Get default parameters
    p = uflacs_default_parameters(optimize=True)

    # Override with uflacs specific parameters if present in given
    # global parameters dict
    for key in p:
        if key in parameters:
            value = parameters[key]
            # Casting done here because main doesn't know about these
            # parameters
            if isinstance(p[key], int):
                value = int(value)
            elif isinstance(p[key], float):
                value = float(value)
            p[key] = value

    # Conditionally disable some optimizations based on integral type,
    # i.e. these options are not valid for certain integral types
    skip_preintegrated = point_integral_types + custom_integral_types
    if integral_type in skip_preintegrated:
        p["enable_preintegration"] = False

    skip_premultiplied = point_integral_types + custom_integral_types
    if integral_type in skip_premultiplied:
        p["enable_premultiplication"] = False

    return p


def build_uflacs_ir(cell, integral_type, entitytype, integrands, tensor_shape,
                    quadrature_rules, parameters):
    # The intermediate representation dict we're building and returning
    # here
    ir = {}

    # Extract uflacs specific optimization and code generation
    # parameters
    p = parse_uflacs_optimization_parameters(parameters, integral_type)

    # Pass on parameters for consumption in code generation
    ir["params"] = p

    # Shared unique tables for all quadrature loops
    ir["unique_tables"] = {}
    ir["unique_table_types"] = {}

    # Shared piecewise expr_ir for all quadrature loops
    ir["piecewise_ir"] = {"factorization": None,
                          "modified_arguments": [],
                          "preintegrated_blocks": {},
                          "premultiplied_blocks": {},
                          "preintegrated_contributions": collections.defaultdict(list),
                          "block_contributions": collections.defaultdict(list)}

    # { num_points: expr_ir for one integrand }
    ir["varying_irs"] = {"factorization": None}

    # Whether we expect the quadrature weight to be applied or not (in
    # some cases it's just set to 1 in ufl integral scaling)
    tdim = cell.topological_dimension()
    expect_weight = (integral_type not in point_integral_types and (entitytype == "cell" or (
        entitytype == "facet" and tdim > 1) or (integral_type in custom_integral_types)))

    # Analyse each num_points/integrand separately
    assert isinstance(integrands, dict)
    all_num_points = sorted(integrands.keys())
    cases = [(num_points, [integrands[num_points]]) for num_points in all_num_points]
    ir["all_num_points"] = all_num_points

    for num_points, expressions in cases:

        assert len(expressions) == 1
        expression = expressions[0]

        # Rebalance order of nested terminal modifiers
        expression = balance_modifiers(expression)

        # Remove QuadratureWeight terminals from expression and replace with 1.0
        expression = replace_quadratureweight(expression)

        # Build initial scalar list-based graph representation
        S = build_scalar_graph(expression)
        S_targets = [i for i, v in S.nodes.items() if v.get('target', False)]
        assert len(S_targets) == 1
        S_target = S_targets[0]

        # Build terminal_data from V here before factorization. Then we
        # can use it to derive table properties for all modified
        # terminals, and then use that to rebuild the scalar graph more
        # efficiently before argument factorization. We can build
        # terminal_data again after factorization if that's necessary.

        initial_terminals = {i: analyse_modified_terminal(v['expression'])
                             for i, v in S.nodes.items()
                             if is_modified_terminal(v['expression'])}

        unique_tables, unique_table_types, unique_table_num_dofs, mt_unique_table_reference = build_optimized_tables(
            num_points,
            quadrature_rules,
            cell,
            integral_type,
            entitytype,
            initial_terminals.values(),
            ir["unique_tables"],
            p["enable_table_zero_compression"],
            rtol=p["table_rtol"],
            atol=p["table_atol"])

        # If there are any 'zero' tables, replace symbolically and rebuild graph
        if 'zeros' in unique_table_types.values():
            for i, mt in initial_terminals.items():
                # Set modified terminals with zero tables to zero
                tr = mt_unique_table_reference.get(mt)
                if tr is not None and tr.ttype == "zeros":
                    S.nodes[i]['expression'] = ufl.as_ufl(0.0)

            # Propagate expression changes using dependency list
            for i, v in S.nodes.items():
                deps = [S.nodes[j]['expression'] for j in S.out_edges[i]]
                if deps:
                    v['expression'] = v['expression']._ufl_expr_reconstruct_(*deps)

            # Rebuild scalar target expressions and graph (this may be
            # overkill and possible to optimize away if it turns out to be
            # costly)
            expression = S.nodes[S_target]['expression']

            # Rebuild scalar list-based graph representation
            S = build_scalar_graph(expression)

        # Output diagnostic graph as pdf
        if parameters['visualise']:
            visualise(S, 'S.pdf')

        # Compute factorization of arguments
        rank = len(tensor_shape)
        F = compute_argument_factorization(S, rank)

        # Get the 'target' nodes that are factors of arguments, and insert in dict
        FV_targets = [i for i, v in F.nodes.items() if v.get('target', False)]
        argument_factorization = {}
        for i in FV_targets:
            for w in F.nodes[i]['target']:
                argument_factorization[w] = i

        # Get list of indices in F which are the arguments (should be at start)
        argkeys = set()
        for w in argument_factorization:
            argkeys = argkeys | set(w)
        argkeys = list(argkeys)

        # Output diagnostic graph as pdf
        if parameters['visualise']:
            visualise(F, 'F.pdf')

        # Build set of modified_terminals for each mt factorized vertex in F
        # and attach tables, if appropriate
        for i, v in F.nodes.items():
            expr = v['expression']
            if is_modified_terminal(expr):
                mt = analyse_modified_terminal(expr)
                F.nodes[i]['mt'] = mt
                tr = mt_unique_table_reference.get(mt)
                if tr is not None:
                    F.nodes[i]['tr'] = tr

        # Attach 'status' to each node: 'inactive', 'piecewise' or 'varying'
        analyse_dependencies(F, mt_unique_table_reference)

        # Save the factorisation graph to the piecewise IR
        ir["piecewise_ir"]["factorization"] = F
        ir["piecewise_ir"]["modified_arguments"] = [F.nodes[i]['mt']
                                                    for i in argkeys]

        # Loop over factorization terms
        block_contributions = collections.defaultdict(list)
        for ma_indices, fi in sorted(argument_factorization.items()):
            # Get a bunch of information about this term
            assert rank == len(ma_indices)
            trs = tuple(F.nodes[ai]['tr'] for ai in ma_indices)

            unames = tuple(tr.name for tr in trs)
            ttypes = tuple(tr.ttype for tr in trs)
            assert not any(tt == "zeros" for tt in ttypes)

            blockmap = tuple(tr.dofmap for tr in trs)

            block_is_uniform = all(tr.is_uniform for tr in trs)

            # Collect relevant restrictions to identify blocks correctly
            # in interior facet integrals
            block_restrictions = []
            for i, ai in enumerate(ma_indices):
                if trs[i].is_uniform:
                    r = None
                else:
                    r = F.nodes[ai]['mt'].restriction
                block_restrictions.append(r)
            block_restrictions = tuple(block_restrictions)

            factor_is_piecewise = F.nodes[fi]['status'] == 'piecewise'

            # TODO: Add separate block modes for quadrature
            # Both arguments in quadrature elements
            """
            for iq
                fw = f*w
                #for i
                #    for j
                #        B[i,j] = fw*U[i]*V[j] = 0 if i != iq or j != iq
                BQ[iq] = B[iq,iq] = fw
            for (iq)
                A[iq+offset0, iq+offset1] = BQ[iq]
            """
            # One argument in quadrature element
            """
            for iq
                fw[iq] = f*w
                #for i
                #    for j
                #        B[i,j] = fw*UQ[i]*V[j] = 0 if i != iq
                for j
                    BQ[iq,j] = fw[iq]*V[iq,j]
            for (iq) for (j)
                A[iq+offset, j+offset] = BQ[iq,j]
            """

            # Decide how to handle code generation for this block
            if p["enable_preintegration"] and (factor_is_piecewise and rank > 0
                                               and "quadrature" not in ttypes):
                # - Piecewise factor is an absolute prerequisite
                # - Could work for rank 0 as well but currently doesn't
                # - Haven't considered how quadrature elements work out
                block_mode = "preintegrated"
            elif p["enable_premultiplication"] and (rank > 0 and all(tt in piecewise_ttypes
                                                                     for tt in ttypes)):
                # Integrate functional in quadloop, scale block after
                # quadloop
                block_mode = "premultiplied"
            elif p["enable_sum_factorization"]:
                if (rank == 2 and any(tt in piecewise_ttypes for tt in ttypes)):
                    # Partial computation in quadloop of f*u[i], compute
                    # (f*u[i])*v[i] outside quadloop, (or with u,v
                    # swapped)
                    block_mode = "partial"
                else:
                    # Full runtime integration of f*u[i]*v[j], can still
                    # do partial computation in quadloop of f*u[i] but
                    # must compute (f*u[i])*v[i] as well inside
                    # quadloop.  (or with u,v swapped)
                    block_mode = "full"
            else:
                # Use full runtime integration with nothing fancy going
                # on
                block_mode = "safe"

            # Carry out decision
            if block_mode == "preintegrated":
                # Add to contributions:
                # P = sum_q weight*u*v;      preintegrated here
                # B[...] = f * P[...];       generated after quadloop
                # A[blockmap] += B[...];     generated after quadloop

                cache = ir["piecewise_ir"]["preintegrated_blocks"]

                block_is_transposed = False
                pname = cache.get(unames)

                # Reuse transpose to save memory
                if p["enable_block_transpose_reuse"] and pname is None and len(unames) == 2:
                    pname = cache.get((unames[1], unames[0]))
                    if pname is not None:
                        # Cache hit on transpose
                        block_is_transposed = True

                if pname is None:
                    # Cache miss, precompute block
                    weights = quadrature_rules[num_points][1]
                    if integral_type == "interior_facet":
                        ptable = integrate_block_interior_facets(
                            weights, unames, ttypes, unique_tables, unique_table_num_dofs)
                    else:
                        ptable = integrate_block(weights, unames, ttypes, unique_tables,
                                                 unique_table_num_dofs)
                    ptable = clamp_table_small_numbers(
                        ptable, rtol=p["table_rtol"], atol=p["table_atol"])

                    pname = "PI%d" % (len(cache, ))
                    cache[unames] = pname
                    unique_tables[pname] = ptable
                    unique_table_types[pname] = "preintegrated"

                assert factor_is_piecewise
                block_unames = (pname, )
                blockdata = block_data_t(
                    block_mode, ttypes, fi, factor_is_piecewise, block_unames,
                    block_restrictions, block_is_transposed, block_is_uniform, pname,
                    None, None)
                block_is_piecewise = True

            elif block_mode == "premultiplied":
                # Add to contributions:
                # P = u*v;                        computed here
                # FI = sum_q weight * f;          generated inside quadloop
                # B[...] = FI * P[...];           generated after quadloop
                # A[blockmap] += B[...];          generated after quadloop

                cache = ir["piecewise_ir"]["premultiplied_blocks"]

                block_is_transposed = False
                pname = cache.get(unames)

                # Reuse transpose to save memory
                if p["enable_block_transpose_reuse"] and pname is None and len(unames) == 2:
                    pname = cache.get((unames[1], unames[0]))
                    if pname is not None:
                        # Cache hit on transpose
                        block_is_transposed = True

                if pname is None:
                    # Cache miss, precompute block
                    if integral_type == "interior_facet":
                        ptable = multiply_block_interior_facets(0, unames, ttypes, unique_tables,
                                                                unique_table_num_dofs)
                    else:
                        ptable = multiply_block(0, unames, ttypes, unique_tables,
                                                unique_table_num_dofs)
                    pname = "PM%d" % (len(cache, ))
                    cache[unames] = pname
                    unique_tables[pname] = ptable
                    unique_table_types[pname] = "premultiplied"

                block_unames = (pname, )
                blockdata = block_data_t(
                    block_mode, ttypes, fi, factor_is_piecewise, block_unames,
                    block_restrictions, block_is_transposed, block_is_uniform, pname, None, None)
                block_is_piecewise = False

#           elif block_mode == "scaled":
#           # TODO: Add mode, block is piecewise but choose not to be premultiplied
#               # Add to contributions:
#               # FI = sum_q weight * f;          generated inside quadloop
#               # B[...] = FI * u * v;            generated after quadloop
#               # A[blockmap] += B[...];          generated after quadloop
#               raise NotImplementedError("scaled block mode not implemented.")
#               # (probably need mostly the same data as
#               # premultiplied, except no P table name or values)
#               block_is_piecewise = False

            elif block_mode in ("partial", "full", "safe"):
                block_is_piecewise = factor_is_piecewise and not expect_weight
                ma_data = []
                for i, ma in enumerate(ma_indices):
                    if not trs[i].is_piecewise:
                        block_is_piecewise = False
                    ma_data.append(ma_data_t(ma, trs[i]))

                block_is_transposed = False  # FIXME: Handle transposes for these block types

                if block_mode == "partial":
                    # Add to contributions:
                    # P[i] = sum_q weight * f * u[i];  generated inside quadloop
                    # B[i,j] = P[i] * v[j];            generated after quadloop (where v is the piecewise ma)
                    # A[blockmap] += B[...];           generated after quadloop

                    # Find first piecewise index TODO: Is last better? just reverse range here
                    for i in range(rank):
                        if trs[i].is_piecewise:
                            piecewise_ma_index = i
                            break
                    assert rank == 2
                    not_piecewise_ma_index = 1 - piecewise_ma_index
                    block_unames = (unames[not_piecewise_ma_index], )
                    blockdata = block_data_t(block_mode, ttypes, fi,
                                             factor_is_piecewise, block_unames,
                                             block_restrictions, block_is_transposed,
                                             None, None, tuple(ma_data), piecewise_ma_index)
                elif block_mode in ("full", "safe"):
                    # Add to contributions:
                    # B[i] = sum_q weight * f * u[i] * v[j];  generated inside quadloop
                    # A[blockmap] += B[i];                    generated after quadloop

                    block_unames = unames
                    blockdata = block_data_t(block_mode, ttypes, fi,
                                             factor_is_piecewise, block_unames,
                                             block_restrictions, block_is_transposed,
                                             None, None, tuple(ma_data), None)
            else:
                raise RuntimeError("Invalid block_mode %s" % (block_mode, ))

            if block_is_piecewise:
                # Insert in piecewise expr_ir
                ir["piecewise_ir"]["block_contributions"][blockmap].append(blockdata)
            else:
                # Insert in varying expr_ir for this quadrature loop
                block_contributions[blockmap].append(blockdata)

        # Figure out which table names are referenced in unstructured
        # partition
        active_table_names = set()
        for i, v in F.nodes.items():
            tr = v.get('tr')
            if tr is not None and F.nodes[i]['status'] != 'inactive':
                active_table_names.add(tr.name)

        # Figure out which table names are referenced in blocks
        for blockmap, contributions in itertools.chain(
                block_contributions.items(), ir["piecewise_ir"]["block_contributions"].items()):
            for blockdata in contributions:
                if blockdata.block_mode in ("preintegrated", "premultiplied"):
                    active_table_names.add(blockdata.name)
                elif blockdata.block_mode in ("partial", "full", "safe"):
                    for mad in blockdata.ma_data:
                        active_table_names.add(mad.tabledata.name)

        # Record all table types before dropping tables
        ir["unique_table_types"].update(unique_table_types)

        # Drop tables not referenced from modified terminals
        # and tables of zeros and ones
        unused_ttypes = ("zeros", "ones", "quadrature")
        keep_table_names = set()
        for name in active_table_names:
            ttype = ir["unique_table_types"][name]
            if ttype not in unused_ttypes:
                if name in unique_tables:
                    keep_table_names.add(name)
        unique_tables = {name: unique_tables[name] for name in keep_table_names}

        # Add to global set of all tables
        for name, table in unique_tables.items():
            tbl = ir["unique_tables"].get(name)
            if tbl is not None and not numpy.allclose(
                    tbl, table, rtol=p["table_rtol"], atol=p["table_atol"]):
                raise RuntimeError("Table values mismatch with same name.")
        ir["unique_tables"].update(unique_tables)

        # Analyse active terminals to check what we'll need to generate code for
        active_mts = []
        for i, v in F.nodes.items():
            mt = v.get('mt', False)
            if mt and F.nodes[i]['status'] != 'inactive':
                active_mts.append(mt)

        # Figure out if we need to access CellCoordinate to avoid
        # generating quadrature point table otherwise
        if integral_type == "cell":
            need_points = any(isinstance(mt.terminal, CellCoordinate) for mt in active_mts)
        elif integral_type in facet_integral_types:
            need_points = any(isinstance(mt.terminal, FacetCoordinate) for mt in active_mts)
        elif integral_type in custom_integral_types:
            need_points = True  # TODO: Always?
        else:
            need_points = False

        # Figure out if we need to access QuadratureWeight to avoid
        # generating quadrature point table otherwise need_weights =
        # any(isinstance(mt.terminal, QuadratureWeight) for mt in
        # active_mts)

        # Count blocks of each mode
        block_modes = collections.defaultdict(int)
        for blockmap, contributions in block_contributions.items():
            for blockdata in contributions:
                block_modes[blockdata.block_mode] += 1

        # Debug output
        summary = "\n".join(
            "  {}\t{}".format(count, mode) for mode, count in sorted(block_modes.items()))
        logger.debug("Blocks of each mode: {}".format(summary))

        # If there are any blocks other than preintegrated we need weights
        if expect_weight and any(mode != "preintegrated" for mode in block_modes):
            need_weights = True
        elif integral_type in custom_integral_types:
            need_weights = True  # TODO: Always?
        else:
            need_weights = False

        # Build IR dict for the given expressions
        # Store final ir for this num_points
        ir["varying_irs"][num_points] = {"factorization": F,
                                         "modified_arguments": [F.nodes[i]['mt'] for i in argkeys],
                                         "block_contributions": block_contributions,
                                         "need_points": need_points,
                                         "need_weights": need_weights}
    return ir


def analyse_dependencies(F, mt_unique_table_reference):
    # Sets 'status' of all nodes to either: 'inactive', 'piecewise' or 'varying'
    # Children of 'target' nodes are either 'piecewise' or 'varying'.
    # All other nodes are 'inactive'.
    # Varying nodes are identified by their tables ('tr'). All their parent
    # nodes are also set to 'varying' - any remaining active nodes are 'piecewise'.

    # Set targets, and dependencies to 'active'
    targets = [i for i, v in F.nodes.items() if v.get('target')]
    for i, v in F.nodes.items():
        v['status'] = 'inactive'

    while targets:
        s = targets.pop()
        F.nodes[s]['status'] = 'active'
        for j in F.out_edges[s]:
            if F.nodes[j]['status'] == 'inactive':
                targets.append(j)

    # Build piecewise/varying markers for factorized_vertices
    varying_ttypes = ("varying", "uniform", "quadrature")
    varying_indices = []
    for i, v in F.nodes.items():
        if v.get('mt') is None:
            continue
        tr = v.get('tr')
        if tr is not None:
            ttype = tr.ttype
            # Check if table computations have revealed values varying over points
            # Note: uniform means entity-wise uniform, varying over points
            if ttype in varying_ttypes:
                varying_indices.append(i)
            else:
                if ttype not in ("fixed", "piecewise", "ones", "zeros"):
                    raise RuntimeError("Invalid ttype %s" % (ttype, ))

        elif not is_cellwise_constant(v['expression']):
            raise RuntimeError("Error")
            # Keeping this check to be on the safe side,
            # not sure which cases this will cover (if any)
            # varying_indices.append(i)

    # Set all parents of active varying nodes to 'varying'
    while varying_indices:
        s = varying_indices.pop()
        if F.nodes[s]['status'] == 'active':
            F.nodes[s]['status'] = 'varying'
            for j in F.in_edges[s]:
                varying_indices.append(j)

    # Any remaining active nodes must be 'piecewise'
    for i, v in F.nodes.items():
        if v['status'] == 'active':
            v['status'] = 'piecewise'


def replace_quadratureweight(expression):
    """Remove any QuadratureWeight terminals and replace with 1.0."""

    r = _find_terminals_in_ufl_expression(expression, QuadratureWeight)
    replace_map = {q: 1.0 for q in r}

    return ufl.algorithms.replace(expression, replace_map)


def _find_terminals_in_ufl_expression(e, etype):
    """Recursively search expression for terminals of type etype."""
    r = []
    for op in e.ufl_operands:
        if is_modified_terminal(op) and isinstance(op, etype):
            r.append(op)
        else:
            r += _find_terminals_in_ufl_expression(op, etype)

    return r
