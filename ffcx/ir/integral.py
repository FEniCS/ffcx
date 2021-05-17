# Copyright (C) 2013-2020 Martin Sandve Alnæs and Michal Habera
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Main algorithm for building the integral intermediate representation."""

import collections
import itertools
import logging

import numpy

import ufl
from ffcx.ir.analysis.factorization import \
    compute_argument_factorization
from ffcx.ir.analysis.graph import build_scalar_graph
from ffcx.ir.analysis.modified_terminals import (
    analyse_modified_terminal, is_modified_terminal)
from ffcx.ir.analysis.visualise import visualise_graph
from ffcx.ir.elementtables import build_optimized_tables
from ufl.algorithms.balancing import balance_modifiers
from ufl.checks import is_cellwise_constant
from ufl.classes import QuadratureWeight

logger = logging.getLogger("ffcx")

ma_data_t = collections.namedtuple("ma_data_t", ["ma_index", "tabledata"])

block_data_t = collections.namedtuple("block_data_t",
                                      ["ttypes",  # list of table types for each block rank
                                       "factor_indices_comp_indices",  # list of tuples (factor index, component index)
                                       "all_factors_piecewise",  # True if all factors for this block are piecewise
                                       "unames",  # list of unique FE table names for each block rank
                                       "restrictions",  # restriction "+" | "-" | None for each block rank
                                       "transposed",  # block is the transpose of another
                                       "is_uniform",
                                       "name",  # used in "preintegrated" and "premultiplied"
                                       "ma_data",  # used in "full", "safe" and "partial"
                                       "piecewise_ma_index",  # used in "partial"
                                       "is_permuted"  # Do quad points on facets need to be permuted?
                                       ])


def compute_integral_ir(cell, integral_type, entitytype, integrands, argument_shape,
                        p, visualise):
    # The intermediate representation dict we're building and returning
    # here
    ir = {}

    # Pass on parameters for consumption in code generation
    ir["params"] = p

    # Shared unique tables for all quadrature loops
    ir["unique_tables"] = {}
    ir["unique_table_types"] = {}

    ir["integrand"] = {}

    ir["table_dofmaps"] = {}
    ir["table_dof_base_transformations"] = {}
    ir["needs_transformation_data"] = 0
    ir["table_needs_transformation_data"] = {}

    for quadrature_rule, integrand in integrands.items():

        expression = integrand

        # Rebalance order of nested terminal modifiers
        expression = balance_modifiers(expression)

        # Remove QuadratureWeight terminals from expression and replace with 1.0
        expression = replace_quadratureweight(expression)

        # Build initial scalar list-based graph representation
        S = build_scalar_graph(expression)

        # Build terminal_data from V here before factorization. Then we
        # can use it to derive table properties for all modified
        # terminals, and then use that to rebuild the scalar graph more
        # efficiently before argument factorization. We can build
        # terminal_data again after factorization if that's necessary.

        initial_terminals = {i: analyse_modified_terminal(v['expression'])
                             for i, v in S.nodes.items()
                             if is_modified_terminal(v['expression'])}

        mt_unique_table_reference = build_optimized_tables(
            quadrature_rule,
            cell,
            integral_type,
            entitytype,
            initial_terminals.values(),
            ir["unique_tables"],
            rtol=p["table_rtol"],
            atol=p["table_atol"])
        unique_tables = {v.name: v.values for v in mt_unique_table_reference.values()}
        unique_table_types = {v.name: v.ttype for v in mt_unique_table_reference.values()}

        for td in mt_unique_table_reference.values():
            ir["table_needs_transformation_data"][td.name] = td.needs_transformation_data
            if td.needs_transformation_data:
                ir["needs_transformation_data"] = 1
            ir["table_dof_base_transformations"][td.name] = td.dof_base_transformations
            ir["table_dofmaps"][td.name] = td.dofmap

        S_targets = [i for i, v in S.nodes.items() if v.get('target', False)]

        if 'zeros' in unique_table_types.values() and len(S_targets) == 1:
            # If there are any 'zero' tables, replace symbolically and rebuild graph
            #
            # TODO: Implement zero table elimination for non-scalar graphs
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
            expression = S.nodes[S_targets[0]]['expression']

            # Rebuild scalar list-based graph representation
            S = build_scalar_graph(expression)

        # Output diagnostic graph as pdf
        if visualise:
            visualise_graph(S, 'S.pdf')

        # Compute factorization of arguments
        rank = len(argument_shape)
        F = compute_argument_factorization(S, rank)

        # Get the 'target' nodes that are factors of arguments, and insert in dict
        FV_targets = [i for i, v in F.nodes.items() if v.get('target', False)]
        argument_factorization = {}

        for fi in FV_targets:
            # Number of blocks using this factor must agree with number of components
            # to which this factor contributes. I.e. there are more blocks iff there are more
            # components
            assert len(F.nodes[fi]['target']) == len(F.nodes[fi]['component'])

            k = 0
            for w in F.nodes[fi]['target']:
                comp = F.nodes[fi]['component'][k]
                argument_factorization[w] = argument_factorization.get(w, [])

                # Store tuple of (factor index, component index)
                argument_factorization[w].append((fi, comp))
                k += 1

        # Get list of indices in F which are the arguments (should be at start)
        argkeys = set()
        for w in argument_factorization:
            argkeys = argkeys | set(w)
        argkeys = list(argkeys)

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

        # Output diagnostic graph as pdf
        if visualise:
            visualise_graph(F, 'F.pdf')

        # Loop over factorization terms
        block_contributions = collections.defaultdict(list)
        for ma_indices, fi_ci in sorted(argument_factorization.items()):
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

            # Check if each *each* factor corresponding to this argument is piecewise
            all_factors_piecewise = all(F.nodes[ifi[0]]["status"] == 'piecewise' for ifi in fi_ci)
            block_is_permuted = False
            for n in unames:
                if unique_tables[n].shape[0] > 1:
                    block_is_permuted = True
            ma_data = []
            for i, ma in enumerate(ma_indices):
                ma_data.append(ma_data_t(ma, trs[i]))

            block_is_transposed = False  # FIXME: Handle transposes for these block types

            block_unames = unames
            blockdata = block_data_t(ttypes, fi_ci,
                                     all_factors_piecewise, block_unames,
                                     block_restrictions, block_is_transposed,
                                     block_is_uniform, None, tuple(ma_data), None, block_is_permuted)

            # Insert in expr_ir for this quadrature loop
            block_contributions[blockmap].append(blockdata)

        # Figure out which table names are referenced
        active_table_names = set()
        for i, v in F.nodes.items():
            tr = v.get('tr')
            if tr is not None and F.nodes[i]['status'] != 'inactive':
                active_table_names.add(tr.name)

        # Figure out which table names are referenced in blocks
        for blockmap, contributions in itertools.chain(
                block_contributions.items()):
            for blockdata in contributions:
                for mad in blockdata.ma_data:
                    active_table_names.add(mad.tabledata.name)

        # Record all table types before dropping tables
        ir["unique_table_types"].update(unique_table_types)

        # Drop tables not referenced from modified terminals
        # and tables of zeros and ones
        unused_ttypes = ("zeros", "ones")
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
        # active_mts = []
        # for i, v in F.nodes.items():
        #     mt = v.get('mt', False)
        #     if mt and F.nodes[i]['status'] != 'inactive':
        #         active_mts.append(mt)

        # Build IR dict for the given expressions
        # Store final ir for this num_points
        ir["integrand"][quadrature_rule] = {"factorization": F,
                                            "modified_arguments": [F.nodes[i]['mt'] for i in argkeys],
                                            "block_contributions": block_contributions}
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
    varying_ttypes = ("varying", "quadrature", "uniform")
    varying_indices = []
    for i, v in F.nodes.items():
        if v.get('mt') is None:
            continue
        tr = v.get('tr')
        if tr is not None:
            ttype = tr.ttype
            # Check if table computations have revealed values varying over points
            if ttype in varying_ttypes:
                varying_indices.append(i)
            else:
                if ttype not in ("fixed", "piecewise", "ones", "zeros"):
                    raise RuntimeError("Invalid ttype %s" % (ttype, ))

        elif not is_cellwise_constant(v['expression']):
            raise RuntimeError("Error " + str(tr))
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

    r = []
    for node in ufl.corealg.traversal.unique_pre_traversal(expression):
        if is_modified_terminal(node) and isinstance(node, QuadratureWeight):
            r.append(node)

    replace_map = {q: 1.0 for q in r}
    return ufl.algorithms.replace(expression, replace_map)
