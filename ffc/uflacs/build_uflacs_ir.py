# -*- coding: utf-8 -*-
# Copyright (C) 2011-2016 Martin Sandve Alnæs
#
# This file is part of UFLACS.
#
# UFLACS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFLACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>.

"""Main algorithm for building the uflacs intermediate representation."""

import numpy

from ufl import product, as_ufl
from ufl.log import error
from ufl.checks import is_cellwise_constant
from ufl.classes import CellCoordinate, FacetCoordinate, QuadratureWeight

from ffc.uflacs.analysis.balancing import balance_modifiers
from ffc.uflacs.analysis.modified_terminals import is_modified_terminal, analyse_modified_terminal
from ffc.uflacs.analysis.graph import build_graph
from ffc.uflacs.analysis.graph_vertices import build_scalar_graph_vertices
from ffc.uflacs.analysis.graph_rebuild import rebuild_with_scalar_subexpressions
from ffc.uflacs.analysis.graph_dependencies import compute_dependencies, mark_active, mark_image
from ffc.uflacs.analysis.graph_ssa import compute_dependency_count, invert_dependencies
#from ffc.uflacs.analysis.graph_ssa import default_cache_score_policy, compute_cache_scores, allocate_registers
from ffc.uflacs.analysis.factorization import compute_argument_factorization
from ffc.uflacs.elementtables import build_optimized_tables


def build_uflacs_ir(cell, integral_type, entitytype,
                    integrands, tensor_shape,
                    coefficient_numbering,
                    quadrature_rules, parameters):
    ir = {}

    # { ufl coefficient: count }
    ir["coefficient_numbering"] = coefficient_numbering

    rank = len(tensor_shape)
    
    # { num_points: expr_ir for one integrand }
    ir["expr_irs"] = {}

    # Build the core uflacs expression ir for each num_points/integrand
    # TODO: Better to compute joint IR for all integrands
    #       and deal with num_points later?
    #       I.e. common_expr_ir = compute_common_expr_ir(integrands)
    #       If we want to adjoint quadrature rules for subterms
    #       automatically anyway, num_points should be advisory.
    #       For now, expecting multiple num_points to be rare.
    for num_points in sorted(integrands.keys()):
        expressions = [integrands[num_points]]

        # TODO: Apply this transformation to integrands earlier?
        expressions = [balance_modifiers(expr) for expr in expressions]

        # Build scalar list-based graph representation
        V, V_deps, V_targets = build_scalar_graph(expressions)


        # Build terminal_data from V here before factorization.
        # Then we can use it to derive table properties for all modified terminals,
        # and then use that to rebuild the scalar graph more efficiently before
        # argument factorization. We can build terminal_data again after factorization
        # if that's necessary.
        initial_terminal_indices = [i for i, v in enumerate(V)
                                    if is_modified_terminal(v)]
        initial_terminal_data = [analyse_modified_terminal(V[i])
                                 for i in initial_terminal_indices]
        unique_tables, mt_table_ranges, table_types = \
            build_optimized_tables(num_points, quadrature_rules,
                cell, integral_type, entitytype, initial_terminal_data, parameters)

        # Build replacement map for modified terminals with zero tables
        z = as_ufl(0.0)
        for i, mt in zip(initial_terminal_indices, initial_terminal_data):
            tr = mt_table_ranges.get(mt)
            if tr is not None:
                uname, begin, end = tr
                ttype = table_types[uname]
                # Any modified terminal with zero table is itself a zero value
                if ttype == "zeros":
                    V[i] = z
        # Propagate expression changes
        # (could possibly use replace() on target expressions instead)
        for i in range(len(V)):
            deps = [V[j] for j in V_deps[i]]
            if deps:
                V[i] = V[i]._ufl_expr_reconstruct_(*deps)

        # Rebuild scalar target expressions and graph
        # (this may be overkill and possible to optimize
        # away if it turns out to be costly)
        expressions = [V[i] for i in V_targets]

        # Rebuild scalar list-based graph representation
        SV, SV_deps, SV_targets = build_scalar_graph(expressions)
        assert all(i < len(SV) for i in SV_targets)


        # Compute factorization of arguments
        (argument_factorizations, modified_arguments,
             FV, FV_deps, FV_targets) = \
            compute_argument_factorization(SV, SV_deps, SV_targets, rank)
        assert len(SV_targets) == len(argument_factorizations)

        # TODO: Still expecting one target variable in code generation
        assert len(argument_factorizations) == 1
        argument_factorization, = argument_factorizations

        # Store modified arguments in analysed form
        for i in range(len(modified_arguments)):
            modified_arguments[i] = analyse_modified_terminal(modified_arguments[i])

        # Build set of modified_terminal indices into factorized_vertices
        modified_terminal_indices = [i for i, v in enumerate(FV)
                                     if is_modified_terminal(v)]

        # Build set of modified terminal ufl expressions
        modified_terminals = [analyse_modified_terminal(FV[i])
                              for i in modified_terminal_indices]

        # Organize table data more, split into arguments and other terminals
        modified_terminal_table_ranges = [mt_table_ranges.get(mt)
                                          for mt in modified_terminals]
        modified_argument_table_ranges = [mt_table_ranges.get(mt)
                                          for mt in modified_arguments]


        # Dependency analysis
        inv_FV_deps, FV_active, FV_piecewise, FV_varying = \
            analyse_dependencies(FV, FV_deps, FV_targets,
                                 modified_terminal_indices,
                                 mt_table_ranges,
                                 table_types)

        # Mark active modified arguments
        #active_modified_arguments = numpy.zeros(len(modified_arguments), dtype=int)
        #for ma_indices in argument_factorization:
        #    for j in ma_indices:
        #        active_modified_arguments[j] = 1


        # Figure out which table names are active
        active_table_names = set()
        for i, tr in zip(modified_terminal_indices, modified_terminal_table_ranges):
            if FV_active[i] and tr is not None:
                active_table_names.add(tr[0])
        for ma_indices in argument_factorization:
            for j in ma_indices:
                tr = modified_argument_table_ranges[j]
                if tr is not None:
                    active_table_names.add(tr[0])

        # Drop tables not referenced from modified terminals
        # and and tables of zeros and ones
        unused_types = ("zeros", "ones", "quadrature")
        used_table_names = set(name for name in active_table_names
                               if name is not None
                                  and table_types[name] not in unused_types)
        unique_tables = { name: unique_tables[name] for name in used_table_names }


        # Analyse active terminals to check what we'll need to generate code for
        active_mts = [mt for i, mt in zip(modified_terminal_indices, modified_terminals)
                      if FV_active[i]]

        # Figure out if we need to access CellCoordinate to
        # avoid generating quadrature point table otherwise
        if integral_type == "cell":
            need_points = any(isinstance(mt.terminal, CellCoordinate)
                              for mt in active_mts)
        elif integral_type in ("interior_facet", "exterior_facet"):
            need_points = any(isinstance(mt.terminal, FacetCoordinate)
                              for mt in active_mts)
        else:
            need_points = False

        # Figure out if we need to access QuadratureWeight to
        # avoid generating quadrature point table otherwise
        need_weights = any(isinstance(mt.terminal, QuadratureWeight)
                           for mt in active_mts)

        # Loop over factorization terms
        from collections import defaultdict
        block_contributions = {
            # TODO: Should not store piecewise blocks inside num_points context
            "piecewise": defaultdict(list),
            "varying": defaultdict(list)
            }
        for ma_indices, fi in sorted(argument_factorization.items()):
            # Get a bunch of information about this term
            rank = len(ma_indices)
            trs = tuple(modified_argument_table_ranges[ai] for ai in ma_indices)
            unames = tuple(tr[0] for tr in trs)
            dofblock = tuple(tr[1:3] for tr in trs)
            ttypes = tuple(table_types[name] for name in unames)
            assert not any(tt == "zeros" for tt in ttypes)

            piecewise_types = ("piecewise", "fixed", "ones")
            if FV_piecewise[fi] and all(tt in piecewise_types for tt in ttypes):
                contributions = block_contributions["piecewise"][dofblock]
            else:
                contributions = block_contributions["varying"][dofblock]

            data = (ma_indices, fi, trs, unames, ttypes)
            contributions.append(data)


        # Build IR dict for the given expressions
        expr_ir = {}

        expr_ir["block_contributions"] = block_contributions
        
        # (array) FV-index -> UFL subexpression
        expr_ir["V"] = FV

        # (array) Flattened input expression component index -> FV-index
        expr_ir["target_variables"] = FV_targets

        ### Result of factorization:
        # (array) MA-index -> UFL expression of modified arguments
        expr_ir["modified_arguments"] = modified_arguments

        # (dict) tuple(MA-indices) -> FV-index of monomial factor
        expr_ir["argument_factorization"] = argument_factorization

        ### Modified terminals
        # (array) list of FV-indices to modified terminals
        expr_ir["modified_terminal_indices"] = modified_terminal_indices

        # Dependency structure of graph:
        # (CRSArray) FV-index -> direct dependency FV-index list
        #expr_ir["dependencies"] = FV_deps

        # (CRSArray) FV-index -> direct dependee FV-index list
        #expr_ir["inverse_dependencies"] = inv_FV_deps

        # Metadata about each vertex
        expr_ir["active"] = FV_active        # (array) FV-index -> bool
        expr_ir["piecewise"] = FV_piecewise  # (array) FV-index -> bool
        expr_ir["varying"] = FV_varying      # (array) FV-index -> bool

        expr_ir["modified_terminal_table_ranges"] = modified_terminal_table_ranges
        expr_ir["modified_argument_table_ranges"] = modified_argument_table_ranges

        # Store table data in FV indexing, this is used in integralgenerator
        expr_ir["table_ranges"] = numpy.empty(len(FV), dtype=object)
        expr_ir["table_ranges"][expr_ir["modified_terminal_indices"]] = \
            expr_ir["modified_terminal_table_ranges"]

        expr_ir["need_points"] = need_points
        expr_ir["need_weights"] = need_weights

        # Store the tables and ranges
        expr_ir["table_types"] = table_types
        expr_ir["unique_tables"] = unique_tables


        # TODO: Some tables are associated with num_points, some are not
        #       (i.e. piecewise constant, averaged and x0).
        #       It will be easier to deal with that if we can join
        #       the expr_ir for all num_points as mentioned above.
        ir["expr_irs"][num_points] = expr_ir

    return ir


def build_scalar_graph(expressions):
    """Build list representation of expression graph covering the given expressions.

    TODO: Renaming, refactoring and cleanup of the graph building algorithms used in here
    """

    # Build the initial coarse computational graph of the expression
    G = build_graph(expressions)

    assert len(expressions) == 1, "FIXME: Multiple expressions in graph building needs more work from this point on."

    # Build more fine grained computational graph of scalar subexpressions
    # TODO: Make it so that
    #   expressions[k] <-> NV[nvs[k][:]],
    #   len(nvs[k]) == value_size(expressions[k])
    scalar_expressions = rebuild_with_scalar_subexpressions(G)

    # Sanity check on number of scalar symbols/components
    assert len(scalar_expressions) == sum(product(expr.ufl_shape) for expr in expressions)

    # Build new list representation of graph where all
    # vertices of V represent single scalar operations
    e2i, V, V_targets = build_scalar_graph_vertices(scalar_expressions)

    # Compute sparse dependency matrix
    V_deps = compute_dependencies(e2i, V)

    return V, V_deps, V_targets


def analyse_dependencies(V, V_deps, V_targets,
                         modified_terminal_indices,
                         mt_table_ranges,
                         table_types):
    # Count the number of dependencies every subexpr has
    V_depcount = compute_dependency_count(V_deps)

    # Build the 'inverse' of the sparse dependency matrix
    inv_deps = invert_dependencies(V_deps, V_depcount)

    # Mark subexpressions of V that are actually needed for final result
    active, num_active = mark_active(V_deps, V_targets)

    # Build piecewise/varying markers for factorized_vertices
    varying_indices = []
    for i in modified_terminal_indices:

        # TODO: Can probably avoid this re-analysis by
        # passing other datastructures in here:
        mt = analyse_modified_terminal(V[i])
        tr = mt_table_ranges.get(mt)
        if tr is not None:
            # Check if table computations have revealed values varying over points
            uname = tr[0]
            ttype = table_types[uname]
            # Note: uniform means entity-wise uniform, varying over points
            if ttype in ("varying", "uniform", "quadrature"):
                varying_indices.append(i)
            else:
                if ttype not in ("fixed", "piecewise", "ones", "zeros"):
                    error("Invalid ttype %s" % (ttype,))

        elif not is_cellwise_constant(V[i]):
            # Keeping this check to be on the safe side,
            # not sure which cases this will cover (if any)
            varying_indices.append(i)

    # Mark every subexpression that is computed
    # from the spatially dependent terminals
    varying, num_varying = mark_image(inv_deps, varying_indices)

    # The rest of the subexpressions are piecewise constant (1-1=0, 1-0=1)
    piecewise = 1 - varying

    # Unmark non-active subexpressions
    varying *= active
    piecewise *= active

    # TODO: Skip literals in both varying and piecewise
    # nonliteral = ...
    # varying *= nonliteral
    # piecewise *= nonliteral

    return inv_deps, active, piecewise, varying


# TODO: Consider comments below and do it or delete them.

""" Old comments:

Work for later::

        - Apply some suitable renumbering of vertices and corresponding arrays prior to returning

        - Allocate separate registers for each partition
          (but e.g. argument[iq][i0] may need to be accessible in other loops)

        - Improve register allocation algorithm

        - Take a list of expressions as input to compile several expressions in one joined graph
          (e.g. to compile a,L,M together for nonlinear problems)

"""


""" # Old comments:

    # TODO: Inspection of varying shows that factorization is
    # needed for effective loop invariant code motion w.r.t. quadrature loop as well.
    # Postphoning that until everything is working fine again.
    # Core ingredients for such factorization would be:
    # - Flatten products of products somehow
    # - Sorting flattened product factors by loop dependency then by canonical ordering
    # Or to keep binary products:
    # - Rebalancing product trees ((a*c)*(b*d) -> (a*b)*(c*d)) to make piecewise quantities 'float' to the top of the list

    # rank = max(len(ma_indices) for ma_indices in argument_factorization)
    # for i,a in enumerate(modified_arguments):
    #    iarg = a.number()
    # ipart = a.part()

    # TODO: More structured MA organization?
    #modified_arguments[rank][block][entry] -> UFL expression of modified argument
    #dofranges[rank][block] -> (begin, end)
    # or
    #modified_arguments[rank][entry] -> UFL expression of modified argument
    #dofrange[rank][entry] -> (begin, end)
    #argument_factorization: (dict) tuple(MA-indices (only relevant ones!)) -> V-index of monomial factor
    # becomes
    #argument_factorization: (dict) tuple(entry for each(!) rank) -> V-index of monomial factor ## doesn't cover intermediate f*u in f*u*v!
"""


"""
def old_code_useful_for_optimization():

    # Use heuristics to mark the usefulness of storing every subexpr in a variable
    scores = compute_cache_scores(V,
                                  active,
                                  dependencies,
                                  inverse_dependencies,
                                  partitions,  # TODO: Rewrite in terms of something else, this doesn't exist anymore
                                  cache_score_policy=default_cache_score_policy)

    # Allocate variables to store subexpressions in
    allocations = allocate_registers(active, partitions, target_variables,
                                     scores, int(parameters["max_registers"]), int(parameters["score_threshold"]))
    target_registers = [allocations[r] for r in target_variables]
    num_registers = sum(1 if x >= 0 else 0 for x in allocations)
    # TODO: If we renumber we can allocate registers separately for each partition, which is probably a good idea.

    expr_oir = {}
    expr_oir["num_registers"] = num_registers
    expr_oir["partitions"] = partitions
    expr_oir["allocations"] = allocations
    expr_oir["target_registers"] = target_registers
    return expr_oir
"""

