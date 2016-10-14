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

from ufl import product
from ufl.checks import is_cellwise_constant
from ufl.classes import CellCoordinate, FacetCoordinate

from ffc.uflacs.analysis.balancing import balance_modifiers
from ffc.uflacs.analysis.modified_terminals import is_modified_terminal, analyse_modified_terminal
from ffc.uflacs.analysis.graph import build_graph
from ffc.uflacs.analysis.graph_vertices import build_scalar_graph_vertices
from ffc.uflacs.analysis.graph_rebuild import rebuild_with_scalar_subexpressions
from ffc.uflacs.analysis.graph_dependencies import compute_dependencies, mark_active, mark_image
from ffc.uflacs.analysis.graph_ssa import compute_dependency_count, invert_dependencies
#from ffc.uflacs.analysis.graph_ssa import default_cache_score_policy, compute_cache_scores, allocate_registers
from ffc.uflacs.analysis.factorization import compute_argument_factorization
from ffc.uflacs.elementtables.terminaltables import build_optimized_tables


def build_uflacs_ir(cell, integral_type, entitytype,
                    integrands, coefficient_numbering,
                    quadrature_rules, parameters):
    uflacs_ir = {}

    # { ufl coefficient: count }
    uflacs_ir["coefficient_numbering"] = coefficient_numbering

    # { num_points: expr_ir for one integrand }
    uflacs_ir["expr_irs"] = {}

    # Build the core uflacs expression ir for each num_points/integrand
    # TODO: Better to compute joint IR for all integrands
    #       and deal with num_points later?
    #       I.e. common_expr_ir = compute_common_expr_ir(integrands)
    #       If we want to adjoint quadrature rules for subterms
    #       automatically anyway, num_points should be advisory.
    #       For now, expecting multiple num_points to be rare.
    for num_points in sorted(integrands.keys()):
        expr_ir = compute_expr_ir(integrands[num_points])

        uflacs_ir["expr_irs"][num_points] = expr_ir

        # Build set of modified terminal ufl expressions
        V = expr_ir["V"]
        modified_terminals = [analyse_modified_terminal(V[i])
                              for i in expr_ir["modified_terminal_indices"]]
        terminal_data = modified_terminals + expr_ir["modified_arguments"]

        # FIXME: For custom integrals, skip table building but set up
        # the necessary table names and classname mappings instead

        # FIXME: Want table information earlier, even before scalar
        # rebuilding! Must split compute_expr_ir to achieve this.
        # FIXME: Store table type as fourth entry in table ranges
        unique_tables, mt_table_ranges, table_types = \
            build_optimized_tables(num_points, quadrature_rules,
                cell, integral_type, entitytype, terminal_data, parameters)

        # Figure out if we need to access CellCoordinate to
        # avoid generating quadrature point table otherwise
        if integral_type == "cell":
            expr_ir["need_points"] = any(isinstance(mt.terminal, CellCoordinate)
                                         for mt in modified_terminals)
        elif integral_type in ("interior_facet", "exterior_facet"):
            expr_ir["need_points"] = any(isinstance(mt.terminal, FacetCoordinate)
                                         for mt in modified_terminals)
        else:
            expr_ir["need_points"] = False


        # Ordered table data
        terminal_table_ranges = [mt_table_ranges.get(mt) for mt in terminal_data]

        # Split into arguments and other terminals before storing in expr_ir
        # TODO: Some tables are associated with num_points, some are not
        #       (i.e. piecewise constant, averaged and x0).
        #       It will be easier to deal with that if we can join
        #       the expr_ir for all num_points as mentioned above.
        n = len(expr_ir["modified_terminal_indices"])
        m = len(expr_ir["modified_arguments"])
        assert len(terminal_data) == n + m
        assert len(terminal_table_ranges) == n + m
        expr_ir["modified_terminal_table_ranges"] = terminal_table_ranges[:n]
        expr_ir["modified_argument_table_ranges"] = terminal_table_ranges[n:]

        # Store table data in V indexing, this is used in integralgenerator
        expr_ir["table_ranges"] = numpy.empty(len(V), dtype=object)
        expr_ir["table_ranges"][expr_ir["modified_terminal_indices"]] = \
            expr_ir["modified_terminal_table_ranges"]

        # FIXME: Drop tables for Real and DG0 elements (all 1.0 for each dof)

        # FIXME: Replace coefficients with empty dofrange with zero (which are these?)
        # FIXME: Propagate constants

        # Drop factorization terms where table dof range is
        # empty for any of the modified arguments
        AF = expr_ir["argument_factorization"]
        MATR = expr_ir["modified_argument_table_ranges"]
        for mas in list(AF.keys()):
            for j in mas:
                dofrange = MATR[j][1:3]
                if dofrange[0] == dofrange[1]:
                    del AF[mas]
                    break
        # FIXME: Propagate dependencies back to remove expressions
        # not used anymore after dropping factorization terms

        # Drop tables not referenced from modified terminals
        # and and tables of zeros and ones
        used_table_names = set()
        for tabledata in terminal_table_ranges:
            if tabledata is not None:
                name, begin, end = tabledata
                if table_types[name] not in ("zeros", "ones"):
                    used_table_names.add(name)
        if None in used_table_names:
            used_table_names.remove(None)
        unique_tables = { name: unique_tables[name] for name in used_table_names }

        # Store the tables and ranges
        expr_ir["table_types"] = table_types
        expr_ir["unique_tables"] = unique_tables

    return uflacs_ir


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

    # Build new list representation of graph where all vertices
    # of V represent single scalar operations
    e2i, V, target_variables = build_scalar_graph_vertices(scalar_expressions)

    return e2i, V, target_variables


def compute_argument_factorization2(expressions):
    # TODO: Can we merge these three calls to something more efficient overall?

    # Build scalar list-based graph representation
    e2i, V, target_variables = build_scalar_graph(expressions)

    # Compute sparse dependency matrix
    dependencies = compute_dependencies(e2i, V)

    # Compute factorization of arguments
    argument_factorization, modified_arguments, V, target_variables, dependencies = \
        compute_argument_factorization(V, target_variables, dependencies)
    return argument_factorization, modified_arguments, V, target_variables, dependencies


def analyse_dependencies(V, target_variables, modified_terminal_indices, dependencies):
    # Count the number of dependencies every subexpr has
    depcount = compute_dependency_count(dependencies)

    # Build the 'inverse' of the sparse dependency matrix
    inverse_dependencies = invert_dependencies(dependencies, depcount)

    # Mark subexpressions of V that are actually needed for final result
    active, num_active = mark_active(dependencies, target_variables)

    # Build piecewise/varying markers for factorized_vertices
    # FIXME: have better measure for spatial dependency now
    spatially_dependent_terminal_indices = [i for i in modified_terminal_indices
                                            if not is_cellwise_constant(V[i])]
    varying, num_spatial = mark_image(inverse_dependencies,
                                      spatially_dependent_terminal_indices)
    piecewise = 1 - varying
    # Skip non-active things
    varying *= active
    piecewise *= active

    # TODO: Skip literals in both varying and piecewise
    # nonliteral = ...
    # varying *= nonliteral
    # piecewise *= nonliteral

    return dependencies, inverse_dependencies, active, piecewise, varying


def compute_expr_ir(expressions):
    """Build the intermediate representation for a list of expressions."""
    
    # Wrap in list if we only get one expression
    if not isinstance(expressions, list):
        expressions = [expressions]

    # TODO: Apply this transformation before calling compute_expr_ir?
    expressions = [balance_modifiers(expr) for expr in expressions]

    argument_factorization, modified_arguments, V, target_variables, dependencies = \
        compute_argument_factorization2(expressions)


    # Store modified arguments in analysed form
    for i in range(len(modified_arguments)):
        modified_arguments[i] = analyse_modified_terminal(modified_arguments[i])

    # Build set of modified_terminal indices into factorized_vertices
    modified_terminal_indices = [i for i, v in enumerate(V)
                                 if is_modified_terminal(v)]


    # Build IR for the given expressions
    expr_ir = {}

    ### Core expression graph:
    # (array) V-index -> UFL subexpression
    expr_ir["V"] = V

    # (array) Flattened input expression component index -> V-index
    expr_ir["target_variables"] = target_variables

    ### Result of factorization:
    # (array) MA-index -> UFL expression of modified arguments
    expr_ir["modified_arguments"] = modified_arguments

    # (dict) tuple(MA-indices) -> V-index of monomial factor
    expr_ir["argument_factorization"] = argument_factorization

    ### Modified terminals
    # (array) list of V-indices to modified terminals
    expr_ir["modified_terminal_indices"] = modified_terminal_indices


    # FIXME: Split function here! Need to get and analyze table data before the below dependency analysis.


    # --- Various dependency analysis ---
    dependencies, inverse_dependencies, active, piecewise, varying = \
        analyse_dependencies(V, target_variables, modified_terminal_indices, dependencies)

    # Dependency structure of graph:
    #expr_ir["dependencies"] = dependencies                           # (CRSArray) V-index -> direct dependency V-index list
    #expr_ir["inverse_dependencies"] = inverse_dependencies           # (CRSArray) V-index -> direct dependee V-index list

    # Metadata about each vertex
    #expr_ir["active"] = active       # (array) V-index -> bool
    expr_ir["piecewise"] = piecewise  # (array) V-index -> bool
    expr_ir["varying"] = varying     # (array) V-index -> bool

    return expr_ir



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

    # rank = max(len(k) for k in argument_factorization.keys())
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

