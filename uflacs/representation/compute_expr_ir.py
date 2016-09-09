# -*- coding: utf-8 -*-
# Copyright (C) 2011-2015 Martin Sandve Alnæs
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
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>

"""Algorithms for the representation phase of the form compilation."""

from ufl import product
from ufl.checks import is_cellwise_constant
from uflacs.analysis.modified_terminals import is_modified_terminal, analyse_modified_terminal

from uflacs.analysis.graph import build_graph
from uflacs.analysis.graph_vertices import build_scalar_graph_vertices
from uflacs.analysis.graph_rebuild import rebuild_with_scalar_subexpressions
from uflacs.analysis.graph_dependencies import compute_dependencies, mark_active, mark_image
from uflacs.analysis.graph_ssa import compute_dependency_count, invert_dependencies
#from uflacs.analysis.graph_ssa import default_cache_score_policy, compute_cache_scores, allocate_registers

from uflacs.analysis.factorization import compute_argument_factorization


def build_scalar_graph(expressions):
    """Build list representation of expression graph covering the given expressions.

    TODO: Renaming, refactoring and cleanup of the graph building algorithms used in here
    """

    # Build the initial coarse computational graph of the expression
    G = build_graph(expressions)

    assert len(expressions) == 1, "Multiple expressions in graph building needs more work from this point on."

    # Build more fine grained computational graph of scalar subexpressions
    # TODO: Make it so that
    #   expressions[k] <-> NV[nvs[k][:]],
    #   len(nvs[k]) == value_size(expressions[k])
    scalar_expressions = rebuild_with_scalar_subexpressions(G)

    assert len(scalar_expressions) == sum(product(expr.ufl_shape) for expr in expressions)

    # Build new list representation of graph where all vertices of V represent single scalar operations
    e2i, V, target_variables = build_scalar_graph_vertices(scalar_expressions)

    return e2i, V, target_variables


def compute_expr_ir(expressions):
    """FIXME: Refactoring in progress!

    TODO:
    Work for later::

        - Apply some suitable renumbering of vertices and corresponding arrays prior to returning

        - Allocate separate registers for each partition
          (but e.g. argument[iq][i0] may need to be accessible in other loops)

        - Improve register allocation algorithm

        - Take a list of expressions as input to compile several expressions in one joined graph
          (e.g. to compile a,L,M together for nonlinear problems)

    """
    # Wrap in list if we only get one expression
    if not isinstance(expressions, list):
        expressions = [expressions]

    # TODO: Can we merge these three calls to something more efficient overall?
    # Build scalar list-based graph representation
    e2i, V, target_variables = build_scalar_graph(expressions)

    # Compute sparse dependency matrix
    dependencies = compute_dependencies(e2i, V)

    # Compute factorization of arguments
    argument_factorization, modified_arguments, V, target_variables, dependencies = \
        compute_argument_factorization(V, target_variables, dependencies)

    # Store modified arguments in analysed form
    for i in range(len(modified_arguments)):
        modified_arguments[i] = analyse_modified_terminal(modified_arguments[i])

    # --- Various dependency analysis ---

    # Count the number of dependencies every subexpr has
    depcount = compute_dependency_count(dependencies)

    # Build the 'inverse' of the sparse dependency matrix
    inverse_dependencies = invert_dependencies(dependencies, depcount)

    # Mark subexpressions of V that are actually needed for final result
    active, num_active = mark_active(dependencies, target_variables)

    # Build set of modified_terminal indices into factorized_vertices
    modified_terminal_indices = [i for i, v in enumerate(V)
                                 if is_modified_terminal(v)]

    # Build piecewise/varying markers for factorized_vertices
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

    # Build IR for the given expressions
    expr_ir = {}

    # Core expression graph:
    expr_ir["V"] = V                               # (array) V-index -> UFL subexpression
    expr_ir["target_variables"] = target_variables  # (array) Flattened input expression component index -> V-index

    # Result of factorization:
    expr_ir["modified_arguments"] = modified_arguments         # (array) MA-index -> UFL expression of modified arguments
    expr_ir["argument_factorization"] = argument_factorization  # (dict) tuple(MA-indices) -> V-index of monomial factor

    # TODO: More structured MA organization?
    #modified_arguments[rank][block][entry] -> UFL expression of modified argument
    #dofranges[rank][block] -> (begin, end)
    # or
    #modified_arguments[rank][entry] -> UFL expression of modified argument
    #dofrange[rank][entry] -> (begin, end)
    #argument_factorization: (dict) tuple(MA-indices (only relevant ones!)) -> V-index of monomial factor
    # becomes
    #argument_factorization: (dict) tuple(entry for each(!) rank) -> V-index of monomial factor ## doesn't cover intermediate f*u in f*u*v!

    # Dependency structure of graph:
    expr_ir["modified_terminal_indices"] = modified_terminal_indices  # (array) list of V-indices to modified terminals
    #expr_ir["dependencies"] = dependencies                           # (CRSArray) V-index -> direct dependency V-index list
    #expr_ir["inverse_dependencies"] = inverse_dependencies           # (CRSArray) V-index -> direct dependee V-index list

    # Metadata about each vertex
    #expr_ir["active"] = active       # (array) V-index -> bool
    expr_ir["piecewise"] = piecewise  # (array) V-index -> bool
    expr_ir["varying"] = varying     # (array) V-index -> bool

    return expr_ir

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

