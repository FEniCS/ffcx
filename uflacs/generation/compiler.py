
import numpy
from uflacs.utils.tictoc import TicToc
from uflacs.utils.log import error, uflacs_assert

from uflacs.analysis.graph import build_graph
from uflacs.analysis.graph_vertices import build_scalar_graph_vertices
from uflacs.analysis.graph_rebuild import rebuild_scalar_e2i
from uflacs.analysis.graph_ssa import (compute_dependencies,
                                         mark_active,
                                         mark_partitions,
                                         compute_dependency_count,
                                         invert_dependencies,
                                         default_cache_score_policy,
                                         compute_cache_scores,
                                         allocate_registers)
from uflacs.analysis.factorization import compute_argument_factorization, rebuild_scalar_graph_from_factorization

from uflacs.codeutils.expr_formatter import ExprFormatter
from uflacs.codeutils.element_tensor_formatter import build_loops
from uflacs.codeutils.format_lines import format_assignments, format_scaled_additions

from uflacs.generation.partitions import build_partition_labels

# FFC compiler calls:
# 1) compile_expression_partitions
# 2) generate_code_from_ssa
# 3) generate_expression_body


def compile_expression_partitions(expressions, parameters):
    """FIXME: Refactoring in progress!

    TODO: assuming more symbolic preprocessing
    - Make caller apply grad->localgrad+jacobianinverse
    - Make caller apply coefficient mapping
    - Make caller apply piola mappings

    TODO:
    - Build modified_argument_factors
    - Build modified_argument_blocks
    - Rewrite build_partition_labels in terms of modified_argument_blocks instead of dofblocks
    - Rewrite mark_partitions to use seeds from partition_labels

    - Place partition_labels, argument_factors, etc. in ir
    - Use new ir information in code generation

    - Apply some suitable renumbering of vertices and corresponding arrays prior to returning
    - Allocate separate registers for each partition
      (but e.g. argument[iq][i0] may need to be accessible in other loops)
    - Improve register allocation algorithm

    - What about conditionals?

    - Take a list of expressions as input to compile several expressions in one joined graph
      (e.g. to compile a,L,M together for nonlinear problems)
    """
    # Timing object
    tt = TicToc('compile_expression_partitions')

    # Wrap in list if we only get one expression
    if not isinstance(expressions, list):
        expressions = [expressions]


    # Build the initial coarse computational graph of the expression
    tt.step('build_graph')
    G = build_graph(expressions)


    # Build more fine grained computational graph of scalar subexpressions
    tt.step('rebuild_scalar_e2i')
    unused_e2i, NV, unused_W, modified_terminals, nvs = rebuild_scalar_e2i(G, DEBUG=False)
    # Target expression is NV[nvs[:]].
    # TODO: Make it so that expressions[k] <-> NV[nvs[k][:]], len(nvs[k]) == value_size(expressions[k])


    # Straigthen out V to represent single operations
    tt.step('build_scalar_graph_vertices')
    se2i, SV, target_variables = build_scalar_graph_vertices([NV[s] for s in nvs])


    # Compute sparse dependency matrix
    tt.step('compute_dependencies')
    dependencies = compute_dependencies(se2i, SV)


    # Compute factorization of arguments
    tt.step('compute_dependencies')
    # TODO: Fix factorization algorithm for multiple target variables! Or is there any point?
    if parameters["enable_factorization"] and len(target_variables) == 1:

        # AV, FV, IM
        argument_factors, factorized_vertices, argument_factorization = \
            compute_argument_factorization(SV, target_variables, dependencies)


        # Rebuild some graphs from factorization
        SV, se2i, dependencies = rebuild_scalar_graph_from_factorization(\
            argument_factors, factorized_vertices, argument_factorization)


        # TODO: target_variables for non-scalar or multiple expressions
        target_variables = [len(SV)-1]


    if 0: # FIXME: Use this!

        # XXX: FIXME: find dofrange for each argument_factors value
        #dofranges[np][iarg] = [...]
        #argument_factors2[np][iarg] = [ma for ma in argument_factors if ma involves iarg and is needed for np] # FIXME: Can we avoid explicit dofranges in here by using ma indices?

        # XXX: FIXME: find dofblock for each argument_factorization item
        #dofblocks[np] = [...]
        #argument_blocks[np] = sorted(argument_factorization.keys()) # FIXME: Can we avoid explicit dofblocks in here by using ma index tuples?

        partition_labels, max_p = build_partition_labels(rank, num_points,
                                                         #argument_factors2, argument_blocks)
                                                         dofranges, dofblocks)

        # TODO: reorder vertices grouped by partition number
        # FIXME: change mark_partitions to use this numbering instead
        # FIXME: change code generation to use this numbering instead


    # Mark subexpressisons that are actually needed for final result
    tt.step('mark_active')
    max_symbol = len(SV)
    active, num_active = mark_active(max_symbol, dependencies, target_variables)

    # TODO: Mark subexpressions with inside/outside quadrature loop


    # TODO: Renumber expressions to group inside/outside and argument factors?


    # Mark subexpressions with which loop they belong inside
    tt.step('mark_partitions')
    # NB! This one assumes that argument mapping has been applied in SV, which it currently hasn't.
    # However this code may be superseeded soon anyway?
    partitions = mark_partitions(SV, active, dependencies, rank)



    # Count the number of dependencies every subexpr has
    tt.step('compute_dependency_count')
    depcount = compute_dependency_count(dependencies)


    # Build the 'inverse' of the sparse dependency matrix
    tt.step('invert_dependencies')
    inverse_dependencies = invert_dependencies(dependencies, depcount)


    # Use heuristics to mark the usefulness of storing every subexpr in a variable
    tt.step('compute_cache_scores')
    scores = compute_cache_scores(SV,
                                  active,
                                  dependencies,
                                  inverse_dependencies,
                                  partitions,
                                  cache_score_policy=default_cache_score_policy)

    if 0:
        print '\n'*5
        print "SV:"
        print '\n'.join('  {}: {}'.format(i,s) for i,s in enumerate(SV))
        print "scores:"
        print '\n'.join('  {}: {}'.format(i,s) for i,s in enumerate(scores))
        print '\n'*5

    # Allocate variables to store subexpressions in
    tt.step('allocate_registers')
    allocations = allocate_registers(active, partitions, target_variables,
                                     scores, int(parameters["max_registers"]), int(parameters["score_threshold"]))
    target_registers = [allocations[r] for r in target_variables]
    num_registers = sum(1 if x >= 0 else 0 for x in allocations)
    # TODO: If we renumber we can allocate registers separately for each partition, which is probably a good idea.

    # Print timing
    tt.stop()
    if parameters["enable_profiling"]:
        print "Profiling results:"
        print tt

    partitions_ir = {}
    partitions_ir["num_registers"] = num_registers
    partitions_ir["SV"] = SV
    partitions_ir["active"] = active
    partitions_ir["partitions"] = partitions
    partitions_ir["allocations"] = allocations
    partitions_ir["target_registers"] = target_registers
    partitions_ir["terminals"] = modified_terminals
    return partitions_ir
