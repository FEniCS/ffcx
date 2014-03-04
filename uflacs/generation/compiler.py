
import numpy
from uflacs.utils.tictoc import TicToc
from uflacs.utils.log import error, uflacs_assert

from uflacs.analysis.graph import build_graph
from uflacs.analysis.graph_vertices import build_scalar_graph_vertices
from uflacs.analysis.graph_rebuild import rebuild_scalar_e2i
from uflacs.analysis.graph_dependencies import (compute_dependencies,
                                                mark_active, mark_image)
from uflacs.analysis.graph_ssa import (mark_partitions,
                                       compute_dependency_count,
                                       invert_dependencies,
                                       default_cache_score_policy,
                                       compute_cache_scores,
                                       allocate_registers)
from uflacs.analysis.factorization import compute_argument_factorization, rebuild_scalar_graph_from_factorization
from uflacs.analysis.modified_terminals import is_modified_terminal
from uflacs.codeutils.expr_formatter import ExprFormatter
from uflacs.codeutils.element_tensor_formatter import build_loops
from uflacs.codeutils.format_lines import format_assignments, format_scaled_additions

from uflacs.generation.partitions import build_partition_labels


# TODO: Move to utils:
def format_sequence(sequence):
    return '\n'.join("{0}".format(v) for v in sequence)
def format_enumerated_sequence(sequence):
    return '\n'.join("{0}: {1}".format(i, v) for i,v in enumerate(sequence))
def format_mapping(mapping):
    return '\n'.join("{0}: {1}".format(k, v) for k,v in mapping.items())


# FFC compiler calls:
# 1) compile_expression_partitions
# 2) generate_code_from_ssa
# 3) generate_expression_body


# TODO: Refactoring and cleanup of the algorithms hidden behind this clean interface
# TODO: Clean up graph building: remove modified_terminals returnee from there, can build that elsewhere
def build_scalar_graph(expressions):
    # Build the initial coarse computational graph of the expression
    G = build_graph(expressions)

    assert len(expressions) == 1, "Multiple expressions in graph building needs more work from this point on."

    # Build more fine grained computational graph of scalar subexpressions
    unused_e2i, NV, unused_W, modified_terminals, nvs = rebuild_scalar_e2i(G, DEBUG=False)

    # Target expression is NV[nvs[:]].
    # TODO: Make it so that expressions[k] <-> NV[nvs[k][:]], len(nvs[k]) == value_size(expressions[k])

    # Get scalar target expressions
    scalar_expressions = [NV[s] for s in nvs]

    # Straigthen out V to represent single operations
    e2i, V, target_variables = build_scalar_graph_vertices(scalar_expressions)

    return e2i, V, target_variables, modified_terminals


def compute_expr_ir(expressions, parameters):
    """FIXME: Refactoring in progress!

    TODO: assuming more symbolic preprocessing
    - Make caller apply grad->localgrad+jacobianinverse
    - Make caller apply coefficient mapping
    - Make caller apply piola mappings

    TODO:
    - Build modified_argument_blocks
    - Use new ir information in code generation

    Work for later:
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


    # Build scalar list-based graph representation
    tt.step('build_scalar_graph')
    e2i, V, target_variables, modified_terminals = build_scalar_graph(expressions)


    # Compute sparse dependency matrix
    tt.step('compute_dependencies')
    dependencies = compute_dependencies(e2i, V)

    # Compute factorization of arguments
    tt.step('compute_argument_factorization')
    argument_factorization, modified_arguments, V, target_variables, dependencies = \
        compute_argument_factorization(V, target_variables, dependencies)


    # Various dependency analysis
    tt.step('various dependency analysis')

    # Count the number of dependencies every subexpr has
    depcount = compute_dependency_count(dependencies)

    # Build the 'inverse' of the sparse dependency matrix
    inverse_dependencies = invert_dependencies(dependencies, depcount)

    # Mark subexpressions of V that are actually needed for final result
    active, num_active = mark_active(dependencies, target_variables)

    # Build set of modified_terminal indices into factorized_vertices
    modified_terminal_indices = [i for i,v in enumerate(V)
                                 if is_modified_terminal(v)]

    # Build piecewise/varying markers for factorized_vertices
    spatially_dependent_terminal_indices = [i for i in modified_terminal_indices
                                   if not V[i].is_cellwise_constant()]
    spatially_dependent_indices, num_spatial = mark_image(inverse_dependencies,
                                                          spatially_dependent_terminal_indices)

    # TODO: Inspection of spatially_dependent_indices shows that factorization is
    # needed for effective loop invariant code motion w.r.t. quadrature loop as well.
    # Postphoning that until everything is working fine again.
    # Core ingredients for such factorization would be:
    # - Flatten products of products somehow
    # - Sorting flattened product factors by loop dependency then by canonical ordering

    #rank = max(len(k) for k in argument_factorization.keys())
    #for i,a in enumerate(modified_arguments):
    #    iarg = a.number()
    #    #ipart = a.part()

    # Print timing
    tt.stop()
    if parameters["enable_profiling"]:
        print "Profiling results:"
        print tt

    # Build IR for the given expressions
    expr_ir = {}

    # Core expression graph:
    expr_ir["V"] = V
    expr_ir["target_variables"] = target_variables

    # Result of factorization:
    expr_ir["modified_arguments"] = modified_arguments
    expr_ir["argument_factorization"] = argument_factorization

    # Metadata and dependency information:
    expr_ir["active"] = active
    expr_ir["dependencies"] = dependencies
    expr_ir["inverse_dependencies"] = inverse_dependencies
    expr_ir["modified_terminal_indices"] = modified_terminal_indices
    expr_ir["spatially_dependent_indices"] = spatially_dependent_indices

    return expr_ir

def old_code_useful_for_optimization():

    # Use heuristics to mark the usefulness of storing every subexpr in a variable
    tt.step('compute_cache_scores')
    scores = compute_cache_scores(V,
                                  active,
                                  dependencies,
                                  inverse_dependencies,
                                  partitions, # TODO: Rewrite in terms of something else, this doesn't exist anymore
                                  cache_score_policy=default_cache_score_policy)

    # Allocate variables to store subexpressions in
    tt.step('allocate_registers')
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
