
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


def build_partition_labels(rank, num_points, dofranges, dofblocks):

    # - build numbering of partitions:
    #   - piecewise: 1
    #   - for each np:
    #     - varying: 1 + conditionals
    #       - conditional: for each conditional expression case?
    #     - argument: for each argument for each dofrange
    #     - integrand: for each dofblock

    # Build a numbering of partitions, where p2 > p1 if p1 occurs before p2 in the same code path
    partition_labels = {}
    p = 1

    # "piecewise": partition of expressions independent of quadrature and argument loops
    partition_labels["piecewise"] = p; p += 1

    # "varying", [np]: partition of expressions dependent on np quadrature but independent of argument loops
    partition_labels["varying"] = {}
    for np in num_points:
        partition_labels["varying"][np] = p; p += 1

    # "argument", [np][iarg][dofrange]: partition depending on this dofrange of argument iarg
    partition_labels["argument"] = {}
    for np in num_points:
        partition_labels["argument"][np] = {}
        for iarg in range(rank):
            partition_labels["argument"][np][iarg] = {}
            for dofrange in dofranges[np][iarg]:
                partition_labels["argument"][np][iarg][dofrange] = p; p += 1

    # "integrand", [np][dofblock]: partition depending on this dofblock
    partition_labels["integrand"] = {}
    for np in num_points:
        partition_labels["integrand"][np] = {}
        for dofblock in dofblocks[np]:
            partition_labels["integrand"][np][dofblock] = p; p += 1

    # TODO: Add partition labels for conditional cases?
    #       Maybe nested label tuples? This is not thought through...
    #partition_labels["conditional"][np] = {}
    #partition_labels["conditional"][np][conditional_true] = p; p += 1
    #partition_labels["conditional"][np][conditional_false] = p; p += 1

    max_p = p
    return partition_labels, max_p

