
import numpy
from uflacs.utils.tictoc import TicToc
from uflacs.utils.log import error, uflacs_assert

from uflacs.analysis.graph import build_graph
from uflacs.analysis.graph_vertices import build_scalar_graph_vertices
from uflacs.analysis.graph_rebuild import rebuild_scalar_e2i
from uflacs.analysis.graph_dependencies import (compute_dependencies,
                                                mark_active)

from uflacs.analysis.graph_ssa import (mark_partitions,
                                         compute_dependency_count,
                                         invert_dependencies,
                                         default_cache_score_policy,
                                         compute_cache_scores,
                                         allocate_registers)
from uflacs.analysis.factorization import compute_argument_factorization, rebuild_scalar_graph_from_factorization
from uflacs.analysis.dependency_handler import DependencyHandler

from uflacs.codeutils.expr_formatter import ExprFormatter
from uflacs.codeutils.element_tensor_formatter import build_loops
from uflacs.codeutils.format_lines import format_assignments, format_scaled_additions


# TODO: Move these to [target_]expr_formatter:
def format_assignment_statement(lhs, rhs):
    #return format_assignments([(lhs,rhs)]) # Use this?
    return "%s = %s;" % (lhs, rhs)

def format_register_variable(p, r):
    return "s[%d]" % (r,) # TODO: Maybe make it "s%d[%d]" % (p, r)


def generate_partition_code(SV, allocations, active, partitions, p, language_formatter):
    # This is a transformer that collects terminal modifiers
    # and delegates formatting to the language_formatter
    expr_formatter = ExprFormatter(language_formatter, {})
    code = []
    # Generate assignment code for each expression with an allocated register in partition p
    for i, v in enumerate(SV):
        vreg = allocations[i]
        if active[i] and partitions[i] == p and vreg >= 0:
            vname = format_register_variable(p, vreg)
            vcode = expr_formatter.visit(v)
            assignment = format_assignment_statement(vname, vcode)
            code.append(assignment)
            expr_formatter.variables[v] = vname
    return code

def generate_code_from_ssa(partitions_ir, language_formatter):
    # Fetch stuff from ir:
    SV = partitions_ir["SV"]
    active = partitions_ir["active"]
    partitions = partitions_ir["partitions"]
    allocations = partitions_ir["allocations"]
    target_registers = partitions_ir["target_registers"]

    # Find partition range, skip negative which means inactive vertices
    min_p = max(0,min(partitions))
    max_p = max(partitions)

    # Handle partitions one at a time, in order
    partition_codes = []
    for p in range(min_p, max_p+1):
        code = generate_partition_code(SV, allocations, active, partitions, p, language_formatter)
        partition_codes.append((p,code))

    p = max_p # TODO: Add partitions to target_registers, or is this fine?
    final_variable_names = [format_register_variable(p, r) for r in target_registers]
    return partition_codes, final_variable_names



# FIXME: Replace generate_expression_body and FFCStatementFormatter with the new IntegralGenerator class
def generate_expression_body(target_statement_formatter, partition_codes, rank,
                             final_variable_names, num_registers):
    "Join partitions with target specific loops and declarations."
    # Use shorter name below
    tfmt = target_statement_formatter

    # Make a shallow copy of dict, we consume the dict entries below
    partition_codes = dict(partition_codes)

    # Build loop structure (intended for tabulate_tensor)
    loops = []
    definitions = []
    partitions = []

    def update_loops(loop, defs, p):
        loops.append(loop)
        definitions.append(defs)
        partitions.append(partition_codes.get(p,""))
        if p in partition_codes:
            del partition_codes[p]

    # TODO: This partition numbering is maybe too "magic"? At least make named constants!
    # --- Partition 0: independent of x and arguments.
    p = 0
    if 1:
        piecewise_defs = []
        piecewise_defs += tfmt.define_output_variables_reset()
        piecewise_defs += tfmt.define_piecewise_geometry()
        piecewise_defs += tfmt.define_piecewise_coefficients()
        # TODO: Make registers separately for each partition level?
        if num_registers:
            piecewise_defs += tfmt.define_registers(num_registers)
        piecewise_defs += [""]
        update_loops(None, piecewise_defs, p)

    # --- Partition 1: x dependency
    # FIXME: Separate cleanly between:
    # - x assumed known
    # - xi assumed known
    # - xi defined by inserted quadrature loop
    # - x defined by inserted quadrature loop
    # ... Probably best to configure tfmt to know this beforehand some place?
    p = 1
    if 1:
        coord_loop = tfmt.define_coord_loop() # Can be None
        coord_dependent_defs = []
        coord_dependent_defs += tfmt.define_coord_vars()
        coord_dependent_defs += tfmt.define_coord_dependent_geometry()
        update_loops(coord_loop, coord_dependent_defs, p)

    # --- Partition 2: coord and coefficient dependency
    p = 2
    if 1:
        w_dependent_defs = []
        w_dependent_defs += tfmt.define_coord_dependent_coefficients()
        update_loops(None, w_dependent_defs, p)

    # --- Partitions 3...3+rank-1: argument function dependency
    poffset = 3
    for ac in range(rank):
        p = poffset + ac
        update_loops(tfmt.define_argument_for_loop(ac),
                     tfmt.define_argument_loop_vars(ac),
                     p)

    # --- Final partition: final assignments
    p = poffset + rank
    if 1:
        assign_to_variables = tfmt.output_variable_names(len(final_variable_names))
        scaling_factor = tfmt.accumulation_scaling_factor()
        if scaling_factor is None:
            final_statements = list(format_assignments(zip(assign_to_variables,
                                                           final_variable_names)))
        else:
            final_statements = list(format_scaled_additions(zip(assign_to_variables,
                                                                final_variable_names),
                                                            scaling_factor))
        update_loops(None, final_statements, p)

    # --- Should be nothing left now
    assert not partition_codes

    # Stitch it together
    code = build_loops(loops, definitions, partitions)

    return code

def generate_expression_code(partitions_ir, form_argument_mapping, object_names,
                             create_language_formatter, create_statement_formatter):
    "Core of toy expression compiler."

    # Create an object to track dependencies across other components
    dependency_handler = DependencyHandler(partitions_ir["terminals"],
                                           form_argument_mapping,
                                           object_names)

    #dependency_handler = DependencyHandler(partitions_ir["terminals"],
    #                                       uflacs_ir["function_replace_map"],
    #                                       uflacs_ir["argument_names"],
    #                                       uflacs_ir["coefficient_names"])

    # This formatter is a multifunction implementing target specific formatting rules
    language_formatter = create_language_formatter(dependency_handler, partitions_ir)

    # Create a formatter for blocks of statements
    statement_formatter = create_statement_formatter(dependency_handler, partitions_ir)

    # Generate code partitions from ir
    partition_codes, final_variable_names = generate_code_from_ssa(partitions_ir, language_formatter)

    # Generate full code from snippets
    #rank = len(target_statement_formatter._dependency_handler.mapped_arguments)
    rank = len(uflacs_ir["argument_names"])
    code = generate_expression_body(statement_formatter,
                                    partition_codes,
                                    rank,
                                    final_variable_names,
                                    partitions_ir["num_registers"])

    # Leave final formatting to the caller
    coefficient_names = sorted(dependency_handler.coefficient_names.values())
    return code, coefficient_names
