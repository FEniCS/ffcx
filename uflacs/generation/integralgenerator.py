

from ufl.common import product
from uflacs.utils.log import debug, info, warning, error, uflacs_assert
from uflacs.codeutils.format_code_structure import format_code_structure, Indented, Block, ForRange
from uflacs.geometry.default_names import names
from uflacs.codeutils.languageformatter import CppStatementFormatterRules
langfmt = CppStatementFormatterRules()


""" # FFC language formatter includes, some may be useful:
from ufl.common import component_to_index
from ufl.permutation import build_component_numbering
from ufl.algorithms import MultiFunction

from uflacs.utils.log import uflacs_assert, warning, error

# TODO: The organization of code utilities is a bit messy...
from uflacs.codeutils.cpp_format import CppFormatterRulesCollection
from uflacs.geometry.default_names import names
from uflacs.backends.ffc.ffc_statement_formatter import langfmt
from uflacs.backends.ffc.ffc_statement_formatter import (format_element_table_access, format_entity_name)
from uflacs.elementtables.table_utils import derivative_listing_to_counts, flatten_component
"""


"""
TODO: Implement all steps in generate_expression_body. See in particular these functions:

Pre quadrature:
  tfmt.define_output_variables_reset()  # DONE: generate_tensor_reset
  tfmt.define_piecewise_geometry()      # FIXME: generate_pre_quadrature_loops
  tfmt.define_piecewise_coefficients()  # FIXME: ditto
  tfmt.define_registers(num_registers)

In quadrature:
  tfmt.define_coord_loop()
  tfmt.define_coord_vars()
  tfmt.define_coord_dependent_geometry()
  tfmt.define_coord_dependent_coefficients()

Argument loops:
  for ac in range(rank):
    tfmt.define_argument_for_loop(ac)
    tfmt.define_argument_loop_vars(ac)

Post quadrature:

"""

class IntegralGenerator(object):
    def __init__(self, ir):
        # Store ir
        self.ir = ir

        # Consistency check on quadrature rules
        nps1 = sorted(ir["uflacs_ir"]["expr_ir"].keys())
        nps2 = sorted(ir["quadrature_rules"].keys())
        if nps1 != nps2:
            uflacs_warning("Got different num_points for expression irs and quadrature rules:\n{0}\n{1}".format(
                nps1, nps2))

        # TODO: Populate these
        self._using_names = set()
        self._includes = set(("#include <cstring>",
                              "#include <cmath>"))

    def generate_using_statements(self):
        return ["using %s;" % name for name in sorted(self._using_names)]

    def get_includes(self):
        return sorted(self._includes)

    def generate(self):
        """Generate entire tabulate_tensor body.

        Assumes that the code returned from here will be wrapped in a context
        that matches a suitable version of the UFC tabulate_tensor signatures.
        """
        parts = []
        parts += [self.generate_using_statements()]
        parts += [self.generate_quadrature_tables()]
        parts += [self.generate_element_tables()]
        parts += [self.generate_tensor_reset()]
        parts += [self.generate_piecewise_partition()]
        parts += [self.generate_quadrature_loops()]
        parts += [self.generate_finishing_statements()]
        return format_code_structure(Indented(parts))

    def generate_quadrature_tables(self):
        "Generate static tables of quadrature points and weights."
        parts = []

        qrs = self.ir["quadrature_rules"]
        if qrs:
            parts += ["// Section for quadrature weights and points"]

        for num_points in sorted(qrs):
            weights = qrs[num_points][0]
            points = qrs[num_points][1]

            # Size of quadrature points depends on context, assume this is correct:
            pdim = len(points[0])

            wname = "%s%d" % (names.weights, num_points)
            pname = "%s%d" % (names.points, num_points)

            # TODO: Improve langfmt.array_decl to handle the {} wrappers here
            weights = "{ %s }" % langfmt.precision_floats(weights)
            points = "{ %s }" % langfmt.precision_floats(x for p in points for x in p)

            parts += [langfmt.array_decl("static const double", wname, num_points, weights)]
            parts += [langfmt.array_decl("static const double", pname, num_points*pdim, points)]
            parts += [""]

        return parts

    def generate_element_tables(self):
        "Generate static tables with precomputed element basis function values in quadrature points."
        parts = []
        parts += [langfmt.comment("Section for precomputed element basis function values")]
        expr_irs = self.ir["uflacs_ir"]["expr_ir"]
        for num_points in sorted(expr_irs):
            tables = expr_irs[num_points]["unique_tables"]
            comment = "Definitions of {0} tables for {0} quadrature points".format(len(tables), num_points)
            parts += [langfmt.comment(comment)]
            for name in sorted(tables):
                table = tables[name]
                if product(table.shape) > 0:
                    # TODO: Move to langfmt:
                    parts += [ArrayDecl("static const double", name, table.shape, table)]
        return parts

    def generate_tensor_reset(self):
        "Generate statements for resetting the element tensor to zero."

        # Compute tensor size
        if self.ir["domain_type"] == "interior_facet":
            A_shape = self.ir["prim_idims"]
        else:
            A_shape = [2*n for n in self.ir["prim_idims"]]
        A_size = product(A_shape)

        # TODO: Move to langfmt:
        def memzero(ptrname, size):
            return "memset({ptrname}, 0, {size} * sizeof(*{ptrname}));".format(ptrname=ptrname, size=size)

        parts = []
        parts += [langfmt.comment("Reset element tensor")]
        parts += [memzero(names.A, A_size)]
        #parts += [langfmt.memzero(names.A, A_size)]
        parts += [""]
        return parts

    def generate_quadrature_loops(self):
        "Generate all quadrature loops."
        parts = []
        for num_points in sorted(self.ir["uflacs"]["expr_ir"]):
            body = self.generate_quadrature_body(num_points)
            parts += [ForRange(names.iq, 0, num_points, body=body)]
        return parts

    def dummy(self):
        expr_ir = self.ir["uflacs_ir"]["expr_ir"][num_points]

        # Core expression graph:
        expr_ir["V"]
        expr_ir["target_variables"]

        # Result of factorization:
        expr_ir["modified_arguments"]
        expr_ir["argument_factorization"]

        # Metadata and dependency information:
        expr_ir["active"]
        expr_ir["dependencies"]
        expr_ir["inverse_dependencies"]
        expr_ir["modified_terminal_indices"]
        expr_ir["spatially_dependent_indices"]

        # Table names and dofranges
        expr_ir["modified_argument_table_ranges"]
        expr_ir["modified_terminal_table_ranges"]

    def generate_quadrature_body(self, num_points):
        """
        """
        parts = []
        parts += ["// Quadrature loop body setup {0}".format(num_points)]
        parts += self.generate_varying_partition(num_points)

        expr_ir = self.ir["uflacs_ir"]["expr_ir"][num_points]
        # tuple(modified_argument_indices) -> code_index
        AF = expr_ir["argument_factorization"]
        # modified_argument_index -> (tablename, dofbegin, dofend)
        MATR = expr_ir["modified_argument_table_ranges"]

        # Generate pre-argument loop computations # FIXME: What is this?
        for mas in sorted(AF):
            dofblock = tuple(MATR[ma][1:3] for ma in mas)
            ssa_index = AF[mas]

            # FIXME: Generate code for f*D and store reference to it with mas
            #f = self._ssa[ssa_index]
            #vname = name_from(mas)
            #vcode = code_from(ssa_index)
            #code += ["%s = (%s) * D;" % (vname, vcode)]
            #factors[mas] = vname

        # Nested argument loops and accumulation into element tensor
        parts += self.generate_quadrature_body_dofblocks(num_points)

        return parts

    def generate_quadrature_body_dofblocks(self, num_points, outer_dofblock=()):
        parts = []

        # The loop level iarg here equals the argument count (in renumbered >= 0 format)
        iarg = len(outer_dofblock)
        if iarg == self.ir["rank"]:
            # At the innermost argument loop level we accumulate into the element tensor
            parts += [self.generate_integrand_accumulation(num_points, outer_dofblock)]
            return parts
        assert iarg < self.ir["rank"]

        expr_ir = self.ir["uflacs_ir"]["expr_ir"][num_points]
        # tuple(modified_argument_indices) -> code_index
        AF = expr_ir["argument_factorization"]

        # modified_argument_index -> (tablename, dofbegin, dofend)
        MATR = expr_ir["modified_argument_table_ranges"]

        # Find dofranges at this loop level iarg starting with outer_dofblock
        dofranges = sorted(MATR[mas[iarg]][1:3] for mas in AF
                           if all(MATR[mas[i]][1:3] == outer_dofblock[i]
                                  for i in xrange(iarg)))

        # Build loops for each dofrange
        for dofrange in dofranges:
            dofblock = outer_dofblock + (dofrange,)
            body = []

            # FIXME: Argument partitions should be computed outside of the dofblock loops,
            #        in separate argument dofrange loops, which reduces n*n computations to n.
            # Generate code partition for dofblock at this loop level
            body += self.generate_argument_partition(num_points, iarg, dofblock[iarg])

            # Generate nested inner loops (only triggers for forms with two or more arguments
            body += self.generate_quadrature_body_dofblocks(num_points, dofblock)

            # Wrap setup, subloops, and accumulation in a loop for this level
            idof = "{name}{num}".format(name=names.ia, num=level)
            parts += [ForRange(idof, dofrange[0], dofrange[1], body=body)]
        return parts

    def generate_piecewise_partition(self): # FIXME: Generate 'piecewise' partition here
        """Generate statements prior to the quadrature loop.

        This mostly includes computations involving piecewise constant geometry and coefficients.
        """
        parts = []
        parts += ["// FIXME: Section for piecewise constant computations"]

        # FIXME: Get partition associated with num_points
        p = None
        #p = self._dofblock_partition_???
        #p = self._partitions["piecewise"]

        if p is not None:
            parts += generate_partition_assignments(p) # FIXME: Generate partition computation here
        else:
            parts += ["// FIXME: sp[...] = ...;"] # TODO: Remove this mock code

        return parts

    def generate_varying_partition(self, num_points): # FIXME: Generate 'varying' partition
        parts = []
        parts += ["// FIXME: Section for geometrically varying computations"]

        # FIXME: Get partition associated with num_points
        p = None
        #p = self._dofblock_partition_???.get(num_points)
        #p = self._partitions["varying"].get(num_points)

        if p is not None:
            parts += generate_partition_assignments(p) # FIXME: Generate partition computation here
        else:
            parts += ["// FIXME: sv[...] = ...;"] # TODO: Remove this mock code

        return parts

    def generate_argument_partition(self, num_points, iarg, dofrange): # FIXME: Define better and implement
        """Generate code for the partition corresponding to arguments 0..iarg within given dofblock."""
        parts = []

        # FIXME: Get partition associated with (num_points, iarg, dofrange)
        p = None
        #p = self._dofblock_partition[num_points][iarg].get(dofrange)
        #p = self._partitions["argument"][num_points][iarg].get(dofrange)

        if p is not None:
            parts += generate_partition_assignments(p) # FIXME: Generate partition computation here
        else:
            parts += ["// FIXME: sa[...] = ...; // {0} x {1}".format(iarg, dofrange)] # TODO: Remove this mock code

        return parts

    def generate_integrand_accumulation(self, num_points, dofblock): # FIXME: Implement!
        parts = []

        # FIXME: Get partition associated with (num_points, dofblock)
        p = None
        #p = self._dofblock_integrand_partition[num_points][dofblock]
        #p = self._partitions["integrand"][num_points].get(dofblock)

        if p is not None:
            parts += generate_partition_accumulations(p) # FIXME: Generate partition computation here
        else:
            parts += ["// FIXME: A[{0}] += f * v * D;".format(dofblock)] # TODO: Remove this mock code

        # Corresponding code from compiler.py:
        #final_variable_names = [format_register_variable(p, r) for r in ir["target_registers"]]
        #assign_to_variables = tfmt.output_variable_names(len(final_variable_names))
        #scaling_factor = tfmt.accumulation_scaling_factor()
        #parts += list(format_scaled_additions(zip(assign_to_variables,
        #                                          final_variable_names),
        #                                          scaling_factor))

        return parts

    def generate_finishing_statements(self):
        """Generate finishing statements.

        This includes assigning to output array if there is no integration.
        """
        parts = []

        if not self.ir["quadrature_rules"]: # Rather check ir["domain_type"]?
            # TODO: Implement for expression support
            error("Expression generation not implemented yet.")
            # TODO: If no integration, assuming we generate an expression, and assign results here
            # Corresponding code from compiler.py:
            #assign_to_variables = tfmt.output_variable_names(len(final_variable_names))
            #parts += list(format_assignments(zip(assign_to_variables, final_variable_names)))

        return parts


class CellIntegralGenerator(IntegralGenerator):
    pass

class ExteriorFacetIntegralGenerator(IntegralGenerator):
    pass

class InteriorFacetIntegralGenerator(IntegralGenerator):
    pass

class DiracIntegralGenerator(IntegralGenerator):
    pass

class QuadratureIntegralGenerator(IntegralGenerator):
    pass
