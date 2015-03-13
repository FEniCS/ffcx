
from six import iteritems, iterkeys
from six.moves import zip
from six.moves import xrange as range

from ufl.common import product
from ufl.classes import ConstantValue

from ffc.log import error

from uflacs.analysis.modified_terminals import analyse_modified_terminal, is_modified_terminal

from uflacs.codeutils.format_code import (format_code, Indented, Comment,
                                          ForRange, Block,
                                          VariableDecl, ArrayDecl, ArrayAccess,
                                          Assign, AssignAdd,
                                          Product)
from uflacs.codeutils.indexmapping import IndexMapping, AxisMapping
from uflacs.codeutils.expr_formatter import ExprFormatter


class IntegralGenerator(object):

    def __init__(self, ir, language_formatter, backend_access, backend_definitions):
        # Store ir
        self.ir = ir

        # Consistency check on quadrature rules
        nps1 = sorted(iterkeys(ir["uflacs"]["expr_ir"]))
        nps2 = sorted(iterkeys(ir["quadrature_rules"]))
        if nps1 != nps2:
            uflacs_warning("Got different num_points for expression irs and quadrature rules:\n{0}\n{1}".format(
                nps1, nps2))

        # Compute shape of element tensor
        if self.ir["integral_type"] == "interior_facet":
            self._A_shape = [2 * n for n in self.ir["prim_idims"]]
        else:
            self._A_shape = self.ir["prim_idims"]

        # TODO: Populate these with only what's needed
        self._using_names = set()
        self._includes = {
            "#include <cstring>",
            "#include <cmath>",
            "#include <boost/math/special_functions.hpp>",
        }

        # Formatter for backend agnostic expressions, delegating to given target language formatter
        self.expr_formatter = ExprFormatter(language_formatter, {})

        # Formatter for defining backend specific variables
        self.backend_definitions = backend_definitions

        # Formatter for accessing backend specific variables
        self.backend_access = backend_access

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
        parts += self.generate_using_statements()
        parts += self.backend_definitions.initial()
        parts += self.generate_quadrature_tables()
        parts += self.generate_element_tables()
        parts += self.generate_tensor_reset()

        # If we have integrals with different number of quadrature points,
        # we wrap each integral in a separate scope, avoiding having to
        # think about name clashes for now. This is a bit wasteful in that
        # piecewise quantities are not shared, but at least it should work.
        expr_irs = self.ir["uflacs"]["expr_ir"]
        all_num_points = sorted(expr_irs)
        for num_points in all_num_points:
            self.expr_formatter.variables = {}
            #self.ast_variables = int_array(FIXME)
            pp = self.generate_piecewise_partition(num_points)
            ql = self.generate_quadrature_loops(num_points)
            if len(all_num_points) > 1:
                parts += [Block([pp, ql])]
            else:
                parts += [pp, ql]

        parts += self.generate_finishing_statements()
        return format_code(Indented(parts))

    def generate_quadrature_tables(self):
        "Generate static tables of quadrature points and weights."
        parts = []

        # No quadrature tables for custom (given argument) or point (evaluation in single vertex)
        if self.ir["integral_type"] in ("custom", "vertex"):
            return parts

        qrs = self.ir["quadrature_rules"]
        if qrs:
            parts += ["// Section for quadrature weights and points"]

        for num_points in sorted(qrs):
            weights = qrs[num_points][0]
            points = qrs[num_points][1]

            # Size of quadrature points depends on context, assume this is correct:
            pdim = len(points[0])

            wname = self.backend_access.weights_array_name(num_points)
            pname = self.backend_access.points_array_name(num_points)

            weights = [self.backend_access.precision_float(w) for w in weights]
            points = [self.backend_access.precision_float(x) for p in points for x in p]

            parts += [ArrayDecl("static const double", wname, num_points, weights)]
            if pdim > 0:
                parts += [ArrayDecl("static const double", pname, num_points * pdim, points)]
            parts += [""]

        return parts

    def generate_element_tables(self):
        "Generate static tables with precomputed element basis function values in quadrature points."
        parts = []
        parts += [Comment("Section for precomputed element basis function values"),
                  Comment("Table dimensions: num_entities, num_points, num_dofs"),
                  ""]
        expr_irs = self.ir["uflacs"]["expr_ir"]
        for num_points in sorted(expr_irs):
            tables = expr_irs[num_points]["unique_tables"]
            comment = "Definitions of {0} tables for {1} quadrature points".format(len(tables), num_points)
            parts += [Comment(comment)]
            for name in sorted(tables):
                table = tables[name]
                if product(table.shape) > 0:
                    parts += [ArrayDecl("static const double", name, table.shape, table), ""]
        return parts

    def generate_tensor_reset(self):
        "Generate statements for resetting the element tensor to zero."

        # Could move this to codeutils or backend
        def memzero(ptrname, size):
            return "memset({ptrname}, 0, {size} * sizeof(*{ptrname}));".format(ptrname=ptrname, size=size)

        # Compute tensor size
        A_size = product(self._A_shape)
        A = self.backend_access.element_tensor_name()

        parts = []
        parts += [Comment("Reset element tensor")]
        parts += [memzero(A, A_size)]
        parts += [""]
        return parts

    def generate_quadrature_loops(self, num_points):
        "Generate all quadrature loops."
        parts = []
        body = self.generate_quadrature_body(num_points)
        iq = self.backend_access.quadrature_loop_index()
        if num_points == 1:
            parts += [Comment("Only 1 quadrature point, no loop"),
                      VariableDecl("const int", iq, 0), # TODO: Inject iq=0 in generated code instead of this line
                      Block(body)] # Wrapping in Block to avoid thinking about scoping issues
        else:
            parts += [ForRange(iq, 0, num_points, body=body)]
        return parts

    def generate_quadrature_body(self, num_points):
        """
        """
        parts = []
        parts += self.generate_varying_partition(num_points)
        if parts:
            parts = [Comment("Quadrature loop body setup (num_points={0})".format(num_points))] + parts + [""]

        # Compute single argument partitions outside of the dofblock loops
        for iarg in range(self.ir["rank"]):
            for dofrange in []:  # TODO: Move f*arg0 out here
                parts += self.generate_argument_partition(num_points, iarg, dofrange)

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

        expr_ir = self.ir["uflacs"]["expr_ir"][num_points]
        # tuple(modified_argument_indices) -> code_index
        AF = expr_ir["argument_factorization"]

        # modified_argument_index -> (tablename, dofbegin, dofend)
        MATR = expr_ir["modified_argument_table_ranges"]

        # Find dofranges at this loop level iarg starting with outer_dofblock
        dofranges = set()
        for mas in AF:
            mas_full_dofblock = tuple(MATR[j][1:3] for j in mas)
            if tuple(mas_full_dofblock[:iarg]) == tuple(outer_dofblock):
                dofrange = mas_full_dofblock[iarg]
                # Skip empty dofranges TODO: Possible to remove these and related code earlier?
                if dofrange[0] != dofrange[1]:
                    dofranges.add(dofrange)
        dofranges = sorted(dofranges)

        # Build loops for each dofrange
        for dofrange in dofranges:
            dofblock = outer_dofblock + (dofrange,)
            body = []

            # Generate nested inner loops (only triggers for forms with two or more arguments
            body += self.generate_quadrature_body_dofblocks(num_points, dofblock)

            # Wrap setup, subloops, and accumulation in a loop for this level
            idof = self.backend_access.argument_loop_index(iarg)
            parts += [ForRange(idof, dofrange[0], dofrange[1], body=body)]
        return parts

    def generate_partition(self, name, V, partition, table_ranges, num_points):  # TODO: Rather take list of vertices, not markers
        terminalcode = []
        assignments = []
        j = 0
        # print "Generating partition ", name, partition
        for i, p in enumerate(partition):
            if p:
                # TODO: Consider optimized ir here with markers for which subexpressions to store in variables.
                # This code just generates variables for _all_ subexpressions marked by p.

                v = V[i]

                if is_modified_terminal(v):
                    mt = analyse_modified_terminal(v)
                    if isinstance(mt.terminal, ConstantValue):
                        # Literal value
                        vaccess = self.expr_formatter.visit(v)
                    else:
                        # Backend specific modified terminal
                        vaccess = self.backend_access(mt.terminal, mt, table_ranges[i], num_points)
                        vdef = self.backend_definitions(mt.terminal, mt, table_ranges[i], vaccess)
                        # Store definitions of terminals in list
                        terminalcode += [vdef]
                else:
                    # Application of operator
                    # Count assignments so we get a new vname each time
                    vaccess = ArrayAccess(name, j)
                    j += 1
                    vcode = self.expr_formatter.visit(v)  # TODO: Generate ASTNode instead of str here?
                    # Store assignments of operator results in list
                    assignments += [Assign(vaccess, vcode)]

                # Store access string, either a variable name or an inlined expression
                # TODO: Can skip format_code if expr_formatter generates ASTNode
                self.expr_formatter.variables[v] = format_code(vaccess)

        parts = []
        # Compute all terminals first
        parts += terminalcode
        if j > 0:
            # Declare array large enough to hold all subexpressions we've emitted
            parts += [ArrayDecl("double", name, j)]
            # Then add all computations
            parts += assignments
        return parts

    # TODO: Rather take list of vertices, not markers
    def alternative_generate_partition(self, name, C, MT, partition, table_ranges, num_points):
        terminalcode = []
        assignments = []

        # C = input CRS representation of expression DAG
        # MT = input list/dict of modified terminals

        self.ast_variables = [None]*len(C) # FIXME: Create outside
        vertices = [i for i, p in enumerate(partition) if p] # TODO: Get this as input instead of partition?

        for i in vertices:
            row = C[i] # FIXME: Get this as input
            if len(row) == 1:
                # Modified terminal
                t, = row
                mt = MT[t] # FIXME: Get this as input

                if isinstance(mt.terminal, ConstantValue):
                    # Format literal value for the chosen language
                    vaccess = modified_literal_to_ast_node[tc](mt) # FIXME: Implement this mapping

                else:
                    # Backend specific modified terminal formatting
                    vaccess = self.backend_access(mt.terminal, mt, table_ranges[i], num_points)
                    vdef = self.backend_definitions(mt.terminal, mt, table_ranges[i], vaccess)

                    # Store definitions of terminals in list
                    terminalcode += [vdef]

            else:
                # Application of operator with typecode tc to operands with indices ops
                tc = mt[0]
                ops = mt[1:]

                # Get operand AST nodes
                opsaccess = [self.ast_variables[k] for k in ops]

                # Generate expression for this operator application
                vcode = typecode2astnode[tc](opsaccess) # FIXME: Implement this mapping

                store_this_in_variable = True # TODO: Don't store all subexpressions
                if store_this_in_variable:
                    # Count assignments so we get a new vname each time
                    vaccess = ArrayAccess(name, len(assignments))
                    # Store assignments of operator results in list
                    assignments += [Assign(vaccess, vcode)]
                else:
                    # Store the inlined expression
                    vaccess = vcode

            # Store access string, either a variable name or an inlined expression
            self.ast_variables[i] = vaccess

        parts = []
        # Compute all terminals first
        parts += terminalcode
        # Then add all computations
        if assignments:
            # Declare array large enough to hold all subexpressions we've emitted
            parts += [ArrayDecl("double", name, len(assignments))]
            parts += assignments
        return parts

    def generate_piecewise_partition(self, num_points):
        """Generate statements prior to the quadrature loop.

        This mostly includes computations involving piecewise constant geometry and coefficients.
        """
        parts = []
        expr_ir = self.ir["uflacs"]["expr_ir"][num_points]
        arrayname = "sp{0}".format(num_points)
        parts += self.generate_partition(arrayname,
                                         expr_ir["V"],
                                         expr_ir["piecewise"],
                                         expr_ir["table_ranges"],
                                         num_points)
        if parts:
            parts = [Comment("Section for piecewise constant computations")] + parts + [""]
        return parts

    def generate_varying_partition(self, num_points):
        parts = []
        expr_ir = self.ir["uflacs"]["expr_ir"][num_points]
        arrayname = "sv{0}".format(num_points)
        parts += self.generate_partition(arrayname,
                                         expr_ir["V"],
                                         expr_ir["varying"],
                                         expr_ir["table_ranges"],
                                         num_points)

        if parts:
            parts = [Comment("Section for geometrically varying computations")] + parts + [""]
        return parts

    def generate_argument_partition(self, num_points, iarg, dofrange):
        """Generate code for the partition corresponding to arguments 0..iarg within given dofblock."""
        parts = []
        # TODO: What do we want to do here? Define!
        # Should this be a single loop over i0, i1 separately outside of the double loop over (i0,i1)?
        return parts

    def generate_integrand_accumulation(self, num_points, dofblock):
        parts = []

        expr_ir = self.ir["uflacs"]["expr_ir"][num_points]
        AF = expr_ir["argument_factorization"]
        V = expr_ir["V"]
        MATR = expr_ir["modified_argument_table_ranges"]
        MA = expr_ir["modified_arguments"]

        idofs = [self.backend_access.argument_loop_index(i) for i in range(self.ir["rank"])]

        # Find the blocks to build: (TODO: This is rather awkward, having to rediscover these relations here)
        arguments_and_factors = sorted(iteritems(expr_ir["argument_factorization"]), key=lambda x: x[0])
        for args, factor_index in arguments_and_factors:
            if not all(tuple(dofblock[iarg]) == tuple(MATR[ma][1:3]) for iarg, ma in enumerate(args)):
                continue

            factors = []

            # Get factor expression
            fcode = self.expr_formatter.visit(V[factor_index])
            #fcode = self.ast_variables[factor_index] # TODO
            if fcode not in ("1", "1.0"):  # TODO: Nicer way to do this
                factors += [fcode]

            # Get table names
            argfactors = []
            for i, ma in enumerate(args):
                access = self.backend_access(MA[ma].terminal, MA[ma], MATR[ma], num_points)
                argfactors += [access]
            factors.extend(argfactors)

            # Format index access to A
            if idofs:
                dofdims = zip(idofs, self._A_shape)
                im = IndexMapping(dict((idof, idim) for idof, idim in dofdims))
                am = AxisMapping(im, [idofs])
                A_ii, = am.format_index_expressions()  # TODO: Integrate this with other code utils
            else:
                A_ii = 0
            A_access = self.backend_access.element_tensor_entry(A_ii)

            # Emit assignment
            parts += [AssignAdd(A_access, Product(factors))]

        return parts

    def generate_finishing_statements(self):
        """Generate finishing statements.

        This includes assigning to output array if there is no integration.
        """
        parts = []

        if not self.ir["quadrature_rules"]:  # Rather check ir["integral_type"]?
            # TODO: Implement for expression support
            error("Expression generation not implemented yet.")
            # TODO: If no integration, assuming we generate an expression, and assign results here
            # Corresponding code from compiler.py:
            # assign_to_variables = tfmt.output_variable_names(len(final_variable_names))
            # parts += list(format_assignments(zip(assign_to_variables, final_variable_names)))

        return parts
