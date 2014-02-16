
from ufl.common import product

from uflacs.utils.log import debug, info, warning, error, uflacs_assert
from uflacs.codeutils.format_code_structure import format_code_structure, Indented, Block, ForRange
from uflacs.geometry.default_names import names

# TODO: Refactor
from uflacs.backends.ffc.ffc_statement_formatter import CppStatementFormatterRules
langfmt = CppStatementFormatterRules()

def create_fake_ir():
    # TODO: This is just fake test data
    ir = {}
    ir["prim_idims"] = [7, 8]
    ir["quadrature_weights"] = {
        1: ( [0.5, 0.5],  [(0.1, 0.2), (0.3, 0.4)] ),
        2: ( [3.5, 3.5],  [(1.1, 1.2), (2.3, 2.4)] ),
        }
    ir["dof_ranges"] = {
        1: [
            ((0,3), (3,6)),
            ((1,2), (4,5))
            ],
        2: [
            ((1,2), (4,6)),
            ((0,5),)
            ],
        }
    return ir

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
    def __init__(self):
        # TODO: This is just fake test data
        ir = create_fake_ir()

        # Parse quadrature rules
        quadrature_rules = ir["quadrature_weights"]
        self._num_points = []
        self._weights = {}
        self._points = {}
        for num_points in sorted(quadrature_rules.keys()):
            self._num_points.append(num_points)
            w, p = quadrature_rules[num_points]
            self._weights[num_points] = w
            self._points[num_points] = p
        if len(quadrature_rules) > 1:
            warning("Multiple quadrature rules not fully implemented, will likely crash somewhere.")

        # Parse argument space dimensions etc.
        self._dof_ranges = ir["dof_ranges"]
        self._argument_space_dimensions = ir["prim_idims"] # FIXME: *2 for dS?
        self._A_size = product(self._argument_space_dimensions)
        self._num_arguments = len(self._argument_space_dimensions)
        self._idofs = ["%s%d" % (names.ia, i) for i in range(self._num_arguments)]

        # TODO: Get this from ir
        # integrand_term_factors: tuple(modified_argument_indices) -> code_index
        self._integrand_term_factors = {}
        self._integrand_term_factors[num_points] = {} # IM in factorization code
        # modified_argument_dofrange: modified_argument_index -> dofrange
        self._modified_argument_dofrange = {}
        self._modified_argument_dofrange[num_points] = [] # TODO: Build from dof ranges of AV in factorization code
        #self._modified_arguments = {}
        #self._modified_arguments[num_points] = [] # AV in factorization code

    def build_datastructures(self):
        # FIXME: Build this representation of integrands, either by looping
        #        like this or by looping over data structures in a more
        #        unstructured fashion. Either way this is the structure we want:

        num_points = [3]

        #dofrange = (begin, end)
        #dofblock = ()  |  (dofrange0,)  |  (dofrange0, dofrange1)

        partitions = {}

        # partitions["piecewise"] = partition of expressions independent of quadrature and argument loops
        partitions["piecewise"] = []

        # partitions["varying"][np] = partition of expressions dependent on np quadrature but independent of argument loops
        partitions["varying"] = dict((np, []) for np in num_points)

        # partitions["argument"][np][iarg][dofrange] = partition depending on this dofrange of argument iarg
        partitions["argument"] = dict((np, [{} for i in range(rank)]) for np in num_points)

        # partitions["integrand"][np][dofblock] = partition for computing final integrand contribution in this dofblock
        partitions["integrand"] = dict((np, {}) for np in num_points)

        #dofblock_partition[np][iarg][dofrange] == partitions["argument"][np][iarg][dofrange]
        #dofblock_integrand_partition[np][dofrange] == partitions["integrand"][np][dofblock]

        argument_dofblocks = {} # { num_points: [dofblock, ...] } # FIXME: Fill this first
        dofblock_partition = {} # { num_points: { iarg: { dofblock: partition } } }
        for num_points in sorted(argument_dofblocks.keys()):
            dofblock_partition[num_points] = {}
            for iarg in range(rank):
                dofblock_partition[num_points][iarg] = {}
                for dofblock in argument_dofblocks[num_points]: # FIXME
                    # FIXME: Make this the partition of code that depends on only argument iarg
                    #        for this dofblock in the integrand corresponding to num_points
                    partition = []
                    dofblock_partition[num_points][iarg][dofblock] = partition

        dofblocks = {} # { num_points: [dofblock, ...] } # FIXME: Fill this first
        dofblock_integrand_partition = {} # { num_points: { dofblock: partition } }
        for num_points in sorted(dofblocks.keys()):
            dofblock_integrand_partition[num_points] = {}
            for dofblock in dofblocks[num_points]:
                # FIXME: Make this the partition of code that depends on this full dofblock
                #        in the integrand corresponding to num_points
                partition = []
                dofblock_integrand_partition[num_points][dofblock] = partition

        self._dofblock_partition = dofblock_partition
        self._dofblock_integrand_partition = dofblock_integrand_partition
        self._dofblocks = dofblocks
        self._dofblocks = dict( (num_points, sorted(db.keys()))
                                for num_points, db in self._dofblock_integrand_partition.iteritems())

    def generate(self):
        """Generate entire tabulate_tensor body.

        Assumes that the code returned from here will be wrapped in a context
        that matches a suitable version of the UFC tabulate_tensor signatures.
        """
        parts = []
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

        if self._num_points:
            parts += ["// Section for quadrature weights and points"]

        for num_points in self._num_points:
            weights = self._weights[num_points]
            points = self._points[num_points]
            pdim = len(points[0])

            weights = "{ %s }" % langfmt.precision_floats(weights)
            points = "{ %s }" % langfmt.precision_floats(x for p in points for x in p)

            #wname = "%s%d" % (names.weights, num_points) # TODO: Use this everywhere
            #pname = "%s%d" % (names.points, num_points) # TODO: Use this everywhere
            wname = "%s" % names.weights
            pname = "%s" % names.points

            parts += [langfmt.array_decl("static const double", wname, num_points, weights)]
            parts += [langfmt.array_decl("static const double", pname, num_points*pdim, points)]
            parts += [""]

        return parts

    def generate_element_tables(self): # FIXME: Generate element tables here
        "Generate static tables with precomputed element basis function values in quadrature points."
        parts = []
        parts += ["// FIXME: Section for precomputed element basis function values"]
        return parts

    def generate_tensor_reset(self):
        "Generate statements for resetting the element tensor to zero."
        memset = "memset(%s, 0, %d * sizeof(%s[0]));" % (names.A, self._A_size, names.A)
        parts = [langfmt.comment("Reset element tensor"), memset, ""]
        return parts

    def generate_quadrature_loops(self):
        "Generate all quadrature loops."
        parts = []
        for num_points in self._num_points:
            body = self.generate_quadrature_body(num_points)
            parts += [ForRange(names.iq, 0, num_points, body=body)]
        return parts

    def generate_quadrature_body(self, num_points):
        """
        """
        parts = []



        # FIXME: Get this from ir in constructor, remove this hack
        # integrand_term_factors: tuple(modified_argument_indices) -> code_index
        self._integrand_term_factors = {}
        self._integrand_term_factors[num_points] = {} # IM in factorization code
        # modified_argument_dofrange: modified_argument_index -> dofrange
        self._modified_argument_dofrange = {}
        self._modified_argument_dofrange[num_points] = []# FIXME: Build from dof ranges of AV in factorization code
        #self._modified_arguments = {}
        #self._modified_arguments[num_points] = [] # AV in factorization code


        parts += ["// Quadrature loop body setup {0}".format(num_points)]

        parts += self.generate_varying_partition(num_points)

        # Generate pre-argument loop computations # FIXME: What is this?
        for mas in self._integrand_term_factors[num_points]:
            dofblock = tuple(self._modified_argument_dofrange[num_points][ma] for ma in mas)
            ssa_index = self._integrand_term_factors[num_points][mas]

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
        if iarg == self._num_arguments:
            # At the innermost argument loop level we accumulate into the element tensor
            parts += [self.generate_integrand_accumulation(num_points, outer_dofblock)]
            return parts
        assert iarg < self._num_arguments

        # integrand_term_factors: tuple(modified_argument_indices) -> code_index
        # modified_argument_dofrange: modified_argument_index -> dofrange
        modified_argument_dofrange = self._modified_argument_dofrange[num_points]

        # Find dofranges at this loop level iarg starting with outer_dofblock
        dofranges = sorted(modified_argument_dofrange[mas[iarg]]
                           for mas in self._integrand_term_factors[num_points]
                           if all(modified_argument_dofrange[mas[i]] == outer_dofblock[i]
                                  for i in xrange(iarg)))

        # Build loops for each dofrange
        for dofrange in dofranges:
            dofblock = outer_dofblock + (dofrange,)
            body = []
            # Generate code partition for dofblock at this loop level
            body += self.generate_argument_partition(num_points, iarg, dofblock)

            # Generate nested inner loops (only triggers for forms with two or more arguments
            body += self.generate_quadrature_body_dofblocks(num_points, dofblock)

            # Wrap setup, subloops, and accumulation in a loop for this level
            parts += [ForRange(self._idofs[level], dofrange[0], dofrange[1], body=body)]
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

    def generate_argument_partition(self, num_points, iarg, dofblock): # FIXME: Define better and implement
        """Generate code for the partition corresponding to arguments 0..iarg within given dofblock."""
        parts = []

        # FIXME: Get partition associated with (num_points, iarg, dofblock)
        p = None
        #p = self._dofblock_partition[num_points][iarg].get(dofblock)
        #p = self._partitions["argument"][num_points][iarg].get(dofblock[iarg])

        if p is not None:
            parts += generate_partition_assignments(p) # FIXME: Generate partition computation here
        else:
            parts += ["// FIXME: sa[...] = ...; // {0} x {1}".format(iarg, dofblock)] # TODO: Remove this mock code

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

        if not self._num_points:
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
