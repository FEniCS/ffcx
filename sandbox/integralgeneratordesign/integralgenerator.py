
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
  tfmt.define_output_variables_reset()
  tfmt.define_piecewise_geometry()
  tfmt.define_piecewise_coefficients()
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

    def generate(self):
        """Generate entire tabulate_tensor body.

        Assumes that the code returned from here will be wrapped in a context
        that matches a suitable version of the UFC tabulate_tensor signatures.
        """
        parts = []
        parts += [self.generate_quadrature_tables()]
        parts += [self.generate_element_tables()]
        parts += [self.generate_tensor_reset()]
        parts += [self.generate_pre_quadrature_loops()]
        parts += [self.generate_quadrature_loops()]
        parts += [self.generate_post_quadrature_loops()]
        return format_code_structure(Indented(parts))

    def generate_quadrature_tables(self):
        "Generate static tables of quadrature points and weights."
        parts = ["// Quadrature weights and points"]

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

    def generate_element_tables(self): # FIXME
        "Generate static tables with precomputed element basis function values in quadrature points."
        parts = ["// FIXME: Precomputed element basis function values"]
        return parts

    def generate_tensor_reset(self):
        "Generate statements for resetting the element tensor to zero."
        memset = "memset(%s, 0, %d * sizeof(%s[0]));" % (names.A, self._A_size, names.A)
        parts = [langfmt.comment("Reset element tensor"), memset, ""]
        return parts

    def generate_pre_quadrature_loops(self): # FIXME
        """Generate statements prior to the quadrature loop.

        This mostly includes computations involving piecewise constant geometry and coefficients.
        """
        parts = []
        parts += ["// Piecewise constant stage"]
        parts += ["// FIXME: Implement this"]
        return parts

    def generate_quadrature_loops(self):
        "Generate all quadrature loops."
        parts = [self.generate_quadrature_loop(num_points) for num_points in self._num_points]
        return parts

    def generate_post_quadrature_loops(self): # TODO
        """Generate finishing statements.

        This includes assigning to output array if there is no integration.
        """
        parts = []

        # TODO: If no integration, assuming we generate an expression, and assign results here
        # Corresponding code from compiler.py:
        #assign_to_variables = tfmt.output_variable_names(len(final_variable_names))
        #parts += list(format_assignments(zip(assign_to_variables, final_variable_names)))

        return parts

    def generate_quadrature_loop(self, num_points):
        "Generate a single quadrature loop."
        return ForRange(names.iq, 0, num_points,
                        body=self.generate_quadrature_body(num_points))

    def generate_quadrature_body_setup(self, num_points): # FIXME
        """
        """
        parts = ["// Quadrature loop body setup {0}".format(num_points)]
        # FIXME: All leftover argument independent computations here
        return parts

    def generate_quadrature_body(self, num_points):
        """
        """
        parts = []
        parts += [self.generate_quadrature_body_setup(num_points)]

        # TODO: Get this from ir
        self._integrand_term_factors = {}
        self._integrand_term_factors[num_points] = {} # IM in factorization code
        self._modified_argument_dofrange = {}
        self._modified_argument_dofrange[num_points] = [] # TODO: Build from dof ranges of AV in factorization code
        self._modified_arguments = {}
        self._modified_arguments[num_points] = [] # AV in factorization code

        # ITF: tuple(modified_argument_indices) -> code_index
        ITF = self._integrand_term_factors[num_points]
        # MAD: modified_argument_index -> dofrange
        MAD = self._modified_argument_dofrange[num_points]

        # Generate pre-argument loop computations
        for mas in ITF:
            db = tuple(MAD[ma] for ma in mas)
            ssa_index = ITF[mas]
            # TODO: Generate code for f*D and store reference to it with mas
            #f = self._ssa[ssa_index]
            #vname = name_from(mas)
            #vcode = code_from(ssa_index)
            #code += ["%s = (%s) * D;" % (vname, vcode)]
            #factors[mas] = vname

        if self._num_arguments == 0: # Functional
            # TODO: Partitions data structure may look something like this:
            #itg_partition = self._partitions["integrand"][num_points][()]

            # Accumulate into element scalar
            parts += [self.generate_integrand_accumulation(num_points, ())]

        elif self._num_arguments == 1: # Linear form
            # Loop over dofranges of argument 0 and accumulate into element vector
            drs0 = sorted(MAD[mas[0]] for mas in ITF)
            for dofrange0 in drs0:
                (b0,e0) = dofrange0
                dofblock = (dofrange0,)

                # TODO: Partitions data structure may look something like this:
                #arg_partition = self._partitions["argument"][num_points][0][dofrange0]
                #itg_partition = self._partitions["integrand"][num_points][dofblock]

                body0 = [self.generate_argument_partition(num_points, 0, dofblock)]
                body0 += [self.generate_integrand_accumulation(num_points, dofblock)]
                loop0 = ForRange(self._idofs[0], b0, e0, body=body0)

                parts += [loop0]

        elif self._num_arguments == 2: # Bilinear form
            # Loop over dofranges of argument 0 (rows of element matrix)
            drs0 = sorted(MAD[mas[0]] for mas in ITF)
            for dofrange0 in drs0:
                (b0,e0) = dofrange0
                dofblock = (dofrange0,)

                # Find dofblocks starting with this dofrange
                dofblocks = sorted(tuple(MAD[ma] for ma in mas) for mas in ITF if MAD[mas[0]] == dofrange0)

                # TODO: Partitions data structure may look something like this:
                #arg_partition = self._partitions["argument"][num_points][0][dofrange0]

                body0 = [self.generate_argument_partition(num_points, 0, dofblock)]

                # Loop over dofranges of argument 1 and accumulate into element matrix
                drs1 = sorted(db[1] for db in dofblocks)
                for dofrange1 in drs1:
                    (b1,e1) = dofrange1
                    dofblock = (dofrange0,dofrange1)

                    # TODO: Partitions data structure may look something like this:
                    #arg_partition = self._partitions["argument"][num_points][1][dofrange1]
                    #itg_partition = self._partitions["integrand"][num_points][dofblock]

                    body1 = [self.generate_argument_partition(num_points, 1, dofblock)]
                    body1 += [self.generate_integrand_accumulation(num_points, dofblock)]
                    loop1 = ForRange(self._idofs[1], b1, e1, body=body1)

                    body0 += [loop1]

                loop0 = ForRange(self._idofs[0], b0, e0, body=body0)

                parts += [loop0]
        else:
            error("Bonus points: Generalize this recursively to support num_arguments > 2.")

        return parts

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
        partitions["argument"] = dict((np, [dict() for i in range(rank)]) for np in num_points)

        # partitions["integrand"][np][dofblock] = partition for computing final integrand contribution in this dofblock
        partitions["integrand"] = dict((np, dict()) for np in num_points)

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
        self._dofblocks = dict( (num_points, db.keys()) for num_points, db in self._dofblock_integrand_partition.iteritems())

    def generate_argument_partition(self, num_points, iarg, dofblock): # FIXME
        parts = []

        # FIXME: Get partition associated with (num_points, iarg, dofblock)
        #p = self._dofblock_partition[num_points][iarg].get(dofblock)
        #p = self._partitions["argument"][num_points][iarg].get(dofrange)
        p = None
        if p is not None:
            # FIXME: Generate partition computation here
            #parts += generate_partition_assignments(p)
            pass
        # TODO: Remove this mock code
        parts += ["s[...] = ...; // {0} x {1}".format(iarg, dofblock)]

        return parts

    def generate_integrand_accumulation(self, num_points, dofblock): # FIXME
        parts = []

        # FIXME: Get partition associated with (num_points, dofblock)
        #p = self._dofblock_integrand_partition[num_points][dofblock]
        #p = self._partitions["integrand"][num_points].get(dofblock)

        # FIXME: Generate accumulation properly
        #parts += generate_partition_accumulations(p)
        parts += ["A[{0}] += f * v * D;".format(dofblock)] # TODO: Remove this mock code

        # Corresponding code from compiler.py:
        #final_variable_names = [format_register_variable(p, r) for r in ir["target_registers"]]
        #assign_to_variables = tfmt.output_variable_names(len(final_variable_names))
        #scaling_factor = tfmt.accumulation_scaling_factor()
        #parts += list(format_scaled_additions(zip(assign_to_variables,
        #                                          final_variable_names),
        #                                          scaling_factor))

        return parts


class CellIntegralGenerator(IntegralGenerator):
    pass

class ExteriorFacetIntegralGenerator(IntegralGenerator):
    pass

class InteriorFacetIntegralGenerator(IntegralGenerator):
    pass

class DiracIntegralGenerator(IntegralGenerator):
    pass

class SubcellIntegralGenerator(IntegralGenerator):
    pass

class CutcellIntegralGenerator(IntegralGenerator):
    pass
