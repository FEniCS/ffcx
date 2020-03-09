# Copyright (C) 2013-2017 Martin Sandve Alnæs
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Controlling algorithm for building the tabulate_tensor source structure from factorized representation."""

import collections
import logging

import ufl
from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.C.format_lines import format_indented_lines
from ffcx.ir.representationutils import initialize_integral_code
from ffcx.ir.uflacs.elementtables import piecewise_ttypes

logger = logging.getLogger(__name__)


def generate_integral_code(ir, parameters):
    """Generate code for integral from intermediate representation."""

    logger.info("Generating code from ffcx.ir.uflacs representation")

    # Create FFCX C backend
    backend = FFCXBackend(ir, parameters)

    # Configure kernel generator
    ig = IntegralGenerator(ir, backend)

    # Generate code ast for the tabulate_tensor body
    parts = ig.generate()

    # Format code as string
    body = format_indented_lines(parts.cs_format(ir.precision), 1)

    # Generate generic ffcx code snippets and add uflacs specific parts
    code = initialize_integral_code(ir, parameters)
    code["tabulate_tensor"] = body

    return code


class IntegralGenerator(object):
    def __init__(self, ir, backend):
        # Store ir
        self.ir = ir

        # Backend specific plugin with attributes
        # - language: for translating ufl operators to target language
        # - symbols: for translating ufl operators to target language
        # - definitions: for defining backend specific variables
        # - access: for accessing backend specific variables
        self.backend = backend

        # Set of operator names code has been generated for,
        # used in the end for selecting necessary includes
        self._ufl_names = set()

        # Initialize lookup tables for variable scopes
        self.init_scopes()

        # Cache
        self.shared_symbols = {}

        # Set of counters used for assigning names to intermediate variables
        self.symbol_counters = collections.defaultdict(int)

        # Contains the variable names to say whether or not each dof in a space needs its direction
        # to be reversed (for vector dofs)
        # If a space's id is not in this, then no reversals are needed
        self.dof_reflections = {}

        # Contains the variable names of the dofmaps for the spaces whose dofs need to be
        # reversed
        # If a space's id is not in this, then no reversals are needed or the dofmap is trivial
        self.table_dofmaps = {}

    def init_scopes(self):
        """Initialize variable scope dicts."""
        # Reset variables, separate sets for quadrature loop
        self.scopes = {num_points: {} for num_points in self.ir.all_num_points}
        self.scopes[None] = {}

    def set_var(self, num_points, v, vaccess):
        """Set a new variable in variable scope dicts.

        Scope is determined by num_points which identifies the
        quadrature loop scope or None if outside quadrature loops.

        v is the ufl expression and vaccess is the CNodes
        expression to access the value in the code.

        """
        self.scopes[num_points][v] = vaccess

    def get_var(self, num_points, v):
        """Lookup ufl expression v in variable scope dicts.

        Scope is determined by num_points which identifies the
        quadrature loop scope or None if outside quadrature loops.

        If v is not found in quadrature loop scope, the piecewise
        scope (None) is checked.

        Returns the CNodes expression to access the value in the code.
        """
        if v._ufl_is_literal_:
            return self.backend.ufl_to_language.get(v)
        f = self.scopes[num_points].get(v)
        if f is None:
            f = self.scopes[None][v]
        return f

    def new_temp_symbol(self, basename):
        """Create a new code symbol named basename + running counter."""
        L = self.backend.language
        name = "%s%d" % (basename, self.symbol_counters[basename])
        self.symbol_counters[basename] += 1
        return L.Symbol(name)

    def get_temp_symbol(self, tempname, key):
        key = (tempname, ) + key
        s = self.shared_symbols.get(key)
        defined = s is not None
        if not defined:
            s = self.new_temp_symbol(tempname)
            self.shared_symbols[key] = s
        return s, defined

    def generate(self):
        """Generate entire tabulate_tensor body.

        Assumes that the code returned from here will be wrapped in a context
        that matches a suitable version of the UFC tabulate_tensor signatures.
        """
        L = self.backend.language

        # Assert that scopes are empty: expecting this to be called only once
        assert not any(d for d in self.scopes.values())

        parts = []

        alignment = self.ir.params['assume_aligned']
        if alignment is not None:
            parts += [L.VerbatimStatement("A = (ufc_scalar_t*)__builtin_assume_aligned(A, {});"
                                          .format(alignment)),
                      L.VerbatimStatement("w = (const ufc_scalar_t*)__builtin_assume_aligned(w, {});"
                                          .format(alignment)),
                      L.VerbatimStatement("c = (const ufc_scalar_t*)__builtin_assume_aligned(c, {});"
                                          .format(alignment)),
                      L.VerbatimStatement(
                          "coordinate_dofs = (const double*)__builtin_assume_aligned(coordinate_dofs, {});"
                          .format(alignment))]

        # Generate array of bools to say whether or not each dof needs to be reversed
        # (for vector valued basis functions)
        parts += self.generate_dof_reflections()

        # Generate dofmaps for spaces whose dofs need to be reversed
        parts += self.generate_table_dofmaps()

        # Generate the tables of quadrature points and weights
        parts += self.generate_quadrature_tables()

        # Generate the tables of basis function values and preintegrated blocks
        parts += self.generate_element_tables()

        # Generate code to compute piecewise constant scalar factors
        parts += self.generate_unstructured_piecewise_partition()

        # Loop generation code will produce parts to go before quadloops,
        # to define the quadloops, and to go after the quadloops
        all_preparts = []
        all_quadparts = []

        # Go through each relevant quadrature loop
        if self.ir.integral_type in ufl.custom_integral_types:
            preparts, quadparts = \
                self.generate_runtime_quadrature_loop()
            all_preparts += preparts
            all_quadparts += quadparts
        else:
            for num_points in self.ir.all_num_points:
                # Generate code to integrate reusable blocks of final element tensor
                preparts, quadparts = self.generate_quadrature_loop(num_points)
                all_preparts += preparts
                all_quadparts += quadparts

        # Generate code to finish computing reusable blocks outside quadloop
        preparts, quadparts = \
            self.generate_dofblock_partition(None)
        all_preparts += preparts
        all_quadparts += quadparts

        # Collect parts before, during, and after quadrature loops
        parts += all_preparts
        parts += all_quadparts

        return L.StatementList(parts)

    def generate_dof_reflections(self):
        """Generate arrays of bool saying whether each dof needs to be reflected."""
        L = self.backend.language
        c_false = L.LiteralBool(False)

        parts = []
        for element, id in self.ir.element_ids.items():
            reflect_dofs = []
            contains_reflections = False
            for dre in self.ir.element_dof_reflection_entities[element]:
                if dre is None:
                    # Dof does not need reflecting, so put false in array
                    reflect_dofs.append(c_false)
                else:
                    # Loop through entities that the direction of the dof depends on to
                    # make a conditional
                    ref = c_false
                    for j in dre:
                        if ref == c_false:
                            # No condition has been added yet, so overwrite false
                            ref = self.backend.symbols.entity_reflection(L, j)
                        else:
                            # This is not the first condition, so XOR
                            ref = L.Conditional(self.backend.symbols.entity_reflection(L, j), L.Not(ref), ref)
                    reflect_dofs.append(ref)
                    if ref != c_false:
                        # Mark this space as needing reflections
                        contains_reflections = True

            # If no dofs need reflecting, don't write any array
            if contains_reflections:
                self.dof_reflections[id] = L.Symbol("ref_dof" + str(id))
                parts.append(L.ArrayDecl(
                    "const bool", self.dof_reflections[id], (len(reflect_dofs), ), values=reflect_dofs))
        return parts

    def generate_table_dofmaps(self):
        """Generate dofmaps for spaces whose dofs need to be reflected."""
        L = self.backend.language

        parts = []
        for tname, dofmap in self.ir.table_dofmaps.items():
            id = self.ir.element_ids[self.ir.table_origins[tname][0]]
            if id in self.dof_reflections:
                # Write the dofmap as it will be needed
                for i, j in enumerate(dofmap):
                    if i != j:
                        # If a dof has been removed, write the data
                        self.table_dofmaps[tname] = L.Symbol(tname + "_dofmap")
                        parts.append(L.ArrayDecl(
                            "const int", self.table_dofmaps[tname], (len(dofmap), ), values=dofmap))
                        break
        return parts

    def generate_quadrature_tables(self):
        """Generate static tables of quadrature points and weights."""
        L = self.backend.language

        parts = []

        # No quadrature tables for custom (given argument)
        # or point (evaluation in single vertex)
        skip = ufl.custom_integral_types + ufl.measure.point_integral_types
        if self.ir.integral_type in skip:
            return parts

        alignas = self.ir.params["alignas"]
        padlen = self.ir.params["padlen"]

        # Loop over quadrature rules
        for num_points in self.ir.all_num_points:
            varying_ir = self.ir.varying_irs[num_points]

            points, weights = self.ir.quadrature_rules[num_points]
            assert num_points == len(weights)
            assert num_points == points.shape[0]

            # Generate quadrature weights array
            if varying_ir["need_weights"]:
                wsym = self.backend.symbols.weights_table(num_points)
                parts += [
                    L.ArrayDecl(
                        "static const double", wsym, num_points, weights, alignas=alignas, padlen=padlen)
                ]

            # Generate quadrature points array
            N = ufl.product(points.shape)
            if varying_ir["need_points"] and N:
                # Flatten array: (TODO: avoid flattening here, it makes padding harder)
                flattened_points = points.reshape(N)
                psym = self.backend.symbols.points_table(num_points)
                parts += [
                    L.ArrayDecl(
                        "static const double", psym, N, flattened_points, alignas=alignas)
                ]

        # Add leading comment if there are any tables
        parts = L.commented_code_list(parts, "Quadrature rules")
        return parts

    def generate_element_tables(self):
        """Generate static tables with precomputed element basis function values in quadrature points."""
        L = self.backend.language
        parts = []

        tables = self.ir.unique_tables
        table_types = self.ir.unique_table_types

        alignas = self.ir.params["alignas"]
        padlen = self.ir.params["padlen"]

        if self.ir.integral_type in ufl.custom_integral_types:
            # Define only piecewise tables
            table_names = [name for name in sorted(tables) if table_types[name] in piecewise_ttypes]
        else:
            # Define all tables
            table_names = sorted(tables)

        for name in table_names:
            table = tables[name]

            decl = L.ArrayDecl(
                "static const double", name, table.shape, table, alignas=alignas, padlen=padlen)
            parts += [decl]

        # Add leading comment if there are any tables
        parts = L.commented_code_list(parts, [
            "Precomputed values of basis functions and precomputations",
            "FE* dimensions: [permutation][entities][points][dofs]",
        ])
        return parts

    def generate_quadrature_loop(self, num_points):
        """Generate quadrature loop with for this num_points."""
        L = self.backend.language

        # Generate unstructured varying partition
        body = self.generate_unstructured_varying_partition(num_points)
        body = L.commented_code_list(
            body, "Quadrature loop body setup (num_points={0})".format(num_points))

        # Generate dofblock parts, some of this
        # will be placed before or after quadloop
        preparts, quadparts = \
            self.generate_dofblock_partition(num_points)
        body += quadparts

        # Wrap body in loop or scope
        if not body:
            # Could happen for integral with everything zero and optimized away
            quadparts = []
        elif num_points == 1:
            # For now wrapping body in Scope to avoid thinking about scoping issues
            quadparts = L.commented_code_list(L.Scope(body), "Only 1 quadrature point, no loop")
        else:
            # Regular case: define quadrature loop
            if num_points == 1:
                iq = 0
            else:
                iq = self.backend.symbols.quadrature_loop_index()
            quadparts = [L.ForRange(iq, 0, num_points, body=body)]

        return preparts, quadparts

    def generate_runtime_quadrature_loop(self):
        """Generate quadrature loop for custom integrals, with physical points given runtime."""
        L = self.backend.language

        assert self.ir.integral_type in ufl.custom_integral_types

        num_points = self.ir.fake_num_points
        chunk_size = self.ir.params["chunk_size"]

        gdim = self.ir.geometric_dimension

        alignas = self.ir.params["alignas"]

        tables = self.ir.unique_tables
        table_types = self.ir.unique_table_types

        # Generate unstructured varying partition
        body = self.generate_unstructured_varying_partition(num_points)
        body = L.commented_code_list(body, [
            "Run-time quadrature loop body setup",
            "(chunk_size={0}, analysis_num_points={1})".format(chunk_size, num_points)
        ])

        # Generate dofblock parts, some of this
        # will be placed before or after quadloop
        preparts, quadparts = \
            self.generate_dofblock_partition(num_points)
        body += quadparts

        # Wrap body in loop
        if not body:
            # Could happen for integral with everything zero and optimized away
            quadparts = []
        else:
            rule_parts = []

            # Define two-level quadrature loop; over chunks then over points in chunk
            iq_chunk = L.Symbol("iq_chunk")
            np = self.backend.symbols.num_custom_quadrature_points()
            num_point_blocks = (np + chunk_size - 1) / chunk_size
            iq = self.backend.symbols.quadrature_loop_index()

            # Not assuming runtime size to be multiple by chunk size
            num_points_in_block = L.Symbol("num_points_in_chunk")
            decl = L.VariableDecl("const int", num_points_in_block,
                                  L.Call("min", (chunk_size, np - iq_chunk * chunk_size)))
            rule_parts.append(decl)

            iq_body = L.ForRange(iq, 0, num_points_in_block, body=body)

            # Preparations for quadrature rules
            #
            varying_ir = self.ir.varying_irs[num_points]

            # Copy quadrature weights for this chunk
            if varying_ir["need_weights"]:
                cwsym = self.backend.symbols.custom_quadrature_weights()
                wsym = self.backend.symbols.custom_weights_table()
                rule_parts += [
                    L.ArrayDecl("ufc_scalar_t", wsym, chunk_size, 0, alignas=alignas),
                    L.ForRange(
                        iq,
                        0,
                        num_points_in_block,
                        body=L.Assign(wsym[iq], cwsym[chunk_size * iq_chunk + iq])),
                ]

            # Copy quadrature points for this chunk
            if varying_ir["need_points"]:
                cpsym = self.backend.symbols.custom_quadrature_points()
                psym = self.backend.symbols.custom_points_table()
                rule_parts += [
                    L.ArrayDecl("ufc_scalar_t", psym, chunk_size * gdim, 0, alignas=alignas),
                    L.ForRange(
                        iq,
                        0,
                        num_points_in_block,
                        body=[
                            L.Assign(psym[iq * gdim + i],
                                     cpsym[chunk_size * iq_chunk * gdim + iq * gdim + i])
                            for i in range(gdim)
                        ])
                ]

            # Add leading comment if there are any tables
            rule_parts = L.commented_code_list(rule_parts, "Quadrature weights and points")

            # Preparations for element tables
            table_parts = []

            # Only declare non-piecewise tables, computed inside chunk loop
            non_piecewise_tables = [
                name for name in sorted(tables) if table_types[name] not in piecewise_ttypes
            ]
            for name in non_piecewise_tables:
                table = tables[name]
                decl = L.ArrayDecl(
                    "ufc_scalar_t", name, (1, chunk_size, table.shape[2]), 0,
                    alignas=alignas)  # padlen=padlen)
                table_parts += [decl]

            table_parts += [L.Comment("FIXME: Fill element tables here")]
            # table_origins

            # Gather all in chunk loop
            chunk_body = rule_parts + table_parts + [iq_body]
            quadparts = [L.ForRange(iq_chunk, 0, num_point_blocks, body=chunk_body)]

        return preparts, quadparts

    def generate_unstructured_piecewise_partition(self):
        L = self.backend.language

        # Get annotated graph of factorisation
        F = self.ir.piecewise_ir["factorization"]

        arraysymbol = L.Symbol("sp")
        num_points = None
        parts = self.generate_partition(arraysymbol, F, "piecewise", num_points)
        parts = L.commented_code_list(parts, "Unstructured piecewise computations")
        return parts

    def generate_unstructured_varying_partition(self, num_points):
        L = self.backend.language

        # Get annotated graph of factorisation
        F = self.ir.varying_irs[num_points]["factorization"]

        arraysymbol = L.Symbol("sv%d" % num_points)
        parts = self.generate_partition(arraysymbol, F, "varying", num_points)
        parts = L.commented_code_list(parts, "Unstructured varying computations for num_points=%d" %
                                      (num_points, ))
        return parts

    def generate_partition(self, symbol, F, mode, num_points):
        L = self.backend.language

        definitions = []
        intermediates = []

        for i, attr in F.nodes.items():
            if attr['status'] != mode:
                continue
            v = attr['expression']
            mt = attr.get('mt')

            if v._ufl_is_literal_:
                vaccess = self.backend.ufl_to_language.get(v)
            elif mt is not None:
                # All finite element based terminals have table data, as well
                # as some, but not all, of the symbolic geometric terminals
                tabledata = attr.get('tr')

                # Backend specific modified terminal translation
                vaccess = self.backend.access.get(mt.terminal, mt, tabledata, num_points)
                vdef = self.backend.definitions.get(mt.terminal, mt, tabledata, num_points, vaccess)

                # Store definitions of terminals in list
                assert isinstance(vdef, list)
                definitions.extend(vdef)
            else:
                # Get previously visited operands
                vops = [self.get_var(num_points, op) for op in v.ufl_operands]

                # get parent operand
                pid = F.in_edges[i][0] if F.in_edges[i] else -1
                if pid and pid > i:
                    parent_exp = F.nodes.get(pid)['expression']
                else:
                    parent_exp = None

                # Mapping UFL operator to target language
                self._ufl_names.add(v._ufl_handler_name_)
                vexpr = self.backend.ufl_to_language.get(v, *vops)

                # Create a new intermediate for each subexpression
                # except boolean conditions and its childs
                if isinstance(parent_exp, ufl.classes.Condition):
                    # Skip intermediates for 'x' and 'y' in x<y
                    # Avoid the creation of complex valued intermediates
                    vaccess = vexpr
                elif isinstance(v, ufl.classes.Condition):
                    # Inline the conditions x < y, condition values
                    # This removes the need to handle boolean intermediate variables.
                    # With tensor-valued conditionals it may not be optimal but we
                    # let the compiler take responsibility for optimizing those cases.
                    vaccess = vexpr
                elif any(op._ufl_is_literal_ for op in v.ufl_operands):
                    # Skip intermediates for e.g. -2.0*x,
                    # resulting in lines like z = y + -2.0*x
                    vaccess = vexpr
                else:
                    # Record assignment of vexpr to intermediate variable
                    j = len(intermediates)
                    if self.ir.params["use_symbol_array"]:
                        vaccess = symbol[j]
                        intermediates.append(L.Assign(vaccess, vexpr))
                    else:
                        vaccess = L.Symbol("%s_%d" % (symbol.name, j))
                        intermediates.append(L.VariableDecl("const ufc_scalar_t", vaccess, vexpr))

            # Store access node for future reference
            self.set_var(num_points, v, vaccess)

        # Join terminal computation, array of intermediate expressions,
        # and intermediate computations
        parts = []
        if definitions:
            parts += definitions
        if intermediates:
            if self.ir.params["use_symbol_array"]:
                alignas = self.ir.params["alignas"]
                padlen = self.ir.params["padlen"]
                parts += [L.ArrayDecl("ufc_scalar_t", symbol, len(intermediates), alignas=alignas, padlen=padlen)]
            parts += intermediates
        return parts

    def generate_dofblock_partition(self, num_points):
        if num_points is None:  # NB! None meaning piecewise partition, not custom integral
            block_contributions = self.ir.piecewise_ir["block_contributions"]
        else:
            block_contributions = self.ir.varying_irs[num_points]["block_contributions"]

        preparts = []
        quadparts = []

        blocks = [(blockmap, blockdata)
                  for blockmap, contributions in sorted(block_contributions.items())
                  for blockdata in contributions]

        for blockmap, blockdata in blocks:

            # Define code for block depending on mode
            block_preparts, block_quadparts = \
                self.generate_block_parts(num_points, blockmap, blockdata)

            # Add definitions
            preparts.extend(block_preparts)

            # Add computations
            quadparts.extend(block_quadparts)

        return preparts, quadparts

    def get_entities(self, blockdata):
        if self.ir.integral_type == "interior_facet":
            # Get the facet entities
            entities = []
            for r in blockdata.restrictions:
                if r is None:
                    entities.append(0)
                else:
                    entities.append(self.backend.symbols.entity(self.ir.entitytype, r))
            if blockdata.transposed:
                return (entities[1], entities[0])
            else:
                return tuple(entities)
        else:
            entity = self.backend.symbols.entity(self.ir.entitytype, None)
            return (entity, )

    def get_permutations(self, blockdata):
        """Get the quadrature permutations for facets."""
        perms = []
        for r in blockdata.restrictions:
            if blockdata.is_permuted:
                qp = self.backend.symbols.quadrature_permutation(0)
                if r == "-":
                    qp = self.backend.symbols.quadrature_permutation(1)
            else:
                qp = 0
            perms.append(qp)
        if blockdata.transposed:
            return (perms[1], perms[0])
        else:
            return tuple(perms)

    def get_arg_factors(self, blockdata, block_rank, num_points, iq, indices):
        arg_factors = []
        for i in range(block_rank):
            mad = blockdata.ma_data[i]
            td = mad.tabledata
            if td.is_piecewise:
                scope = self.ir.piecewise_ir["modified_arguments"]
            else:
                scope = self.ir.varying_irs[num_points]["modified_arguments"]
            mt = scope[mad.ma_index]

            # Translate modified terminal to code
            # TODO: Move element table access out of backend?
            #       Not using self.backend.access.argument() here
            #       now because it assumes too much about indices.

            table = self.backend.symbols.element_table(td, self.ir.entitytype, mt.restriction)

            assert td.ttype != "zeros"

            if td.ttype == "ones":
                arg_factor = self.get_vector_reflection(td.name, indices)
            elif td.ttype == "quadrature":  # TODO: Revisit all quadrature ttype checks
                assert self.get_vector_reflection(td.name, iq) == 1
                arg_factor = table[iq]
            else:
                # Assuming B sparsity follows element table sparsity
                arg_factor = self.get_vector_reflection(td.name, indices) * table[indices[i]]
            arg_factors.append(arg_factor)
        return arg_factors

    def generate_block_parts(self, num_points, blockmap, blockdata):
        """Generate and return code parts for a given block.

        Returns parts occuring before, inside, and after
        the quadrature loop identified by num_points.

        Should be called with num_points=None for quadloop-independent blocks.
        """
        L = self.backend.language

        # The parts to return
        preparts = []
        quadparts = []

        block_rank = len(blockmap)
        blockdims = tuple(len(dofmap) for dofmap in blockmap)

        ttypes = blockdata.ttypes
        if "zeros" in ttypes:
            raise RuntimeError("Not expecting zero arguments to be left in dofblock generation.")

        if num_points is None:
            iq = None
        elif num_points == 1:
            iq = 0
        else:
            iq = self.backend.symbols.quadrature_loop_index()

        # Override dof index with quadrature loop index for arguments with
        # quadrature element, to index B like B[iq*num_dofs + iq]
        arg_indices = tuple(self.backend.symbols.argument_loop_index(i) for i in range(block_rank))
        B_indices = []
        for i in range(block_rank):
            if ttypes[i] == "quadrature":
                B_indices.append(iq)
            else:
                B_indices.append(arg_indices[i])
        B_indices = list(B_indices)

        # Get factor expression
        if blockdata.factor_is_piecewise:
            F = self.ir.piecewise_ir["factorization"]
        else:
            F = self.ir.varying_irs[num_points]["factorization"]

        if len(blockdata.factor_indices_comp_indices) > 1:
            raise RuntimeError("Code generation for non-scalar integrals unsupported")

        # We have scalar integrand here, take just the factor index
        factor_index = blockdata.factor_indices_comp_indices[0][0]

        v = F.nodes[factor_index]['expression']
        f = self.get_var(num_points, v)

        # Quadrature weight was removed in representation, add it back now
        if num_points is None:
            weight = L.LiteralFloat(1.0)
        elif self.ir.integral_type in ufl.custom_integral_types:
            weights = self.backend.symbols.custom_weights_table()
            weight = weights[iq]
        else:
            weights = self.backend.symbols.weights_table(num_points)
            weight = weights[iq]

        # Define fw = f * weight
        assert not blockdata.transposed, "Not handled yet"

        # Fetch code to access modified arguments
        arg_factors = self.get_arg_factors(blockdata, block_rank, num_points, iq, B_indices)

        fw_rhs = L.float_product([f, weight])
        if not isinstance(fw_rhs, L.Product):
            fw = fw_rhs
        else:
            # Define and cache scalar temp variable
            key = (num_points, factor_index, blockdata.factor_is_piecewise)
            fw, defined = self.get_temp_symbol("fw", key)
            if not defined:
                quadparts.append(L.VariableDecl("const ufc_scalar_t", fw, fw_rhs))

        # Naively accumulate integrand for this block in the innermost loop
        assert not blockdata.transposed
        A_shape = self.ir.tensor_shape

        Asym = self.backend.symbols.element_tensor()
        A = L.FlattenedArray(Asym, dims=A_shape)

        B_rhs = L.float_product([fw] + arg_factors)
        A_indices = []

        for i, bm in enumerate(blockmap):
            offset = blockmap[i][0]
            A_indices.append(arg_indices[i] + offset)

        body = L.AssignAdd(A[A_indices], B_rhs)

        for i in reversed(range(block_rank)):
            body = L.ForRange(B_indices[i], 0, blockdims[i], body=body)
        quadparts += [body]

        return preparts, quadparts

    def get_vector_reflection(self, pname, indices):
        """Get the vector reflection for entry the table pname accessed using indices."""
        L = self.backend.language
        origin = self.ir.table_origins[pname]
        if isinstance(origin[0], str):
            # If the table is preintegrated, then origin will be a tuple of strings
            # to identify which tables were used for each dof dimension
            tablenames = origin
        else:
            # Otherwise, there is only one tablename; put it in a list so we can iterate
            tablenames = [pname]
        used_indices = indices[-len(tablenames):]

        output = 1
        for tablename, index in zip(tablenames, used_indices):
            element = self.ir.table_origins[tablename][0]
            id = self.ir.element_ids[element]
            if id in self.dof_reflections:
                # If at least one vector dof needs reflecting, return a conditional that gives -1
                # if the dof needs negating
                if tablename in self.table_dofmaps:
                    index = self.table_dofmaps[tablename][index]
                output *= L.Conditional(self.dof_reflections[id][index], 1, -1)
        return output
