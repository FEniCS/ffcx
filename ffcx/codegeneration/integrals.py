# Copyright (C) 2015-2020 Martin Sandve Alnæs and Michal Habera
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import collections
import itertools
import logging

import ufl
from ffcx.codegeneration import geometry
from ffcx.codegeneration import integrals_template as ufc_integrals
from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.C.format_lines import format_indented_lines
from ffcx.codegeneration.utils import apply_transformations_to_data
from ffcx.ir.elementtables import piecewise_ttypes

logger = logging.getLogger("ffcx")


def generator(ir, parameters):
    logger.info("Generating code for integral:")
    logger.info(f"--- type: {ir.integral_type}")
    logger.info(f"--- name: {ir.name}")

    """Generate code for an integral."""
    factory_name = ir.name
    integral_type = ir.integral_type

    # Format declaration
    if integral_type == "custom":
        declaration = ufc_integrals.custom_declaration.format(factory_name=factory_name)
    else:
        declaration = ufc_integrals.declaration.format(factory_name=factory_name)

    # Create FFCx C backend
    backend = FFCXBackend(ir, parameters)

    # Configure kernel generator
    ig = IntegralGenerator(ir, backend)

    # Generate code ast for the tabulate_tensor body
    parts = ig.generate()

    # Format code as string
    body = format_indented_lines(parts.cs_format(ir.precision), 1)

    # Generate generic FFCx code snippets and add specific parts
    code = {}
    code["class_type"] = ir.integral_type + "_integral"
    code["name"] = ir.name
    code["members"] = ""
    code["constructor"] = ""
    code["constructor_arguments"] = ""
    code["initializer_list"] = ""
    code["destructor"] = ""

    L = backend.language
    if len(ir.enabled_coefficients) > 0:
        code["enabled_coefficients_init"] = L.ArrayDecl(
            "bool", f"enabled_coefficients_{ir.name}",
            values=ir.enabled_coefficients, sizes=len(ir.enabled_coefficients))
        code["enabled_coefficients"] = f"enabled_coefficients_{ir.name}"
    else:
        code["enabled_coefficients_init"] = ""
        code["enabled_coefficients"] = L.Null()

    code["additional_includes_set"] = set()  # FIXME: Get this out of code[]
    code["tabulate_tensor"] = body

    if parameters["tabulate_tensor_void"]:
        code["tabulate_tensor"] = ""

    implementation = ufc_integrals.factory.format(
        factory_name=factory_name,
        enabled_coefficients=code["enabled_coefficients"],
        enabled_coefficients_init=code["enabled_coefficients_init"],
        tabulate_tensor=code["tabulate_tensor"],
        needs_transformation_data=ir.needs_transformation_data)

    return declaration, implementation


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

    def init_scopes(self):
        """Initialize variable scope dicts."""
        # Reset variables, separate sets for each quadrature rule
        self.scopes = {quadrature_rule: {} for quadrature_rule in self.ir.integrand.keys()}
        self.scopes[None] = {}

    def set_var(self, quadrature_rule, v, vaccess):
        """Set a new variable in variable scope dicts.

        Scope is determined by quadrature_rule which identifies the
        quadrature loop scope or None if outside quadrature loops.

        v is the ufl expression and vaccess is the CNodes
        expression to access the value in the code.

        """
        self.scopes[quadrature_rule][v] = vaccess

    def get_var(self, quadrature_rule, v):
        """Lookup ufl expression v in variable scope dicts.

        Scope is determined by quadrature rule which identifies the
        quadrature loop scope or None if outside quadrature loops.

        If v is not found in quadrature loop scope, the piecewise
        scope (None) is checked.

        Returns the CNodes expression to access the value in the code.
        """
        if v._ufl_is_literal_:
            return self.backend.ufl_to_language.get(v)
        f = self.scopes[quadrature_rule].get(v)
        if f is None:
            f = self.scopes[None].get(v)
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
        if alignment != -1:
            parts += [L.VerbatimStatement(f"A = (ufc_scalar_t*)__builtin_assume_aligned(A, {alignment});"),
                      L.VerbatimStatement(f"w = (const ufc_scalar_t*)__builtin_assume_aligned(w, {alignment});"),
                      L.VerbatimStatement(f"c = (const ufc_scalar_t*)__builtin_assume_aligned(c, {alignment});"),
                      L.VerbatimStatement(
                          f"coordinate_dofs = (const double*)__builtin_assume_aligned(coordinate_dofs, {alignment});")]

        # Generate the tables of quadrature points and weights
        parts += self.generate_quadrature_tables()

        # Generate the tables of basis function values and preintegrated
        # blocks
        parts += self.generate_element_tables()

        # Generate the tables of geometry data that are needed
        parts += self.generate_geometry_tables()

        # Loop generation code will produce parts to go before
        # quadloops, to define the quadloops, and to go after the
        # quadloops
        all_preparts = []
        all_quadparts = []

        for rule in self.ir.integrand.keys():
            # Generate code to compute piecewise constant scalar factors
            all_preparts += self.generate_piecewise_partition(rule)

            # Generate code to integrate reusable blocks of final
            # element tensor
            preparts, quadparts = self.generate_quadrature_loop(rule)
            all_preparts += preparts
            all_quadparts += quadparts

        # Collect parts before, during, and after quadrature loops
        parts += all_preparts
        parts += all_quadparts

        return L.StatementList(parts)

    def generate_quadrature_tables(self):
        """Generate static tables of quadrature points and weights."""
        L = self.backend.language

        parts = []

        # No quadrature tables for custom (given argument) or point
        # (evaluation in single vertex)
        skip = ufl.custom_integral_types + ufl.measure.point_integral_types
        if self.ir.integral_type in skip:
            return parts

        padlen = self.ir.params["padlen"]

        # Loop over quadrature rules
        for quadrature_rule, integrand in self.ir.integrand.items():

            num_points = quadrature_rule.weights.shape[0]
            # Generate quadrature weights array
            wsym = self.backend.symbols.weights_table(quadrature_rule)
            parts += [
                L.ArrayDecl(
                    "static const double", wsym, num_points,
                    quadrature_rule.weights, padlen=padlen)
            ]

        # Add leading comment if there are any tables
        parts = L.commented_code_list(parts, "Quadrature rules")
        return parts

    def generate_geometry_tables(self):
        """Generate static tables of geometry data."""
        L = self.backend.language

        ufl_geometry = {
            ufl.geometry.FacetEdgeVectors: "facet_edge_vertices",
            ufl.geometry.CellFacetJacobian: "reference_facet_jacobian",
            ufl.geometry.ReferenceCellVolume: "reference_cell_volume",
            ufl.geometry.ReferenceFacetVolume: "reference_facet_volume",
            ufl.geometry.ReferenceCellEdgeVectors: "reference_edge_vectors",
            ufl.geometry.ReferenceFacetEdgeVectors: "facet_reference_edge_vectors",
            ufl.geometry.ReferenceNormal: "reference_facet_normals",
            ufl.geometry.FacetOrientation: "facet_orientation"
        }
        cells = {t: set() for t in ufl_geometry.keys()}

        for integrand in self.ir.integrand.values():
            for attr in integrand["factorization"].nodes.values():
                mt = attr.get("mt")
                if mt is not None:
                    t = type(mt.terminal)
                    if t in ufl_geometry:
                        cells[t].add(mt.terminal.ufl_domain().ufl_cell().cellname())

        parts = []
        for i, cell_list in cells.items():
            for c in cell_list:
                parts.append(geometry.write_table(L, ufl_geometry[i], c))

        return parts

    def generate_element_tables(self):
        """Generate static tables with precomputed element basis function values in quadrature points."""
        L = self.backend.language
        parts = []

        tables = self.ir.unique_tables
        table_types = self.ir.unique_table_types

        padlen = self.ir.params["padlen"]

        if self.ir.integral_type in ufl.custom_integral_types:
            # Define only piecewise tables
            table_names = [name for name in sorted(tables) if table_types[name] in piecewise_ttypes]
        else:
            # Define all tables
            table_names = sorted(tables)

        for name in table_names:
            table = tables[name]
            parts += self.declare_table(name, table, padlen)

        # Add leading comment if there are any tables
        parts = L.commented_code_list(parts, [
            "Precomputed values of basis functions and precomputations",
            "FE* dimensions: [permutation][entities][points][dofs]",
        ])
        return parts

    def declare_table(self, name, table, padlen):
        """Declare a table.

        If the dof dimensions of the table have dof rotations, apply
        these rotations.

        """
        L = self.backend.language

        if not self.ir.table_needs_transformation_data[name]:
            return [L.ArrayDecl(
                "static const double", name, table.shape, table, padlen=padlen)]

        out = [L.ArrayDecl(
            "double", name, table.shape, table, padlen=padlen)]

        dummy_vars = tuple(0 if j == 1 else L.Symbol(f"i{i}") for i, j in enumerate(table.shape[:-1]))
        ranges = tuple((dummy_vars[i], 0, j) for i, j in enumerate(table.shape[:-1]) if j != 1)
        apply_transformations = apply_transformations_to_data(
            L, self.ir.table_dof_base_transformations[name], self.ir.cell_shape, L.Symbol(name),
            indices=lambda dof: dummy_vars + (dof, ), ranges=ranges)
        if len(apply_transformations) > 0:
            out += ["{"] + apply_transformations + ["}"]
        return out

    def generate_quadrature_loop(self, quadrature_rule):
        """Generate quadrature loop with for this num_points."""
        L = self.backend.language

        # Generate varying partition
        body = self.generate_varying_partition(quadrature_rule)
        body = L.commented_code_list(
            body, f"Quadrature loop body setup for quadrature rule {quadrature_rule.id()}")

        # Generate dofblock parts, some of this will be placed before or
        # after quadloop
        preparts, quadparts = \
            self.generate_dofblock_partition(quadrature_rule)
        body += quadparts

        # Wrap body in loop or scope
        if not body:
            # Could happen for integral with everything zero and
            # optimized away
            quadparts = []
        else:
            num_points = quadrature_rule.points.shape[0]
            iq = self.backend.symbols.quadrature_loop_index()
            quadparts = [L.ForRange(iq, 0, num_points, body=body)]

        return preparts, quadparts

    def generate_piecewise_partition(self, quadrature_rule):
        L = self.backend.language

        # Get annotated graph of factorisation
        F = self.ir.integrand[quadrature_rule]["factorization"]

        arraysymbol = L.Symbol(f"sp_{quadrature_rule.id()}")
        parts = self.generate_partition(arraysymbol, F, "piecewise", None)
        parts = L.commented_code_list(
            parts, f"Quadrature loop independent computations for quadrature rule {quadrature_rule.id()}")

        return parts

    def generate_varying_partition(self, quadrature_rule):
        L = self.backend.language

        # Get annotated graph of factorisation
        F = self.ir.integrand[quadrature_rule]["factorization"]

        arraysymbol = L.Symbol(f"sv_{quadrature_rule.id()}")
        parts = self.generate_partition(arraysymbol, F, "varying", quadrature_rule)
        parts = L.commented_code_list(
            parts, f"Varying computations for quadrature rule {quadrature_rule.id()}")
        return parts

    def generate_partition(self, symbol, F, mode, quadrature_rule):
        L = self.backend.language

        definitions = []
        intermediates = []

        use_symbol_array = True

        for i, attr in F.nodes.items():
            if attr['status'] != mode:
                continue
            v = attr['expression']
            mt = attr.get('mt')

            # Generate code only if the expression is not already in
            # cache
            if not self.get_var(quadrature_rule, v):
                if v._ufl_is_literal_:
                    vaccess = self.backend.ufl_to_language.get(v)
                elif mt is not None:
                    # All finite element based terminals have table
                    # data, as well as some, but not all, of the
                    # symbolic geometric terminals
                    tabledata = attr.get('tr')

                    # Backend specific modified terminal translation
                    vaccess = self.backend.access.get(mt.terminal, mt, tabledata, quadrature_rule)
                    vdef = self.backend.definitions.get(mt.terminal, mt, tabledata, quadrature_rule, vaccess)

                    # Store definitions of terminals in list
                    assert isinstance(vdef, list)
                    definitions.extend(vdef)
                else:
                    # Get previously visited operands
                    vops = [self.get_var(quadrature_rule, op) for op in v.ufl_operands]

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
                        # This removes the need to handle boolean
                        # intermediate variables. With tensor-valued
                        # conditionals it may not be optimal but we let
                        # the compiler take responsibility for
                        # optimizing those cases.
                        vaccess = vexpr
                    elif any(op._ufl_is_literal_ for op in v.ufl_operands):
                        # Skip intermediates for e.g. -2.0*x,
                        # resulting in lines like z = y + -2.0*x
                        vaccess = vexpr
                    else:
                        # Record assignment of vexpr to intermediate variable
                        j = len(intermediates)
                        if use_symbol_array:
                            vaccess = symbol[j]
                            intermediates.append(L.Assign(vaccess, vexpr))
                        else:
                            vaccess = L.Symbol("%s_%d" % (symbol.name, j))
                            intermediates.append(L.VariableDecl("const ufc_scalar_t", vaccess, vexpr))

                # Store access node for future reference
                self.set_var(quadrature_rule, v, vaccess)

        # Join terminal computation, array of intermediate expressions,
        # and intermediate computations
        parts = []
        if definitions:
            parts += definitions
        if intermediates:
            if use_symbol_array:
                padlen = self.ir.params["padlen"]
                parts += [L.ArrayDecl("ufc_scalar_t", symbol, len(intermediates), padlen=padlen)]
            parts += intermediates
        return parts

    def generate_dofblock_partition(self, quadrature_rule):
        block_contributions = self.ir.integrand[quadrature_rule]["block_contributions"]
        preparts = []
        quadparts = []
        blocks = [(blockmap, blockdata)
                  for blockmap, contributions in sorted(block_contributions.items())
                  for blockdata in contributions]

        for blockmap, blockdata in blocks:

            # Define code for block depending on mode
            block_preparts, block_quadparts = \
                self.generate_block_parts(quadrature_rule, blockmap, blockdata)

            # Add definitions
            preparts.extend(block_preparts)

            # Add computations
            quadparts.extend(block_quadparts)

        return preparts, quadparts

    def get_entities(self, blockdata):
        L = self.backend.language

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
            # Get the current cell or facet entity
            if blockdata.is_uniform:
                # uniform, i.e. constant across facets
                entity = L.LiteralInt(0)
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

    def get_arg_factors(self, blockdata, block_rank, quadrature_rule, iq, indices):
        arg_factors = []
        for i in range(block_rank):
            mad = blockdata.ma_data[i]
            td = mad.tabledata
            scope = self.ir.integrand[quadrature_rule]["modified_arguments"]
            mt = scope[mad.ma_index]

            # Translate modified terminal to code
            # TODO: Move element table access out of backend?
            #       Not using self.backend.access.argument() here
            #       now because it assumes too much about indices.

            table = self.backend.symbols.element_table(td, self.ir.entitytype, mt.restriction)

            assert td.ttype != "zeros"

            if td.ttype == "ones":
                arg_factor = 1
            else:
                # Assuming B sparsity follows element table sparsity
                arg_factor = table[indices[i]]
            arg_factors.append(arg_factor)
        return arg_factors

    def generate_block_parts(self, quadrature_rule, blockmap, blockdata):
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

        iq = self.backend.symbols.quadrature_loop_index()

        # Override dof index with quadrature loop index for arguments
        # with quadrature element, to index B like B[iq*num_dofs + iq]
        arg_indices = tuple(self.backend.symbols.argument_loop_index(i) for i in range(block_rank))
        B_indices = []
        for i in range(block_rank):
            B_indices.append(arg_indices[i])
        B_indices = list(B_indices)

        # Get factor expression
        F = self.ir.integrand[quadrature_rule]["factorization"]

        if len(blockdata.factor_indices_comp_indices) > 1:
            raise RuntimeError("Code generation for non-scalar integrals unsupported")

        # We have scalar integrand here, take just the factor index
        factor_index = blockdata.factor_indices_comp_indices[0][0]

        v = F.nodes[factor_index]['expression']
        f = self.get_var(quadrature_rule, v)

        # Quadrature weight was removed in representation, add it back now
        if self.ir.integral_type in ufl.custom_integral_types:
            weights = self.backend.symbols.custom_weights_table()
            weight = weights[iq]
        else:
            weights = self.backend.symbols.weights_table(quadrature_rule)
            weight = weights[iq]

        # Define fw = f * weight
        assert not blockdata.transposed, "Not handled yet"

        fw_rhs = L.float_product([f, weight])
        if not isinstance(fw_rhs, L.Product):
            fw = fw_rhs
        else:
            # Define and cache scalar temp variable
            key = (quadrature_rule, factor_index, blockdata.all_factors_piecewise)
            fw, defined = self.get_temp_symbol("fw", key)
            if not defined:
                quadparts.append(L.VariableDecl("const ufc_scalar_t", fw, fw_rhs))

        # Naively accumulate integrand for this block in the innermost
        # loop
        assert not blockdata.transposed
        A_shape = self.ir.tensor_shape

        Asym = self.backend.symbols.element_tensor()
        A = L.FlattenedArray(Asym, dims=A_shape)

        # Check if DOFs in dofrange are equally spaced
        expand_loop = False
        for i, bm in enumerate(blockmap):
            for a, b in zip(bm[1:-1], bm[2:]):
                if b - a != bm[1] - bm[0]:
                    expand_loop = True
                    break
            else:
                continue
            break

        if expand_loop:
            # If DOFs in dofrange are not equally spaced, then expand
            # out the for loop
            for A_indices, B_indices in zip(itertools.product(*blockmap),
                                            itertools.product(*[range(len(b)) for b in blockmap])):
                quadparts += [
                    L.AssignAdd(
                        A[A_indices],
                        L.float_product([fw] + self.get_arg_factors(
                            blockdata, block_rank,
                            quadrature_rule, iq, B_indices)
                        )
                    )
                ]
        else:
            # Fetch code to access modified arguments
            arg_factors = self.get_arg_factors(blockdata, block_rank, quadrature_rule, iq, B_indices)

            B_rhs = L.float_product([fw] + arg_factors)
            A_indices = []

            for bm, index in zip(blockmap, arg_indices):
                # TODO: switch order here? (optionally)
                offset = bm[0]
                if len(bm) == 1:
                    A_indices.append(index + offset)
                else:
                    block_size = bm[1] - bm[0]
                    A_indices.append(block_size * index + offset)

            body = L.AssignAdd(A[A_indices], B_rhs)

            for i in reversed(range(block_rank)):
                body = L.ForRange(B_indices[i], 0, blockdims[i], body=body)
            quadparts += [body]

        return preparts, quadparts
