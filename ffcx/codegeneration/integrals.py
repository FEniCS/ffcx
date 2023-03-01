# Copyright (C) 2015-2021 Martin Sandve Alnæs, Michal Habera, Igor Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import collections
import logging
from typing import Any, Dict, List, Set, Tuple
import numpy
import copy

import ufl
from ffcx.codegeneration import geometry
from ffcx.codegeneration import integrals_template as ufcx_integrals
from ffcx.codegeneration.indices import create_dof_index, create_quadrature_index
from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.optimise import fuse_loops, sum_factorise
from ffcx.codegeneration.C.cnodes import BinOp, CNode
from ffcx.codegeneration.C.format_lines import format_indented_lines
from ffcx.ir.elementtables import piecewise_ttypes
from ffcx.ir.integral import BlockData
from ffcx.ir.representationutils import QuadratureRule
from ffcx.naming import cdtype_to_numpy, scalar_to_value_type

logger = logging.getLogger("ffcx")


def generator(ir, options):
    logger.info("Generating code for integral:")
    logger.info(f"--- type: {ir.integral_type}")
    logger.info(f"--- name: {ir.name}")

    """Generate code for an integral."""
    factory_name = ir.name

    # Format declaration
    declaration = ufcx_integrals.declaration.format(factory_name=factory_name)

    # Create FFCx C backend
    backend = FFCXBackend(ir, options)

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

    lang = backend.language
    if len(ir.enabled_coefficients) > 0:
        code["enabled_coefficients_init"] = lang.ArrayDecl(
            "bool", f"enabled_coefficients_{ir.name}", values=ir.enabled_coefficients,
            sizes=len(ir.enabled_coefficients))
        code["enabled_coefficients"] = f"enabled_coefficients_{ir.name}"
    else:
        code["enabled_coefficients_init"] = ""
        code["enabled_coefficients"] = lang.Null()

    code["additional_includes_set"] = set()  # FIXME: Get this out of code[]
    code["tabulate_tensor"] = body

    if options["tabulate_tensor_void"]:
        code["tabulate_tensor"] = ""

    code["result_needs_permuting"] = 0
    code["result_permutations"] = "NULL"

    implementation = ufcx_integrals.factory.format(
        factory_name=factory_name,
        enabled_coefficients=code["enabled_coefficients"],
        enabled_coefficients_init=code["enabled_coefficients_init"],
        tabulate_tensor=code["tabulate_tensor"],
        needs_facet_permutations="true" if ir.needs_facet_permutations else "false",
        scalar_type=options["scalar_type"],
        geom_type=scalar_to_value_type(options["scalar_type"]),
        np_scalar_type=cdtype_to_numpy(options["scalar_type"]),
        coordinate_element=lang.AddressOf(lang.Symbol(ir.coordinate_element)))

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

        # Set of operator names code has been generated for, used in the
        # end for selecting necessary includes
        self._ufl_names = set()

        # Initialize lookup tables for variable scopes
        self.init_scopes()

        # Cache
        self.shared_symbols = {}

        # Literals
        self.literals = {}

        # Set of counters used for assigning names to intermediate
        # variables
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
        lang = self.backend.language
        batch_size = self.backend.access.options["batch_size"]
        if batch_size > 1:
            if v._ufl_is_literal_:
                if v not in self.literals:
                    self.literals[v] = lang.Symbol("literal" + str(len(self.literals)))
                return self.literals[v]
        else:
            if v._ufl_is_literal_:
                return self.backend.ufl_to_language.get(v)

        f = self.scopes[quadrature_rule].get(v)
        if f is None:
            f = self.scopes[None].get(v)
        return f

    def new_temp_symbol(self, basename):
        """Create a new code symbol named basename + running counter."""
        lang = self.backend.language
        name = "%s%d" % (basename, self.symbol_counters[basename])
        self.symbol_counters[basename] += 1
        return lang.Symbol(name)

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

        Assumes that the code returned from here will be wrapped in a
        context that matches a suitable version of the UFC
        tabulate_tensor signatures.
        """
        lang = self.backend.language

        # Assert that scopes are empty: expecting this to be called only
        # once
        assert not any(d for d in self.scopes.values())

        parts = []
        scalar_type = self.backend.access.options["scalar_type"]
        value_type = scalar_to_value_type(scalar_type)
        alignment = self.ir.options['assume_aligned']
        if alignment != -1:
            scalar_type = self.backend.access.options["scalar_type"]
            parts += [lang.VerbatimStatement(f"A = ({scalar_type}*)__builtin_assume_aligned(A, {alignment});"),
                      lang.VerbatimStatement(f"w = (const {scalar_type}*)__builtin_assume_aligned(w, {alignment});"),
                      lang.VerbatimStatement(f"c = (const {scalar_type}*)__builtin_assume_aligned(c, {alignment});"),
                      lang.VerbatimStatement(f"coordinate_dofs = (const {value_type}*)__builtin_assume_aligned(coordinate_dofs, {alignment});")]  # noqa

        # Generate the tables of quadrature points and weights
        parts += self.generate_quadrature_tables(value_type)

        # Generate the tables of basis function values and
        # pre-integrated blocks
        parts += self.generate_element_tables(value_type)

        # Generate the tables of geometry data that are needed
        parts += self.generate_geometry_tables(value_type)

        # Loop generation code will produce parts to go before
        # quadloops, to define the quadloops, and to go after the
        # quadloops
        all_preparts = []
        all_quadparts = []

        # Pre-definitions are collected across all quadrature loops to
        # improve re-use and avoid name clashes
        all_predefinitions = dict()
        for rule in self.ir.integrand.keys():
            # Generate code to compute piecewise constant scalar factors
            all_preparts += self.generate_piecewise_partition(rule)

            # Generate code to integrate reusable blocks of final
            # element tensor
            pre_definitions, preparts, quadparts = self.generate_quadrature_loop(rule)

            all_preparts += preparts
            all_quadparts += quadparts
            all_predefinitions.update(pre_definitions)

        pre_definitions = fuse_loops(lang, all_predefinitions)
        parts += lang.commented_code_list(pre_definitions,
                                          "Pre-definitions of modified terminals to enable unit-stride access")

        for literal in self.literals.keys():
            scalar_type = self.backend.access.options["scalar_type"]
            batch_size = self.backend.access.options["batch_size"]
            if batch_size > 1:
                scalar_type += str(batch_size)
            values = self.backend.ufl_to_language.get(literal)
            init_list = lang.as_symbol(lang.build_1d_initializer_list(numpy.array([values] * batch_size), str))
            all_preparts.insert(0, lang.VariableDecl(f"const {scalar_type}", self.literals[literal], init_list))

        # Collect parts before, during, and after quadrature loops
        parts += all_preparts
        parts += all_quadparts

        return lang.StatementList(parts)

    def generate_quadrature_tables(self, value_type: str) -> List[str]:
        """Generate static tables of quadrature points and weights."""
        lang = self.backend.language

        parts: List[str] = []

        # No quadrature tables for custom (given argument) or point
        # (evaluation in single vertex)
        skip = ufl.custom_integral_types + ufl.measure.point_integral_types
        if self.ir.integral_type in skip:
            return parts

        padlen = self.ir.options["padlen"]

        # Loop over quadrature rules
        for quadrature_rule, integrand in self.ir.integrand.items():
            num_points = quadrature_rule.weights.shape[0]
            weights = quadrature_rule.weights
            if quadrature_rule.has_tensor_factors:
                weights = quadrature_rule.tensor_factors[0][1]
                num_points = weights.shape[0]

            # Generate quadrature weights array
            wsym = self.backend.symbols.weights_table(quadrature_rule)
            parts += [
                lang.ArrayDecl(
                    f"static const {value_type}", wsym, num_points, weights, padlen=padlen)
            ]

        # Add leading comment if there are any tables
        parts = lang.commented_code_list(parts, "Quadrature rules")
        return parts

    def generate_geometry_tables(self, float_type: str):
        """Generate static tables of geometry data."""
        lang = self.backend.language

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
        cells: Dict[Any, Set[Any]] = {t: set() for t in ufl_geometry.keys()}

        for integrand in self.ir.integrand.values():
            for attr in integrand["factorization"].nodes.values():
                mt = attr.get("mt")
                if mt is not None:
                    t = type(mt.terminal)
                    if t in ufl_geometry:
                        cells[t].add(ufl.domain.extract_unique_domain(mt.terminal).ufl_cell().cellname())

        parts = []
        for i, cell_list in cells.items():
            for c in cell_list:
                parts.append(geometry.write_table(lang, ufl_geometry[i], c))

        return parts

    def generate_element_tables(self, float_type):
        """Generate static tables with precomputed element basis function values in quadrature points."""
        lang = self.backend.language
        parts = []
        tables = self.ir.unique_tables
        table_types = self.ir.unique_table_types
        padlen = self.ir.options["padlen"]
        if self.ir.integral_type in ufl.custom_integral_types:
            # Define only piecewise tables
            table_names = [name for name in sorted(tables) if table_types[name] in piecewise_ttypes]
        else:
            # Define all tables
            table_names = sorted(tables)

        for name in table_names:
            table = tables[name]
            parts += self.declare_table(name, table, padlen, float_type)

        # Add leading comment if there are any tables
        parts = lang.commented_code_list(parts, [
            "Precomputed values of basis functions and precomputations",
            "FE* dimensions: [permutation][entities][points][dofs]"])
        return parts

    def declare_table(self, name, table, padlen, value_type: str):
        """Declare a table.

        If the dof dimensions of the table have dof rotations, apply
        these rotations.

        """
        lang = self.backend.language

        return [lang.ArrayDecl(
            f"static const {value_type}", name, table.shape, table, padlen=padlen)]

    def generate_quadrature_loop(self, quadrature_rule: QuadratureRule):
        """Generate quadrature loop with for this quadrature_rule."""
        lang = self.backend.language
        iq = create_quadrature_index(lang, quadrature_rule)

        quadrature_values = []

        # Generate varying partition
        pre_definitions, body, intermediates = self.generate_varying_partition(quadrature_rule)
        quadrature_values.append(intermediates)

        body = lang.commented_code_list(
            body, f"Quadrature loop body setup for quadrature rule {quadrature_rule.id()}")

        # quadrature_values = self.generate_quadrature_values(quadrature_rule)

        # Generate dofblock parts, some of this will be placed before or
        # after quadloop
        preparts, quadparts, intermediates = \
            self.generate_dofblock_partition(quadrature_rule)
        quadrature_values.append(intermediates)
        if iq.dim > 1:
            quadrature_values = [lang.NestedForRange([iq], quadrature_values)]
        body += quadrature_values
        body += quadparts

        # Wrap body in loop or scope
        if not body:
            # Could happen for integral with everything zero and
            # optimized away
            quadparts = []
        else:
            quadparts = body
            if iq.dim == 1:
                quadparts = [lang.NestedForRange([iq], body)]

        return pre_definitions, preparts, quadparts

    def generate_piecewise_partition(self, quadrature_rule):
        lang = self.backend.language

        # Get annotated graph of factorisation
        F = self.ir.integrand[quadrature_rule]["factorization"]

        arraysymbol = lang.Symbol(f"sp_{quadrature_rule.id()}")
        pre_definitions, parts, intermediates = self.generate_partition(arraysymbol, F, "piecewise", None)

        assert len(pre_definitions) == 0, "Quadrature independent code should have not pre-definitions"

        code = [parts, intermediates]
        code = lang.commented_code_list(
            code, f"Quadrature loop independent computations for quadrature rule {quadrature_rule.id()}")

        return code

    def generate_varying_partition(self, quadrature_rule):
        lang = self.backend.language

        # Get annotated graph of factorisation
        F = self.ir.integrand[quadrature_rule]["factorization"]

        arraysymbol = lang.Symbol(f"sv_{quadrature_rule.id()}")
        pre_definitions, parts, intermediates = self.generate_partition(arraysymbol, F, "varying", quadrature_rule)
        parts = lang.commented_code_list(
            parts, f"Varying computations for quadrature rule {quadrature_rule.id()}")

        return pre_definitions, parts, intermediates

    def generate_partition(self, symbol, F, mode, quadrature_rule):
        lang = self.backend.language

        definitions = dict()
        pre_definitions = dict()
        intermediates = []
        quadrature_values = []

        batch_size = self.backend.access.options["batch_size"]
        scalar_type = self.backend.access.options["scalar_type"]
        if batch_size > 1:
            scalar_type += str(batch_size)

        use_symbol_array = True

        iq = create_quadrature_index(lang, quadrature_rule)

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
                    quadrature_values.append(vaccess)

                    predef, vdef = self.backend.definitions.get(mt.terminal, mt, tabledata, quadrature_rule, vaccess)
                    if predef:
                        access = predef[0].symbol.name
                        pre_definitions[str(access)] = predef

                    # Store definitions of terminals in list
                    assert isinstance(vdef, list)
                    definitions[str(vaccess)] = vdef
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
                        if isinstance(vexpr, lang.Call) and batch_size > 1:
                            j = len(intermediates)
                            vaccess = symbol[j]
                            for b in range(batch_size):
                                argument = copy.deepcopy(vexpr.arguments[0])
                                new_vexpr = copy.deepcopy(vexpr)
                                new_vexpr.arguments[0] = argument[b]
                                intermediates.append(lang.Assign(vaccess[b], new_vexpr))
                        else:
                            if use_symbol_array:
                                j = len(intermediates)
                                vaccess = symbol[j]
                                intermediates.append(lang.Assign(vaccess, vexpr))
                            else:
                                vaccess = lang.Symbol("%s_%d" % (symbol.name, j))
                                intermediates.append(lang.VariableDecl(f"const {scalar_type}", vaccess, vexpr))

                # Store access node for future reference
                self.set_var(quadrature_rule, v, vaccess)

        # Join terminal computation, array of intermediate expressions,
        # and intermediate computations
        parts = []
        parts += fuse_loops(lang, definitions)

        if intermediates:
            if use_symbol_array:
                padlen = self.ir.options["padlen"]
                scalar_type = self.backend.access.options["scalar_type"]
                batch_size = self.backend.access.options["batch_size"]
                if batch_size > 1:
                    scalar_type += str(batch_size)
                declaration = [lang.ArrayDecl(scalar_type, symbol, len(intermediates))]
                intermediates.insert(0, declaration)

        def substitute_acess(intermediates, quadrature_values):
            if (isinstance(intermediates, list)):
                new_intermediates = []
                for intermediate in intermediates:
                    expression = substitute_acess(intermediate, quadrature_values)
                    new_intermediates.append(expression)
                return new_intermediates
            elif isinstance(intermediates, lang.Call):
                return intermediates
            elif hasattr(intermediates, 'rhs'):
                intermediates.rhs = substitute_acess(intermediates.rhs, quadrature_values)
                intermediates.lhs = substitute_acess(intermediates.lhs, quadrature_values)
                return intermediates
            else:
                if intermediates in quadrature_values:
                    if hasattr(intermediates, 'global_idx'):
                        return intermediates[intermediates.global_idx]
                    else:
                        return intermediates[iq.global_idx()]
                else:
                    return intermediates
        if iq.dim > 1:
            intermediates = substitute_acess(intermediates, quadrature_values)

        return pre_definitions, parts, intermediates

    def generate_dofblock_partition(self, quadrature_rule: QuadratureRule):
        block_contributions = self.ir.integrand[quadrature_rule]["block_contributions"]
        preparts = []
        quadparts = []
        blocks = [(blockmap, blockdata)
                  for blockmap, contributions in sorted(block_contributions.items())
                  for blockdata in contributions]

        block_groups = collections.defaultdict(list)

        # Group loops by blockmap, in Vector elements each component has
        # a different blockmap
        for blockmap, blockdata in blocks:
            scalar_blockmap = []
            assert len(blockdata.ma_data) == len(blockmap)
            for i, b in enumerate(blockmap):
                bs = blockdata.ma_data[i].tabledata.block_size
                offset = blockdata.ma_data[i].tabledata.offset
                b = tuple([(idx - offset) // bs for idx in b])
                scalar_blockmap.append(b)
            block_groups[tuple(scalar_blockmap)].append(blockdata)

        for blockmap in block_groups:
            block_preparts, block_quadparts, intermediates = \
                self.generate_block_parts(quadrature_rule, blockmap, block_groups[blockmap])

            # Add definitions
            preparts.extend(block_preparts)

            # Add computations
            quadparts.extend(block_quadparts)

        return preparts, quadparts, intermediates

    def get_arg_factors(self, blockdata, block_rank, quadrature_rule, iq, indices):

        scope = self.ir.integrand[quadrature_rule]["modified_arguments"]
        entitytype = self.ir.entitytype

        arg_factors = []
        for i in range(block_rank):
            mad = blockdata.ma_data[i]
            td = mad.tabledata
            mt = scope[mad.ma_index]

            assert td.ttype != "zeros"
            if td.ttype == "ones":
                arg_factor = 1
            else:
                arg_factor = self.backend.symbols.table_access(
                    td, entitytype, mt.restriction, iq, indices[i])
            arg_factors.append(arg_factor)

        return arg_factors

    def generate_block_parts(self, quadrature_rule: QuadratureRule, blockmap: Tuple, blocklist: List[BlockData]):
        """Generate and return code parts for a given block.

        Returns parts occurring before, inside, and after the quadrature
        loop identified by the quadrature rule.

        Should be called with quadrature_rule=None for
        quadloop-independent blocks.
        """
        lang = self.backend.language

        # The parts to return
        preparts: List[CNode] = []
        quadparts: List[CNode] = []
        intermediates: List[CNode] = []

        # RHS expressions grouped by LHS "dofmap"
        rhs_expressions = collections.defaultdict(list)

        block_rank = len(blockmap)
        blockdims = tuple(len(dofmap) for dofmap in blockmap)

        scalar_type = self.backend.access.options["scalar_type"]
        iq = create_quadrature_index(lang, quadrature_rule)

        for blockdata in blocklist:
            # Override dof index with quadrature loop index for arguments
            # with quadrature element, to index B like B[iq*num_dofs + iq]
            B_indices = []
            args = ["i", "j", "k", "l"]
            for i in range(block_rank):
                table_ref = blockdata.ma_data[i].tabledata
                B_indices.append(create_dof_index(lang, table_ref, args[i]))

            ttypes = blockdata.ttypes
            if "zeros" in ttypes:
                raise RuntimeError("Not expecting zero arguments to be left in dofblock generation.")

            if len(blockdata.factor_indices_comp_indices) > 1:
                raise RuntimeError("Code generation for non-scalar integrals unsupported")

            # We have scalar integrand here, take just the factor index
            factor_index = blockdata.factor_indices_comp_indices[0][0]

            # Get factor expression
            F = self.ir.integrand[quadrature_rule]["factorization"]
            v = F.nodes[factor_index]['expression']
            f = self.get_var(quadrature_rule, v)

            # Quadrature weight was removed in representation, add it back now
            if self.ir.integral_type in ufl.custom_integral_types:
                weights = self.backend.symbols.custom_weights_table()
                weight = weights[iq.global_idx()]
            else:
                weights = self.backend.symbols.weights_table(quadrature_rule)
                weight = [weights[iq.local_idx(i)] for i in range(iq.dim)]
                weight = lang.Product(weight)

            # Define fw = f * weight
            fw_rhs = lang.float_product([f, weight])
            if not isinstance(fw_rhs, lang.Product):
                fw = fw_rhs
            else:
                # Define and cache scalar temp variable
                key = (quadrature_rule, factor_index, blockdata.all_factors_piecewise)
                fw, defined = self.get_temp_symbol("fw", key)
                if not defined:
                    scalar_type = self.backend.access.options["scalar_type"]
                    batch_size = self.backend.access.options["batch_size"]
                    if batch_size > 1:
                        scalar_type += str(batch_size)
                    if iq.dim > 1:
                        preparts += [lang.ArrayDecl(scalar_type, fw, iq.ranges, values=[0.0])]
                        intermediates += [lang.Assign(fw[iq], fw_rhs)]
                    else:
                        preparts += [lang.VariableDecl(scalar_type, fw, value=0.0)]
                        intermediates += [lang.Assign(fw, fw_rhs)]

            assert not blockdata.transposed, "Not handled yet"
            A_shape = self.ir.tensor_shape

            Asym = self.backend.symbols.element_tensor()
            A = lang.FlattenedArray(Asym, dims=A_shape)

            # Fetch code to access modified arguments
            arg_factors = self.get_arg_factors(blockdata, block_rank, quadrature_rule, iq, B_indices)

            if iq.dim > 1:
                B_rhs = lang.float_product([fw[iq]] + arg_factors)
            else:
                B_rhs = lang.float_product([fw] + arg_factors)

            A_indices = []
            for i in range(block_rank):
                index = B_indices[i]
                tabledata = blockdata.ma_data[i].tabledata
                offset = tabledata.offset
                if len(blockmap[i]) == 1:
                    A_indices.append(index.global_idx() + offset)
                else:
                    block_size = tabledata.block_size
                    A_indices.append(block_size * index.global_idx() + offset)
            rhs_expressions[tuple(A_indices)].append(B_rhs)

        # List of statements to keep in the inner loop
        keep = collections.defaultdict(list)
        # List of temporary array declarations
        pre_loop: List[CNode] = []
        # List of loop invariant expressions to hoist
        hoist: List[BinOp] = []

        for indices in rhs_expressions:
            keep[indices] = rhs_expressions[indices]

        hoist_code: List[CNode] = [lang.ForRange(B_indices[0], 0, blockdims[0], body=hoist)] if hoist else []

        body: List[CNode] = []

        for indices in keep:
            body.append(lang.AssignAdd(A[indices], lang.Sum(keep[indices])))

        if iq.dim > 1:
            B_indices.insert(0, iq)

        body = sum_factorise(lang, lang.NestedForRange(B_indices, body), scalar_type)

        quadparts += pre_loop
        quadparts += hoist_code
        quadparts += [body]

        return preparts, quadparts, intermediates
