# Copyright (C) 2015-2023 Martin Sandve Alnæs, Michal Habera, Igor Baratta, Chris Richardson
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import collections
import logging
from typing import Any, Dict, List, Set, Tuple

import ufl
from ffcx.codegeneration import geometry
from ffcx.ir.elementtables import piecewise_ttypes
from ffcx.ir.integral import BlockDataT
import ffcx.codegeneration.lnodes as L
from ffcx.codegeneration.lnodes import LNode
from ffcx.ir.representationutils import QuadratureRule
from ffcx.codegeneration.definitions import create_quadrature_index, create_dof_index

logger = logging.getLogger("ffcx")


class IntegralGenerator(object):
    def __init__(self, ir, backend):
        # Store ir
        self.ir = ir

        # Backend specific plugin with attributes
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
        self.temp_symbols = {}

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

        v is the ufl expression and vaccess is the LNodes
        expression to access the value in the code.

        """
        self.scopes[quadrature_rule][v] = vaccess

    def get_var(self, quadrature_rule, v):
        """Lookup ufl expression v in variable scope dicts.

        Scope is determined by quadrature rule which identifies the
        quadrature loop scope or None if outside quadrature loops.

        If v is not found in quadrature loop scope, the piecewise
        scope (None) is checked.

        Returns the LNodes expression to access the value in the code.
        """
        if v._ufl_is_literal_:
            return L.ufl_to_lnodes(v)
        f = self.scopes[quadrature_rule].get(v)
        if f is None:
            f = self.scopes[None].get(v)
        return f

    def new_temp_symbol(self, basename):
        """Create a new code symbol named basename + running counter."""
        name = "%s%d" % (basename, self.symbol_counters[basename])
        self.symbol_counters[basename] += 1
        return L.Symbol(name, dtype=L.DataType.SCALAR)

    def get_temp_symbol(self, tempname, key):
        key = (tempname,) + key
        s = self.temp_symbols.get(key)
        defined = s is not None
        if not defined:
            s = self.new_temp_symbol(tempname)
            self.temp_symbols[key] = s
        return s, defined

    def generate(self):
        """Generate entire tabulate_tensor body.

        Assumes that the code returned from here will be wrapped in a
        context that matches a suitable version of the UFC
        tabulate_tensor signatures.
        """
        # Assert that scopes are empty: expecting this to be called only
        # once
        assert not any(d for d in self.scopes.values())

        parts = []

        # Generate the tables of quadrature points and weights
        parts += self.generate_quadrature_tables()

        # Generate the tables of basis function values and
        # pre-integrated blocks
        parts += self.generate_element_tables()

        # Generate the tables of geometry data that are needed
        parts += self.generate_geometry_tables()

        # Loop generation code will produce parts to go before
        # quadloops, to define the quadloops, and to go after the
        # quadloops
        all_preparts = []
        all_quadparts = []

        # Pre-definitions are collected across all quadrature loops to
        # improve re-use and avoid name clashes
        all_predefinitions = []
        for rule in self.ir.integrand.keys():
            # Generate code to compute piecewise constant scalar factors
            all_preparts += self.generate_piecewise_partition(rule)

            # Generate code to integrate reusable blocks of final
            # element tensor
            pre_definitions, preparts, quadparts = self.generate_quadrature_loop(rule)
            all_preparts += preparts
            all_quadparts += quadparts
            all_predefinitions += pre_definitions

        all_predefinitions = list(set(all_predefinitions))

        # Collect pre-definitions of type ArrayDecl
        declarations = [p for p in all_predefinitions if isinstance(p, L.ArrayDecl)]
        for_loops = [p for p in all_predefinitions if isinstance(p, L.ForRange)]

        # remove repeated for loops (should be done with section)
        # FIXME: no two sections with same inputs and outputs should be generated
        for_loops_dict = {}
        for for_loop in for_loops:
            for_loop_key = (for_loop.body.statements[0].expr.lhs.array)
            if for_loop_key not in for_loops_dict:
                for_loops_dict[for_loop_key] = for_loop
        for_loops = list(for_loops_dict.values())

        parts += L.commented_code_list(declarations + for_loops,
                                       "Pre-definitions of modified terminals to enable unit-stride access")

        # Collect parts before, during, and after quadrature loops
        parts += all_preparts
        parts += all_quadparts

        return L.StatementList(parts)

    def generate_quadrature_tables(self):
        """Generate static tables of quadrature points and weights."""
        parts = []

        # No quadrature tables for custom (given argument) or point
        # (evaluation in single vertex)
        skip = ufl.custom_integral_types + ufl.measure.point_integral_types
        if self.ir.integral_type in skip:
            return parts

        # Loop over quadrature rules
        for quadrature_rule, integrand in self.ir.integrand.items():
            # Generate quadrature weights array
            wsym = self.backend.symbols.weights_table(quadrature_rule)
            parts += [L.ArrayDecl(wsym, values=quadrature_rule.weights, const=True)]

        # Add leading comment if there are any tables
        parts = L.commented_code_list(parts, "Quadrature rules")
        return parts

    def generate_geometry_tables(self):
        """Generate static tables of geometry data."""
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
                parts.append(geometry.write_table(ufl_geometry[i], c))

        return parts

    def generate_element_tables(self):
        """Generate static tables with precomputed element basisfunction values in quadrature points."""
        parts = []
        tables = self.ir.unique_tables
        table_types = self.ir.unique_table_types
        if self.ir.integral_type in ufl.custom_integral_types:
            # Define only piecewise tables
            table_names = [name for name in sorted(tables) if table_types[name] in piecewise_ttypes]
        else:
            # Define all tables
            table_names = sorted(tables)

        for name in table_names:
            table = tables[name]
            parts += self.declare_table(name, table)

        # Add leading comment if there are any tables
        parts = L.commented_code_list(parts, [
            "Precomputed values of basis functions and precomputations",
            "FE* dimensions: [permutation][entities][points][dofs]"])
        return parts

    def declare_table(self, name, table):
        """Declare a table.

        If the dof dimensions of the table have dof rotations, apply
        these rotations.

        """
        table_symbol = L.Symbol(name, dtype=L.DataType.REAL)
        self.backend.symbols.element_tables[name] = table_symbol
        return [L.ArrayDecl(table_symbol, values=table, const=True)]

    def generate_quadrature_loop(self, quadrature_rule: QuadratureRule):
        """Generate quadrature loop with for this quadrature_rule."""
        # Generate varying partition
        _, body = self.generate_varying_partition(quadrature_rule)

        body = L.commented_code_list(body, f"Quadrature loop body setup for quadrature rule {quadrature_rule.id()}")

        # Generate dofblock parts, some of this will be placed before or
        # after quadloop
        preparts, quadparts = self.generate_dofblock_partition(quadrature_rule)
        body += quadparts

        # Wrap body in loop or scope
        if not body:
            # Could happen for integral with everything zero and
            # optimized away
            quadparts = []
        else:
            iq_symbol = self.backend.symbols.quadrature_loop_index
            iq = create_quadrature_index(quadrature_rule, iq_symbol)
            quadparts = [L.create_nested_for_loops([iq], body)]

        return [], preparts, quadparts

    def generate_piecewise_partition(self, quadrature_rule):
        # Get annotated graph of factorisation
        F = self.ir.integrand[quadrature_rule]["factorization"]

        arraysymbol = L.Symbol(f"sp_{quadrature_rule.id()}", dtype=L.DataType.SCALAR)
        pre_definitions, parts = self.generate_partition(arraysymbol, F, "piecewise", None)
        assert len(pre_definitions) == 0, "Quadrature independent code should have not pre-definitions"
        parts = L.commented_code_list(
            parts, f"Quadrature loop independent computations for quadrature rule {quadrature_rule.id()}")

        return parts

    def generate_varying_partition(self, quadrature_rule):

        # Get annotated graph of factorisation
        F = self.ir.integrand[quadrature_rule]["factorization"]
        arraysymbol = L.Symbol(f"sv_{quadrature_rule.id()}", dtype=L.DataType.SCALAR)
        pre_definitions, parts = self.generate_partition(arraysymbol, F, "varying", quadrature_rule)
        parts = L.commented_code_list(parts, f"Varying computations for quadrature rule {quadrature_rule.id()}")

        return pre_definitions, parts

    def generate_partition(self, symbol, F, mode, quadrature_rule):

        definitions = []
        intermediates = []

        for i, attr in F.nodes.items():
            if attr['status'] != mode:
                continue
            v = attr['expression']
            mt = attr.get('mt')

            # Generate code only if the expression is not already in cache
            if not self.get_var(quadrature_rule, v):
                if v._ufl_is_literal_:
                    vaccess = L.ufl_to_lnodes(v)
                elif mt:
                    assert mt is not None
                    tabledata = attr.get('tr')

                    # Backend specific modified terminal translation
                    vaccess = self.backend.access.get(mt, tabledata, quadrature_rule)
                    vdef = self.backend.definitions.get(mt, tabledata, quadrature_rule, vaccess)

                    if isinstance(vdef, L.Section):
                        # Store definitions of terminals in list
                        definitions += [vdef]
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
                    vexpr = L.ufl_to_lnodes(v, *vops)

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
                        vaccess = L.Symbol("%s_%d" % (symbol.name, j), dtype=L.DataType.SCALAR)
                        intermediates.append(L.VariableDecl(vaccess, vexpr))

                # Store access node for future reference
                self.set_var(quadrature_rule, v, vaccess)

        # Join terminal computation, array of intermediate expressions,
        # and intermediate computations
        parts = []
        parts += definitions

        if intermediates:
            intermediates = [L.Section(f"Intermediate computations for {mode} computations", intermediates)]
            parts += intermediates

        return [], parts

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
            block_preparts, block_quadparts = self.generate_block_parts(
                quadrature_rule, blockmap, block_groups[blockmap])

            # Add definitions
            preparts.extend(block_preparts)

            # Add computations
            quadparts.extend(block_quadparts)

        return preparts, quadparts

    def get_arg_factors(self, blockdata, block_rank, quadrature_rule, iq, indices):
        arg_factors = []
        tables = []
        for i in range(block_rank):
            mad = blockdata.ma_data[i]
            td = mad.tabledata
            scope = self.ir.integrand[quadrature_rule]["modified_arguments"]
            mt = scope[mad.ma_index]
            arg_tables = []

            # Translate modified terminal to code
            # TODO: Move element table access out of backend?
            #       Not using self.backend.access.argument() here
            #       now because it assumes too much about indices.

            assert td.ttype != "zeros"

            if td.ttype == "ones":
                arg_factor = 1
            else:
                # Assuming B sparsity follows element table sparsity
                arg_factor, arg_tables = self.backend.access.table_access(
                    td, self.ir.entitytype, mt.restriction, iq, indices[i])

            tables += arg_tables
            arg_factors.append(arg_factor)

        return arg_factors, tables

    def generate_block_parts(self, quadrature_rule: QuadratureRule, blockmap: Tuple, blocklist: List[BlockDataT]):
        """Generate and return code parts for a given block.

        Returns parts occurring before, inside, and after the quadrature
        loop identified by the quadrature rule.

        Should be called with quadrature_rule=None for
        quadloop-independent blocks.
        """
        # The parts to return
        preparts: List[LNode] = []
        quadparts: List[LNode] = []

        # RHS expressions grouped by LHS "dofmap"
        rhs_expressions = collections.defaultdict(list)

        block_rank = len(blockmap)
        iq_symbol = self.backend.symbols.quadrature_loop_index
        iq = create_quadrature_index(quadrature_rule, iq_symbol)

        for blockdata in blocklist:

            B_indices = []
            for i in range(block_rank):
                table_ref = blockdata.ma_data[i].tabledata
                symbol = self.backend.symbols.argument_loop_index(i)
                index = create_dof_index(table_ref, symbol)
                B_indices.append(index)

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
                weights = self.backend.symbols.custom_weights_table
                weight = weights[iq.global_index]
            else:
                weights = self.backend.symbols.weights_table(quadrature_rule)
                weight = weights[iq.global_index]

            # Define fw = f * weight
            fw_rhs = L.float_product([f, weight])
            if not isinstance(fw_rhs, L.Product):
                fw = fw_rhs
            else:
                # Define and cache scalar temp variable
                key = (quadrature_rule, factor_index, blockdata.all_factors_piecewise)
                fw, defined = self.get_temp_symbol("fw", key)
                if not defined:
                    input = [f, weight]
                    # filter only L.Symbol in input
                    input = [i for i in input if isinstance(i, L.Symbol)]
                    output = [fw]

                    # assert input and output are Symbol objects
                    assert all(isinstance(i, L.Symbol) for i in input)
                    assert all(isinstance(o, L.Symbol) for o in output)

                    quad_sec = L.Section(f"Definition of {fw}", [L.VariableDecl(fw, fw_rhs)],
                                         input=input, output=output)

                    quadparts.append(quad_sec)

            assert not blockdata.transposed, "Not handled yet"

            # Fetch code to access modified arguments
            arg_factors, tables = self.get_arg_factors(blockdata, block_rank, quadrature_rule, iq, B_indices)

            # Define B_rhs = fw * arg_factors
            B_rhs = L.float_product([fw] + arg_factors)

            A_indices = []
            for i in range(block_rank):
                index = B_indices[i]
                tabledata = blockdata.ma_data[i].tabledata
                offset = tabledata.offset
                if len(blockmap[i]) == 1:
                    A_indices.append(index.global_index + offset)
                else:
                    block_size = blockdata.ma_data[i].tabledata.block_size
                    A_indices.append(block_size * index.global_index + offset)
            rhs_expressions[tuple(A_indices)].append(B_rhs)

        # List of statements to keep in the inner loop
        keep = collections.defaultdict(list)

        for indices in rhs_expressions:
            keep[indices] = rhs_expressions[indices]

        body: List[LNode] = []

        A = self.backend.symbols.element_tensor
        A_shape = self.ir.tensor_shape
        for indices in keep:
            multi_index = L.MultiIndex(list(indices), A_shape)
            body.append(L.AssignAdd(A[multi_index], L.Sum(keep[indices])))

        # reverse B_indices
        B_indices = B_indices[::-1]
        body = [L.create_nested_for_loops(B_indices, body)]
        input = [fw, *tables]
        output = [A]

        # assert input and output are Symbol objects
        assert all(isinstance(i, L.Symbol) for i in input)
        assert all(isinstance(o, L.Symbol) for o in output)

        quadparts += [L.Section("Tensor Computation", body, input=input, output=output)]

        return preparts, quadparts
