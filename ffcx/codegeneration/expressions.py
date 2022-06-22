# Copyright (C) 2019 Michal Habera
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import collections
import logging
from typing import Dict, Set, Any, DefaultDict
from itertools import product

import ufl

from ffcx.codegeneration import geometry
from ffcx.codegeneration import expressions_template
from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.C.format_lines import format_indented_lines
from ffcx.codegeneration.C.cnodes import CNode
from ffcx.ir.representation import ir_expression
from ffcx.naming import cdtype_to_numpy

logger = logging.getLogger("ffcx")


def generator(ir, parameters):
    """Generate UFC code for an expression."""
    logger.info("Generating code for expression:")
    logger.info(f"--- points: {ir.points}")
    logger.info(f"--- name: {ir.name}")

    factory_name = ir.name

    # Format declaration
    declaration = expressions_template.declaration.format(
        factory_name=factory_name, name_from_uflfile=ir.name_from_uflfile)

    backend = FFCXBackend(ir, parameters)
    L = backend.language
    eg = ExpressionGenerator(ir, backend)

    d = {}
    d["name_from_uflfile"] = ir.name_from_uflfile
    d["factory_name"] = ir.name

    parts = eg.generate()

    body = format_indented_lines(parts.cs_format(), 1)
    d["tabulate_expression"] = body

    if len(ir.original_coefficient_positions) > 0:
        d["original_coefficient_positions"] = f"original_coefficient_positions_{ir.name}"
        d["original_coefficient_positions_init"] = L.ArrayDecl(
            "static int", f"original_coefficient_positions_{ir.name}",
            values=ir.original_coefficient_positions, sizes=len(ir.original_coefficient_positions))
    else:
        d["original_coefficient_positions"] = L.Null()
        d["original_coefficient_positions_init"] = ""

    d["points_init"] = L.ArrayDecl(
        "static double", f"points_{ir.name}", values=ir.points.flatten(), sizes=ir.points.size)
    d["points"] = L.Symbol(f"points_{ir.name}")

    if len(ir.expression_shape) > 0:
        d["value_shape_init"] = L.ArrayDecl(
            "static int", f"value_shape_{ir.name}", values=ir.expression_shape, sizes=len(ir.expression_shape))
        d["value_shape"] = f"value_shape_{ir.name}"
    else:
        d["value_shape_init"] = ""
        d["value_shape"] = L.Null()

    d["num_components"] = len(ir.expression_shape)
    d["num_coefficients"] = len(ir.coefficient_numbering)
    d["num_constants"] = len(ir.constant_names)
    d["num_points"] = ir.points.shape[0]
    d["topological_dimension"] = ir.points.shape[1]
    d["scalar_type"] = parameters["scalar_type"]
    d["np_scalar_type"] = cdtype_to_numpy(parameters["scalar_type"])
    d["rank"] = len(ir.tensor_shape)

    if len(ir.coefficient_names) > 0:
        d["coefficient_names_init"] = L.ArrayDecl(
            "static const char*", f"coefficient_names_{ir.name}", values=ir.coefficient_names,
            sizes=len(ir.coefficient_names))
        d["coefficient_names"] = f"coefficient_names_{ir.name}"
    else:
        d["coefficient_names_init"] = ""
        d["coefficient_names"] = L.Null()

    if len(ir.constant_names) > 0:
        d["constant_names_init"] = L.ArrayDecl(
            "static const char*", f"constant_names_{ir.name}", values=ir.constant_names,
            sizes=len(ir.constant_names))
        d["constant_names"] = f"constant_names_{ir.name}"
    else:
        d["constant_names_init"] = ""
        d["constant_names"] = L.Null()

    code = []

    # FIXME: Should be handled differently, revise how
    # ufcx_function_space is generated (also for ufcx_form)
    for (name, (element, dofmap, cmap_family, cmap_degree)) in ir.function_spaces.items():
        code += [f"static ufcx_function_space function_space_{name}_{ir.name_from_uflfile} ="]
        code += ["{"]
        code += [f".finite_element = &{element},"]
        code += [f".dofmap = &{dofmap},"]
        code += [f".geometry_family = \"{cmap_family}\","]
        code += [f".geometry_degree = {cmap_degree}"]
        code += ["};"]

    d["function_spaces_alloc"] = L.StatementList(code)
    d["function_spaces"] = ""

    if len(ir.function_spaces) > 0:
        d["function_spaces"] = f"function_spaces_{ir.name}"
        d["function_spaces_init"] = L.ArrayDecl("ufcx_function_space*", f"function_spaces_{ir.name}", values=[
                                                L.AddressOf(L.Symbol(f"function_space_{name}_{ir.name_from_uflfile}"))
                                                for (name, _) in ir.function_spaces.items()],
                                                sizes=len(ir.function_spaces))
    else:
        d["function_spaces"] = L.Null()
        d["function_spaces_init"] = ""

    # Check that no keys are redundant or have been missed
    from string import Formatter
    fields = [fname for _, fname, _, _ in Formatter().parse(expressions_template.factory) if fname]

    assert set(fields) == set(d.keys()), "Mismatch between keys in template and in formattting dict"

    # Format implementation code
    implementation = expressions_template.factory.format_map(d)

    return declaration, implementation


class ExpressionGenerator:
    def __init__(self, ir: ir_expression, backend: FFCXBackend):

        if len(list(ir.integrand.keys())) != 1:
            raise RuntimeError("Only one set of points allowed for expression evaluation")

        self.ir = ir
        self.backend = backend
        self.scope: Dict[Any, CNode] = {}
        self._ufl_names: Set[Any] = set()
        self.symbol_counters: DefaultDict[Any, int] = collections.defaultdict(int)
        self.shared_symbols: Dict[Any, Any] = {}
        self.quadrature_rule = list(self.ir.integrand.keys())[0]

    def generate(self):
        L = self.backend.language

        parts = []

        parts += self.generate_element_tables()
        # Generate the tables of geometry data that are needed
        parts += self.generate_geometry_tables()
        parts += self.generate_piecewise_partition()

        all_preparts = []
        all_quadparts = []

        preparts, quadparts = self.generate_quadrature_loop()
        all_preparts += preparts
        all_quadparts += quadparts

        # Collect parts before, during, and after quadrature loops
        parts += all_preparts
        parts += all_quadparts

        return L.StatementList(parts)

    def generate_geometry_tables(self):
        """Generate static tables of geometry data."""
        L = self.backend.language

        # Currently we only support circumradius
        ufl_geometry = {
            ufl.geometry.ReferenceCellVolume: "reference_cell_volume"
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
        """Generate tables of FE basis evaluated at specified points."""
        L = self.backend.language
        parts = []

        tables = self.ir.unique_tables

        padlen = self.ir.params["padlen"]
        table_names = sorted(tables)

        scalar_type = self.backend.access.parameters["scalar_type"]

        for name in table_names:
            table = tables[name]
            decl = L.ArrayDecl(
                f"static const {scalar_type}", name, table.shape, table, padlen=padlen)
            parts += [decl]

        # Add leading comment if there are any tables
        parts = L.commented_code_list(parts, [
            "Precomputed values of basis functions",
            "FE* dimensions: [entities][points][dofs]",
        ])
        return parts

    def generate_quadrature_loop(self):
        """Generate quadrature loop for this quadrature rule.

        In the context of expressions quadrature loop is not accumulated.

        """
        L = self.backend.language

        # Generate varying partition
        body = self.generate_varying_partition()
        body = L.commented_code_list(
            body, f"Points loop body setup quadrature loop {self.quadrature_rule.id()}")

        # Generate dofblock parts, some of this
        # will be placed before or after quadloop
        preparts, quadparts = \
            self.generate_dofblock_partition()
        body += quadparts

        # Wrap body in loop or scope
        if not body:
            # Could happen for integral with everything zero and optimized away
            quadparts = []
        else:
            iq = self.backend.symbols.quadrature_loop_index()
            num_points = self.quadrature_rule.points.shape[0]
            quadparts = [L.ForRange(iq, 0, num_points, body=body)]

        return preparts, quadparts

    def generate_varying_partition(self):
        """Generate factors of blocks which are not cellwise constant."""
        L = self.backend.language

        # Get annotated graph of factorisation
        F = self.ir.integrand[self.quadrature_rule]["factorization"]

        arraysymbol = L.Symbol(f"sv_{self.quadrature_rule.id()}")
        parts = self.generate_partition(arraysymbol, F, "varying")
        parts = L.commented_code_list(
            parts, f"Unstructured varying computations for quadrature rule {self.quadrature_rule.id()}")
        return parts

    def generate_piecewise_partition(self):
        """Generate factors of blocks which are constant (i.e. do not depent on quadrature points)."""
        L = self.backend.language

        # Get annotated graph of factorisation
        F = self.ir.integrand[self.quadrature_rule]["factorization"]

        arraysymbol = L.Symbol("sp")
        parts = self.generate_partition(arraysymbol, F, "piecewise")
        parts = L.commented_code_list(parts, "Unstructured piecewise computations")
        return parts

    def generate_dofblock_partition(self):
        """Generate assignments of blocks multiplied with their factors into final tensor A."""
        block_contributions = self.ir.integrand[self.quadrature_rule]["block_contributions"]

        preparts = []
        quadparts = []

        blocks = [(blockmap, blockdata)
                  for blockmap, contributions in sorted(block_contributions.items())
                  for blockdata in contributions]

        for blockmap, blockdata in blocks:

            # Define code for block depending on mode
            block_preparts, block_quadparts = \
                self.generate_block_parts(blockmap, blockdata)

            # Add definitions
            preparts.extend(block_preparts)

            # Add computations
            quadparts.extend(block_quadparts)

        return preparts, quadparts

    def generate_block_parts(self, blockmap, blockdata):
        """Generate and return code parts for a given block."""
        L = self.backend.language

        # The parts to return
        preparts = []
        quadparts = []

        block_rank = len(blockmap)
        blockdims = tuple(len(dofmap) for dofmap in blockmap)

        ttypes = blockdata.ttypes
        if "zeros" in ttypes:
            raise RuntimeError("Not expecting zero arguments to be left in dofblock generation.")

        arg_indices = tuple(self.backend.symbols.argument_loop_index(i) for i in range(block_rank))

        F = self.ir.integrand[self.quadrature_rule]["factorization"]

        assert not blockdata.transposed, "Not handled yet"
        components = ufl.product(self.ir.expression_shape)

        num_points = self.quadrature_rule.points.shape[0]
        A_shape = self.ir.tensor_shape
        Asym = self.backend.symbols.element_tensor()
        A = L.FlattenedArray(Asym, dims=[num_points, components] + A_shape)

        iq = self.backend.symbols.quadrature_loop_index()

        # Check if DOFs in dofrange are equally spaced.
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
            # If DOFs in dofrange are not equally spaced, then expand out the for loop
            for A_indices, B_indices in zip(product(*blockmap),
                                            product(*[range(len(b)) for b in blockmap])):
                B_indices = tuple([iq] + list(B_indices))
                A_indices = tuple([iq] + A_indices)
                for fi_ci in blockdata.factor_indices_comp_indices:
                    f = self.get_var(F.nodes[fi_ci[0]]["expression"])
                    arg_factors = self.get_arg_factors(blockdata, block_rank, B_indices)
                    Brhs = L.float_product([f] + arg_factors)
                    quadparts.append(L.AssignAdd(A[(A_indices[0], fi_ci[1]) + A_indices[1:]], Brhs))
        else:

            # Prepend dimensions of dofmap block with free index
            # for quadrature points and expression components
            B_indices = tuple([iq] + list(arg_indices))

            # Fetch code to access modified arguments
            # An access to FE table data
            arg_factors = self.get_arg_factors(blockdata, block_rank, B_indices)

            # TODO: handle non-contiguous dof ranges

            A_indices = []
            for bm, index in zip(blockmap, arg_indices):
                # TODO: switch order here? (optionally)
                offset = bm[0]
                if len(bm) == 1:
                    A_indices.append(index + offset)
                else:
                    block_size = bm[1] - bm[0]
                    A_indices.append(block_size * index + offset)
            A_indices = tuple([iq] + A_indices)

            # Multiply collected factors
            # For each component of the factor expression
            # add result inside quadloop
            body = []

            for fi_ci in blockdata.factor_indices_comp_indices:
                f = self.get_var(F.nodes[fi_ci[0]]["expression"])
                Brhs = L.float_product([f] + arg_factors)
                body.append(L.AssignAdd(A[(A_indices[0], fi_ci[1]) + A_indices[1:]], Brhs))

            for i in reversed(range(block_rank)):
                body = L.ForRange(
                    B_indices[i + 1], 0, blockdims[i], body=body)
            quadparts += [body]

        return preparts, quadparts

    def get_arg_factors(self, blockdata, block_rank, indices):
        """Get argument factors (i.e. blocks).

        Parameters
        ----------
        blockdata
        block_rank
        indices
            Indices used to index element tables

        """
        L = self.backend.language

        arg_factors = []
        for i in range(block_rank):
            mad = blockdata.ma_data[i]
            td = mad.tabledata
            mt = self.ir.integrand[self.quadrature_rule]["modified_arguments"][mad.ma_index]

            table = self.backend.symbols.element_table(td, self.ir.entitytype, mt.restriction)

            assert td.ttype != "zeros"

            if td.ttype == "ones":
                arg_factor = L.LiteralFloat(1.0)
            else:
                arg_factor = table[indices[i + 1]]
            arg_factors.append(arg_factor)
        return arg_factors

    def new_temp_symbol(self, basename):
        """Create a new code symbol named basename + running counter."""
        L = self.backend.language
        name = "%s%d" % (basename, self.symbol_counters[basename])
        self.symbol_counters[basename] += 1
        return L.Symbol(name)

    def get_var(self, v):
        if v._ufl_is_literal_:
            return self.backend.ufl_to_language.get(v)
        f = self.scope.get(v)
        return f

    def generate_partition(self, symbol, F, mode):
        """Generate computations of factors of blocks."""
        L = self.backend.language

        definitions = []
        pre_definitions = dict()
        intermediates = []

        use_symbol_array = True

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
                vaccess = self.backend.access.get(mt.terminal, mt, tabledata, 0)

                predef, vdef = self.backend.definitions.get(mt.terminal, mt, tabledata, 0, vaccess)
                if predef:
                    pre_definitions[str(predef[0].symbol.name)] = predef

                # Store definitions of terminals in list
                assert isinstance(vdef, list)
                definitions.extend(vdef)
            else:
                # Get previously visited operands
                vops = [self.get_var(op) for op in v.ufl_operands]

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
                    if use_symbol_array:
                        vaccess = symbol[j]
                        intermediates.append(L.Assign(vaccess, vexpr))
                    else:
                        scalar_type = self.backend.access.parameters["scalar_type"]
                        vaccess = L.Symbol("%s_%d" % (symbol.name, j))
                        intermediates.append(L.VariableDecl(f"const {scalar_type}", vaccess, vexpr))

            # Store access node for future reference
            self.scope[v] = vaccess

        # Join terminal computation, array of intermediate expressions,
        # and intermediate computations
        parts = []

        for _, definition in pre_definitions.items():
            parts += definition

        if definitions:
            parts += definitions

        if intermediates:
            if use_symbol_array:
                scalar_type = self.backend.access.parameters["scalar_type"]
                parts += [L.ArrayDecl(scalar_type, symbol, len(intermediates))]
            parts += intermediates
        return parts
