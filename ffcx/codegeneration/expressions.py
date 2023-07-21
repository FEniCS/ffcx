# Copyright (C) 2019 Michal Habera
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import collections
import logging
from itertools import product
from typing import Any, DefaultDict, Dict, Set

import ufl
from ffcx.codegeneration import expressions_template, geometry
from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.c_implementation import c_format
from ffcx.codegeneration.C.format_lines import format_indented_lines
from ffcx.ir.representation import ExpressionIR
from ffcx.naming import cdtype_to_numpy, scalar_to_value_type
import ffcx.codegeneration.lnodes as L

logger = logging.getLogger("ffcx")


def generator(ir, options):
    """Generate UFC code for an expression."""
    logger.info("Generating code for expression:")
    logger.info(f"--- points: {ir.points}")
    logger.info(f"--- name: {ir.name}")

    factory_name = ir.name

    # Format declaration
    declaration = expressions_template.declaration.format(
        factory_name=factory_name, name_from_uflfile=ir.name_from_uflfile
    )

    backend = FFCXBackend(ir, options)
    eg = ExpressionGenerator(ir, backend)

    d = {}
    d["name_from_uflfile"] = ir.name_from_uflfile
    d["factory_name"] = ir.name

    parts = eg.generate()

    body = format_indented_lines(c_format(parts), 1)
    d["tabulate_expression"] = body

    if len(ir.original_coefficient_positions) > 0:
        d[
            "original_coefficient_positions"
        ] = f"original_coefficient_positions_{ir.name}"
        n = len(ir.original_coefficient_positions)
        originals = ", ".join(str(i) for i in ir.original_coefficient_positions)
        d[
            "original_coefficient_positions_init"
        ] = f"static int original_coefficient_positions_{ir.name}[{n}] = {{{originals}}};"

    else:
        d["original_coefficient_positions"] = "NULL"
        d["original_coefficient_positions_init"] = ""

    points = ", ".join(str(p) for p in ir.points.flatten())
    n = ir.points.size
    d["points_init"] = f"static double points_{ir.name}[{n}] = {{{points}}};"
    d["points"] = f"points_{ir.name}"

    if len(ir.expression_shape) > 0:
        n = len(ir.expression_shape)
        shape = ", ".join(str(i) for i in ir.expression_shape)
        d["value_shape_init"] = f"static int value_shape_{ir.name}[{n}] = {{{shape}}};"
        d["value_shape"] = f"value_shape_{ir.name}"
    else:
        d["value_shape_init"] = ""
        d["value_shape"] = "NULL"

    d["num_components"] = len(ir.expression_shape)
    d["num_coefficients"] = len(ir.coefficient_numbering)
    d["num_constants"] = len(ir.constant_names)
    d["num_points"] = ir.points.shape[0]
    d["topological_dimension"] = ir.points.shape[1]
    d["scalar_type"] = options["scalar_type"]
    d["geom_type"] = scalar_to_value_type(options["scalar_type"])
    d["np_scalar_type"] = cdtype_to_numpy(options["scalar_type"])

    d["rank"] = len(ir.tensor_shape)

    if len(ir.coefficient_names) > 0:
        names = ", ".join(f'"{name}"' for name in ir.coefficient_names)
        n = len(ir.coefficient_names)
        d[
            "coefficient_names_init"
        ] = f"static const char* coefficient_names_{ir.name}[{n}] = {{{names}}};"

        d["coefficient_names"] = f"coefficient_names_{ir.name}"
    else:
        d["coefficient_names_init"] = ""
        d["coefficient_names"] = "NULL"

    if len(ir.constant_names) > 0:
        names = ", ".join(f'"{name}"' for name in ir.constant_names)
        n = len(ir.constant_names)
        d[
            "constant_names_init"
        ] = f"static const char* constant_names_{ir.name}[{n}] = {{{names}}};"
        d["constant_names"] = f"constant_names_{ir.name}"
    else:
        d["constant_names_init"] = ""
        d["constant_names"] = "NULL"

    code = []

    # FIXME: Should be handled differently, revise how
    # ufcx_function_space is generated (also for ufcx_form)
    for name, (element, dofmap, cmap_family, cmap_degree) in ir.function_spaces.items():
        code += [
            f"static ufcx_function_space function_space_{name}_{ir.name_from_uflfile} ="
        ]
        code += ["{"]
        code += [f".finite_element = &{element},"]
        code += [f".dofmap = &{dofmap},"]
        code += [f'.geometry_family = "{cmap_family}",']
        code += [f".geometry_degree = {cmap_degree}"]
        code += ["};"]

    d["function_spaces_alloc"] = "\n".join(code)
    d["function_spaces"] = ""

    if len(ir.function_spaces) > 0:
        d["function_spaces"] = f"function_spaces_{ir.name}"
        fs_list = ", ".join(
            f"&function_space_{name}_{ir.name_from_uflfile}"
            for (name, _) in ir.function_spaces.items()
        )
        n = len(ir.function_spaces.items())
        d[
            "function_spaces_init"
        ] = f"ufcx_function_space* function_spaces_{ir.name}[{n}] = {{{fs_list}}};"
        print(d["function_spaces_init"])
    else:
        d["function_spaces"] = "NULL"
        d["function_spaces_init"] = ""

    # Check that no keys are redundant or have been missed
    from string import Formatter

    fields = [
        fname
        for _, fname, _, _ in Formatter().parse(expressions_template.factory)
        if fname
    ]
    assert set(fields) == set(
        d.keys()
    ), "Mismatch between keys in template and in formatting dict"

    # Format implementation code
    implementation = expressions_template.factory.format_map(d)

    return declaration, implementation


class ExpressionGenerator:
    def __init__(self, ir: ExpressionIR, backend: FFCXBackend):
        if len(list(ir.integrand.keys())) != 1:
            raise RuntimeError(
                "Only one set of points allowed for expression evaluation"
            )

        self.ir = ir
        self.backend = backend
        self.scope: Dict = {}
        self._ufl_names: Set[Any] = set()
        self.symbol_counters: DefaultDict[Any, int] = collections.defaultdict(int)
        self.shared_symbols: Dict[Any, Any] = {}
        self.quadrature_rule = list(self.ir.integrand.keys())[0]

    def generate(self):
        parts = []
        scalar_type = self.backend.access.options["scalar_type"]
        value_type = scalar_to_value_type(scalar_type)

        parts += self.generate_element_tables(value_type)
        # Generate the tables of geometry data that are needed
        parts += self.generate_geometry_tables(value_type)
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

    def generate_geometry_tables(self, float_type: str):
        """Generate static tables of geometry data."""
        # Currently we only support circumradius
        ufl_geometry = {
            ufl.geometry.ReferenceCellVolume: "reference_cell_volume",
        }
        cells: Dict[Any, Set[Any]] = {t: set() for t in ufl_geometry.keys()}

        for integrand in self.ir.integrand.values():
            for attr in integrand["factorization"].nodes.values():
                mt = attr.get("mt")
                if mt is not None:
                    t = type(mt.terminal)
                    if t in ufl_geometry:
                        cells[t].add(
                            ufl.domain.extract_unique_domain(mt.terminal)
                            .ufl_cell()
                            .cellname()
                        )

        parts = []
        for i, cell_list in cells.items():
            for c in cell_list:
                parts.append(geometry.write_table(L, ufl_geometry[i], c, float_type))

        return parts

    def generate_element_tables(self, float_type: str):
        """Generate tables of FE basis evaluated at specified points."""
        parts = []

        tables = self.ir.unique_tables
        table_names = sorted(tables)

        for name in table_names:
            table = tables[name]
            decl = L.ArrayDecl(f"static const {float_type}", name, table.shape, table)
            parts += [decl]

        # Add leading comment if there are any tables
        parts = L.commented_code_list(
            parts,
            [
                "Precomputed values of basis functions",
                "FE* dimensions: [entities][points][dofs]",
            ],
        )
        return parts

    def generate_quadrature_loop(self):
        """Generate quadrature loop for this quadrature rule.

        In the context of expressions quadrature loop is not accumulated.

        """
        # Generate varying partition
        body = self.generate_varying_partition()
        body = L.commented_code_list(
            body, f"Points loop body setup quadrature loop {self.quadrature_rule.id()}"
        )

        # Generate dofblock parts, some of this
        # will be placed before or after quadloop
        preparts, quadparts = self.generate_dofblock_partition()
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
        # Get annotated graph of factorisation
        F = self.ir.integrand[self.quadrature_rule]["factorization"]

        arraysymbol = L.Symbol(f"sv_{self.quadrature_rule.id()}")
        parts = self.generate_partition(arraysymbol, F, "varying")
        parts = L.commented_code_list(
            parts,
            f"Unstructured varying computations for quadrature rule {self.quadrature_rule.id()}",
        )
        return parts

    def generate_piecewise_partition(self):
        """Generate factors of blocks which are constant (i.e. do not depend on quadrature points)."""
        # Get annotated graph of factorisation
        F = self.ir.integrand[self.quadrature_rule]["factorization"]

        arraysymbol = L.Symbol("sp")
        parts = self.generate_partition(arraysymbol, F, "piecewise")
        parts = L.commented_code_list(parts, "Unstructured piecewise computations")
        return parts

    def generate_dofblock_partition(self):
        """Generate assignments of blocks multiplied with their factors into final tensor A."""
        block_contributions = self.ir.integrand[self.quadrature_rule][
            "block_contributions"
        ]

        preparts = []
        quadparts = []

        blocks = [
            (blockmap, blockdata)
            for blockmap, contributions in sorted(block_contributions.items())
            for blockdata in contributions
        ]

        for blockmap, blockdata in blocks:
            # Define code for block depending on mode
            block_preparts, block_quadparts = self.generate_block_parts(
                blockmap, blockdata
            )

            # Add definitions
            preparts.extend(block_preparts)

            # Add computations
            quadparts.extend(block_quadparts)

        return preparts, quadparts

    def generate_block_parts(self, blockmap, blockdata):
        """Generate and return code parts for a given block."""
        # The parts to return
        preparts = []
        quadparts = []

        block_rank = len(blockmap)
        blockdims = tuple(len(dofmap) for dofmap in blockmap)

        ttypes = blockdata.ttypes
        if "zeros" in ttypes:
            raise RuntimeError(
                "Not expecting zero arguments to be left in dofblock generation."
            )

        arg_indices = tuple(
            self.backend.symbols.argument_loop_index(i) for i in range(block_rank)
        )

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
            for A_indices, B_indices in zip(
                product(*blockmap), product(*[range(len(b)) for b in blockmap])
            ):
                B_indices = tuple([iq] + list(B_indices))
                A_indices = tuple([iq] + A_indices)
                for fi_ci in blockdata.factor_indices_comp_indices:
                    f = self.get_var(F.nodes[fi_ci[0]]["expression"])
                    arg_factors = self.get_arg_factors(blockdata, block_rank, B_indices)
                    Brhs = L.float_product([f] + arg_factors)
                    quadparts.append(
                        L.AssignAdd(A[(A_indices[0], fi_ci[1]) + A_indices[1:]], Brhs)
                    )
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
                body.append(
                    L.AssignAdd(A[(A_indices[0], fi_ci[1]) + A_indices[1:]], Brhs)
                )

            for i in reversed(range(block_rank)):
                body = L.ForRange(B_indices[i + 1], 0, blockdims[i], body=body)
            quadparts += [body]

        return preparts, quadparts

    def get_arg_factors(self, blockdata, block_rank, indices):
        """Get argument factors (i.e. blocks).

        Options
        ----------
        blockdata
        block_rank
        indices
            Indices used to index element tables

        """
        arg_factors = []
        for i in range(block_rank):
            mad = blockdata.ma_data[i]
            td = mad.tabledata
            mt = self.ir.integrand[self.quadrature_rule]["modified_arguments"][
                mad.ma_index
            ]

            table = self.backend.symbols.element_table(
                td, self.ir.entitytype, mt.restriction
            )

            assert td.ttype != "zeros"

            if td.ttype == "ones":
                arg_factor = L.LiteralFloat(1.0)
            else:
                arg_factor = table[indices[i + 1]]
            arg_factors.append(arg_factor)
        return arg_factors

    def new_temp_symbol(self, basename):
        """Create a new code symbol named basename + running counter."""
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
        definitions = []
        pre_definitions = dict()
        intermediates = []

        use_symbol_array = True

        for i, attr in F.nodes.items():
            if attr["status"] != mode:
                continue
            v = attr["expression"]
            mt = attr.get("mt")

            if v._ufl_is_literal_:
                vaccess = self.backend.ufl_to_language.get(v)
            elif mt is not None:
                # All finite element based terminals have table data, as well
                # as some, but not all, of the symbolic geometric terminals
                tabledata = attr.get("tr")

                # Backend specific modified terminal translation
                vaccess = self.backend.access.get(mt.terminal, mt, tabledata, 0)

                predef, vdef = self.backend.definitions.get(
                    mt.terminal, mt, tabledata, 0, vaccess
                )
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
                    parent_exp = F.nodes.get(pid)["expression"]
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
                        scalar_type = self.backend.access.options["scalar_type"]
                        vaccess = L.Symbol("%s_%d" % (symbol.name, j))
                        intermediates.append(
                            L.VariableDecl(f"const {scalar_type}", vaccess, vexpr)
                        )

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
                scalar_type = self.backend.access.options["scalar_type"]
                parts += [L.ArrayDecl(scalar_type, symbol, len(intermediates))]
            parts += intermediates
        return parts
