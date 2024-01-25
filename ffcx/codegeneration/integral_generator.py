# Copyright (C) 2015-2023 Martin Sandve Alnæs, Michal Habera, Igor Baratta, Chris Richardson
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import collections
import logging
from typing import List, Tuple

import ffcx.codegeneration.lnodes as L
import ufl
from ffcx.ir.integral import BlockDataT
from ffcx.ir.representationutils import QuadratureRule
from ffcx.codegeneration.optimizer import optimize
# from numbers import Integral
import ffcx.codegeneration.generator as gen

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

        # Cache
        self.temp_symbols = {}

    def get_temp_symbol(self, tempname, key):
        key = (tempname,) + key
        s = self.temp_symbols.get(key)
        defined = s is not None
        if not defined:
            name = f"{tempname}{len(self.temp_symbols)}"
            s = L.Symbol(name, dtype=L.DataType.SCALAR)
            self.temp_symbols[key] = s
        return s, defined

    def generate(self):
        """Generate entire tabulate_tensor body.

        Assumes that the code returned from here will be wrapped in a
        context that matches a suitable version of the UFC
        tabulate_tensor signatures.
        """
        parts = []

        # Generate the tables of quadrature points and weights
        parts += gen.generate_quadrature_tables(self.ir, self.backend)

        # Generate the tables of basis function values and
        # pre-integrated blocks
        parts += gen.generate_element_tables(self.ir, self.backend)

        # Generate the tables of geometry data that are needed
        parts += gen.generate_geometry_tables(self.ir, self.backend)

        # Loop generation code will produce parts to go before
        # quadloops, to define the quadloops, and to go after the
        # quadloops
        all_preparts = []
        all_quadparts = []

        # Pre-definitions are collected across all quadrature loops to
        # improve re-use and avoid name clashes
        for rule in self.ir.integrand.keys():
            # Generate code to compute piecewise constant scalar factors
            all_preparts += gen.generate_piecewise_partition(self.ir, self.backend, rule)

            # Generate code to integrate reusable blocks of final
            # element tensor
            all_quadparts += self.generate_quadrature_loop(rule)

        # Collect parts before, during, and after quadrature loops
        parts += all_preparts
        parts += all_quadparts

        return L.StatementList(parts)

    def generate_quadrature_loop(self, quadrature_rule: QuadratureRule):
        """Generate quadrature loop with for this quadrature_rule."""
        # Generate varying partition
        definitions, intermediates_0 = gen.generate_varying_partition(self.ir, self.backend, quadrature_rule)

        # Generate dofblock parts, some of this will be placed before or
        # after quadloop
        tensor_comp, intermediates_fw = self.generate_dofblock_partition(quadrature_rule)
        assert all([isinstance(tc, L.Section) for tc in tensor_comp])

        # Check if we only have Section objects
        inputs = []
        for definition in definitions:
            assert isinstance(definition, L.Section)
            inputs += definition.output

        # Create intermediates section
        output = []
        declarations = []
        for fw in intermediates_fw:
            assert isinstance(fw, L.VariableDecl)
            output += [fw.symbol]
            declarations += [L.VariableDecl(fw.symbol, 0)]
            intermediates_0 += [L.Assign(fw.symbol, fw.value)]
        intermediates = [L.Section("Intermediates", intermediates_0, declarations, inputs, output)]

        iq = self.backend.quadrature_indices[quadrature_rule]
        code = definitions + intermediates + tensor_comp
        code = optimize(code, quadrature_rule)

        return [L.create_nested_for_loops([iq], code)]

    def generate_dofblock_partition(self, quadrature_rule: QuadratureRule):
        block_contributions = self.ir.integrand[quadrature_rule]["block_contributions"]
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

        intermediates = []
        for blockmap in block_groups:
            for block in block_groups[blockmap]:
                block_quadparts, intermediate = self.generate_block_parts(
                    quadrature_rule, blockmap, block)
                intermediates += intermediate
                quadparts.extend(block_quadparts)

        return quadparts, intermediates

    def generate_block_parts(self,
                             quadrature_rule: QuadratureRule,
                             blockmap: Tuple,
                             blockdata: BlockDataT):
        """Generate and return code parts for a given block.

        Returns parts occurring before, inside, and after the quadrature
        loop identified by the quadrature rule.

        Should be called with quadrature_rule=None for
        quadloop-independent blocks.
        """
        # The parts to return
        quadparts: List[L.LNode] = []
        intermediates: List[L.LNode] = []
        tables = []
        vars = []

        # RHS expressions grouped by LHS "dofmap"
        rhs_expressions = collections.defaultdict(list)

        block_rank = len(blockmap)
        iq = self.backend.quadrature_indices[quadrature_rule]

        B_indices = []
        for i in range(block_rank):
            table_ref = blockdata.ma_data[i].tabledata
            symbol = self.backend.symbols.argument_loop_index(i)
            index = self.backend.create_dof_index(table_ref, symbol)
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

        v: ufl.core.Expr = F.nodes[factor_index]['expression']
        f: L.Symbol = self.backend.get_var(quadrature_rule, v)

        # Quadrature weight was removed in representation, add it back now
        weights = self.backend.symbols.weights_table(quadrature_rule)
        weight = weights[iq.global_index]

        # Define fw = f * weight
        fw_rhs = L.float_product([f, weight])

        # Define and cache scalar temp variable
        key = (quadrature_rule, factor_index, blockdata.all_factors_piecewise)
        fw, defined = self.get_temp_symbol("fw", key)
        if not defined:
            input = [f, weight]
            input = [i for i in input if isinstance(i, L.Symbol)]
            output = [fw]

            # assert input and output are Symbol objects
            assert all(isinstance(i, L.Symbol) for i in input)
            assert all(isinstance(o, L.Symbol) for o in output)

            intermediates += [L.VariableDecl(fw, fw_rhs)]

        var = fw if isinstance(fw, L.Symbol) else fw.array
        vars += [var]
        assert not blockdata.transposed, "Not handled yet"

        # Fetch code to access modified arguments
        arg_factors, table = gen.get_arg_factors(
            self.backend, blockdata, block_rank, quadrature_rule, iq, B_indices)
        tables += table

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

        body: List[L.LNode] = []

        A = self.backend.symbols.element_tensor
        A_shape = self.ir.tensor_shape
        for indices in keep:
            multi_index = L.MultiIndex(list(indices), A_shape)
            for expression in keep[indices]:
                body.append(L.AssignAdd(A[multi_index], expression))

        # reverse B_indices
        B_indices = B_indices[::-1]
        body = [L.create_nested_for_loops(B_indices, body)]
        input = [*vars, *tables]
        output = [A]

        # Make sure we don't have repeated symbols in input
        input = list(set(input))

        # assert input and output are Symbol objects
        assert all(isinstance(i, L.Symbol) for i in input)
        assert all(isinstance(o, L.Symbol) for o in output)

        annotations = []
        if len(B_indices) > 1:
            annotations.append(L.Annotation.licm)

        quadparts += [L.Section("Tensor Computation", body, [], input, output, annotations)]

        return quadparts, intermediates
