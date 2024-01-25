# Copyright (C) 2019 Michal Habera
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import logging
from itertools import product

import ffcx.codegeneration.lnodes as L
import ufl
from ffcx.codegeneration.backend import ReprManager
from ffcx.ir.representation import ExpressionIR
import ffcx.codegeneration.generator as gen

logger = logging.getLogger("ffcx")


class ExpressionGenerator:
    def __init__(self, ir: ExpressionIR, backend: ReprManager):

        if len(list(ir.integrand.keys())) != 1:
            raise RuntimeError("Only one set of points allowed for expression evaluation")

        self.ir = ir
        self.backend = backend
        self.quadrature_rule = list(self.ir.integrand.keys())[0]

    def generate(self):
        parts = []
        parts += gen.generate_element_tables(self.ir, self.backend)
        parts += gen.generate_geometry_tables(self.ir, self.backend)
        parts += gen.generate_piecewise_partition(self.ir, self.backend, self.quadrature_rule)

        piecewise, varying = self.generate_quadrature_loop()
        parts += piecewise + varying

        return L.StatementList(parts)

    def generate_quadrature_loop(self):
        """Generate quadrature loop for this quadrature rule.

        In the context of expressions quadrature loop is not accumulated.

        """
        # Generate varying partition
        definitions, inter = gen.generate_varying_partition(self.ir, self.backend, self.quadrature_rule)
        body = L.commented_code_list(
            definitions + inter, f"Points loop body setup quadrature loop {self.quadrature_rule.id()}")

        # Generate dofblock parts, some of this
        # will be placed before or after quadloop
        preparts, quadparts = \
            self.generate_dofblock_partition()
        body += quadparts

        # Wrap body in loop or scopes
        if not body:
            # Could happen for integral with everything zero and optimized away
            quadparts = []
        else:
            iq = self.backend.symbols.quadrature_loop_index
            num_points = self.quadrature_rule.points.shape[0]
            quadparts = [L.ForRange(iq, 0, num_points, body=body)]
        return preparts, quadparts

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
        A_shape = [num_points, components] + self.ir.tensor_shape
        A = self.backend.symbols.element_tensor
        iq = self.backend.symbols.quadrature_loop_index

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
                    f = self.backend.get_var(self.quadrature_rule, F.nodes[fi_ci[0]]["expression"])
                    arg_factors, _  = gen.get_arg_factors(self.backend, blockdata, block_rank, B_indices)
                    Brhs = L.float_product([f] + arg_factors)
                    multi_index = L.MultiIndex([A_indices[0], fi_ci[1]] + A_indices[1:], A_shape)
                    quadparts.append(L.AssignAdd(A[multi_index], Brhs))
        else:

            # Prepend dimensions of dofmap block with free index
            # for quadrature points and expression components
            B_indices = tuple([iq] + list(arg_indices))
            iq_ = self.backend.quadrature_indices[self.quadrature_rule]
            # Fetch code to access modified arguments
            # An access to FE table data
            arg_factors, _ = gen.get_arg_factors(self.backend, blockdata, block_rank, self.quadrature_rule, iq_, B_indices)

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
                f = self.backend.get_var(self.quadrature_rule, F.nodes[fi_ci[0]]["expression"])
                f = f[iq] if (hasattr(f, "size") and f.size() > 1) else f
                Brhs = L.float_product([f] + arg_factors)
                indices = [A_indices[0], fi_ci[1]] + list(A_indices[1:])
                multi_index = L.MultiIndex(indices, A_shape)
                body.append(L.AssignAdd(A[multi_index], Brhs))

            for i in reversed(range(block_rank)):
                body = L.ForRange(
                    B_indices[i + 1], 0, blockdims[i], body=body)
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
            mt = self.ir.integrand[self.quadrature_rule]["modified_arguments"][mad.ma_index]

            table = self.backend.symbols.element_table(td, self.ir.entitytype, mt.restriction)

            assert td.ttype != "zeros"

            if td.ttype == "ones":
                arg_factor = L.LiteralFloat(1.0)
            else:
                arg_factor = table[indices[i + 1]]
            arg_factors.append(arg_factor)
        return arg_factors
