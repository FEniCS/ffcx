# Copyright (C) 2026 FEniCS Project
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Code generation for the interpolations evaluated in a kernel.

The expression an interpolation interpolates is compiled into an expression
kernel of its own, evaluated at the interpolation points of the target element,
see `ffcx.analysis.analyze_ufl_objects`. Integral and expression kernels call it
before anything that reads the interpolation:

- A proxy coefficient takes the place of an interpolation of an expression. Its
  degrees of freedom are the interpolation matrix of the target element applied
  to the values the kernel writes.
- An interpolation of an expression that is only linear in an argument has its
  element table built for each cell, by contracting the values the kernel writes
  for each degree of freedom of the argument with the target element table.
"""

from __future__ import annotations

import typing

import numpy as np

import ffcx.codegeneration.lnodes as L

if typing.TYPE_CHECKING:
    from ffcx.analysis import ProxyCoefficient
    from ffcx.codegeneration.backend import FFCXBackend
    from ffcx.ir.elementtables import InterpolatedTableT
    from ffcx.ir.representation import InterpolationIR, InterpolationKernelIR


def generate_interpolations(
    ir: InterpolationIR,
    backend: FFCXBackend,
    interpolated_tables: dict[str, InterpolatedTableT],
) -> list[L.LNode]:
    """Generate the evaluation of the interpolations in a kernel.

    Args:
        ir: The interpolations evaluated in the kernel.
        backend: The backend of the kernel.
        interpolated_tables: The element tables built for each cell, by name.

    Returns:
        Code to place before anything that reads a proxy coefficient or an
        element table built for each cell.
    """
    return _proxy_coefficients(ir.proxies, backend) + _interpolated_tables(
        ir.tables, backend, interpolated_tables
    )


def _call_kernel(
    kernel: InterpolationKernelIR, label: str, backend: FFCXBackend
) -> tuple[L.Symbol, list[L.Declaration], list[L.LNode]]:
    """Evaluate an interpolated expression by calling its expression kernel.

    The coefficients and constants that the kernel reads are packed into arrays
    of their own, in the order that the kernel reads them.

    Args:
        kernel: The kernel to call.
        label: Prefix of the names of the arrays declared.
        backend: The backend of the calling kernel.

    Returns:
        The array the kernel writes, ``[point][component][dof]...``, and the
        declarations and statements that fill it.
    """
    symbols = backend.symbols
    index = L.Symbol("ai", dtype=L.DataType.INT)
    declarations: list[L.Declaration] = []
    statements: list[L.LNode] = []

    def pack(name: str, source: L.Symbol, blocks: list[tuple[int, int]]) -> L.Symbol:
        """Copy blocks ``(offset, size)`` of ``source`` one after the other."""
        array = L.Symbol(name, dtype=L.DataType.SCALAR)
        size = sum(size for _, size in blocks)
        # A zero length array is not valid C, and is never read
        declarations.append(L.ArrayDecl(array, sizes=max(size, 1), values=None if size else [0.0]))
        position = 0
        for offset, size in blocks:
            copy = L.Assign(array[position + index], source[offset + index])
            statements.append(L.ForRange(index, 0, size, [copy]))
            position += size
        return array

    coefficients = pack(
        f"{label}_w",
        symbols.coefficients,
        [(symbols.coefficient_offsets[c], c.ufl_element().dim) for c in kernel.coefficients],
    )
    constants = pack(
        f"{label}_c",
        symbols.constants,
        [
            (symbols.original_constant_offsets[c], int(np.prod(c.ufl_shape, dtype=int)))
            for c in kernel.constants
        ],
    )

    values = L.Symbol(f"{label}_values", dtype=L.DataType.SCALAR)
    size = int(np.prod(kernel.shape))
    declarations.append(L.ArrayDecl(values, sizes=size))
    # An expression kernel adds to the values it writes
    statements.append(L.ForRange(index, 0, size, [L.Assign(values[index], 0.0)]))
    statements.append(
        L.Statement(
            L.CallOp(
                kernel.name,
                (
                    values,
                    coefficients,
                    constants,
                    symbols.coordinate_dofs,
                    symbols.entity_local_index,
                    symbols.quadrature_permutation,
                    L.Symbol("custom_data", dtype=L.DataType.SCALAR),
                ),
            )
        )
    )
    return values, declarations, statements


def _proxy_coefficients(
    kernels: list[InterpolationKernelIR], backend: FFCXBackend
) -> list[L.LNode]:
    """Generate the degrees of freedom of the proxy coefficients.

    They are written where they are read from, see
    `ffcx.codegeneration.symbols.FFCXBackendSymbols.proxy_coefficient_dof_access`.
    """
    if not kernels:
        return []
    pw = backend.symbols.proxy_coefficients
    offsets = backend.symbols.proxy_coefficient_offsets
    size = max(offset + proxy.ufl_element().dim for proxy, offset in offsets.items())

    parts: list[L.LNode] = [L.ArrayDecl(pw, sizes=size)]
    for i, kernel in enumerate(kernels):
        label = f"proxy{i}"
        values, declarations, statements = _call_kernel(kernel, label, backend)
        proxy = typing.cast("ProxyCoefficient", kernel.interpolation)
        statements += _interpolate(proxy, values, pw, offsets[proxy], label, declarations)
        parts.append(
            L.Section(f"Proxy coefficient {i}", statements, declarations, input=[], output=[])
        )
    return parts


def _interpolate(
    proxy: ProxyCoefficient,
    values: L.Symbol,
    dofs: L.Symbol,
    offset: int,
    label: str,
    declarations: list[L.Declaration],
) -> list[L.LNode]:
    """Apply the interpolation matrix of the target element.

    Args:
        proxy: The proxy coefficient.
        values: The values of the interpolated expression, ``[point][component]``.
        dofs: The array to write the degrees of freedom of the proxy to.
        offset: Where the degrees of freedom of the proxy start in ``dofs``.
        label: Prefix of the names of the arrays declared.
        declarations: Declarations to add to.
    """
    element = proxy.ufl_element()
    # For a blocked element, this is the scalar sub-element that each block
    # component is interpolated with.
    basix_element = element.basix_element
    block_size = element.block_size
    j = L.Symbol("aj", dtype=L.DataType.INT)
    if basix_element.interpolation_is_identity:
        return [L.ForRange(j, 0, element.dim, [L.Assign(dofs[offset + j], values[j])])]

    # The matrix is ordered [dof][component][point], while the values are
    # [point][component][block component].
    matrix = basix_element.interpolation_matrix
    num_dofs = matrix.shape[0]
    value_size = int(np.prod(basix_element.value_shape, dtype=int))
    matrix = matrix.reshape(num_dofs, value_size, -1).transpose(0, 2, 1).reshape(num_dofs, -1)
    M = L.Symbol(f"{label}_M", dtype=L.DataType.REAL)
    declarations.append(L.ArrayDecl(M, values=matrix, const=True))

    k = L.Symbol("ak", dtype=L.DataType.INT)
    if block_size == 1:
        dof, value = dofs[offset + j], values[k]
    else:
        b = L.Symbol("ab", dtype=L.DataType.INT)
        dof, value = dofs[offset + j * block_size + b], values[k * block_size + b]
    body: list[L.LNode] = [
        L.Assign(dof, 0.0),
        L.ForRange(k, 0, matrix.shape[1], [L.AssignAdd(dof, M[j][k] * value)]),
    ]
    if block_size > 1:
        body = [L.ForRange(b, 0, block_size, body)]
    return [L.ForRange(j, 0, num_dofs, body)]


def _interpolated_tables(
    kernels: list[InterpolationKernelIR],
    backend: FFCXBackend,
    interpolated_tables: dict[str, InterpolatedTableT],
) -> list[L.LNode]:
    """Generate the element tables that are built for each cell.

    Interpolating an expression that is only linear in an argument brings
    coefficients into the interpolation, so the argument's table depends on the
    cell. The expression kernel evaluates the interpolated expression at the
    target element's interpolation points for each dof of the argument, and
    that is contracted with the target element table times the interpolation
    matrix.
    """
    needed = {table.proxy for table in interpolated_tables.values()}
    parts: list[L.LNode] = []
    evaluated = {}
    for i, kernel in enumerate(kernels):
        if kernel.interpolation not in needed:
            continue
        values, declarations, statements = _call_kernel(kernel, f"interpolation{i}", backend)
        evaluated[kernel.interpolation] = (values, kernel.shape)
        parts.append(
            L.Section(
                f"Evaluate interpolated expression {i}",
                statements,
                declarations,
                input=[],
                output=[],
            )
        )

    for name, data in sorted(interpolated_tables.items()):
        values, (_, value_size, *dof_dims) = evaluated[data.proxy]
        num_quadrature_points, num_points = data.contraction.shape

        table = L.Symbol(name, dtype=L.DataType.SCALAR)
        backend.symbols.element_tables[name] = table
        contraction = L.Symbol(f"{name}_M", dtype=L.DataType.REAL)

        aq = L.Symbol("aq", dtype=L.DataType.INT)
        ap = L.Symbol("ap", dtype=L.DataType.INT)
        dof_indices = [L.Symbol(f"ad{axis}", dtype=L.DataType.INT) for axis in range(len(dof_dims))]
        table_sizes: tuple[int, ...]
        # The kernel writes the expression as [point][component][dof]... .
        block = int(np.prod(dof_dims))
        flat_index = ap * (value_size * block) + data.flat_component * block
        stride = block
        for index, dim in zip(dof_indices, dof_dims):
            stride //= dim
            flat_index = flat_index + index * stride

        if len(dof_dims) == 1:
            # Read back through `table_access`, which indexes a table by
            # permutation and entity first. Both are singletons here: the
            # table is the same for every entity and is not permuted.
            entry = table[0][0][aq]
            table_sizes = (1, 1, num_quadrature_points, *dof_dims)
        else:
            # A dense element tensor block, read directly in `get_arg_factors`
            # by quadrature point and one index per argument.
            entry = table[aq]
            table_sizes = (num_quadrature_points, *dof_dims)
        for index in dof_indices:
            entry = entry[index]

        body: list[L.LNode] = [
            L.Assign(entry, 0.0),
            L.ForRange(
                ap, 0, num_points, [L.AssignAdd(entry, contraction[aq][ap] * values[flat_index])]
            ),
        ]
        for index, dim in zip(reversed(dof_indices), reversed(dof_dims)):
            body = [L.ForRange(index, 0, dim, body)]
        parts.append(
            L.Section(
                name=f"Build interpolated table {name}",
                statements=[L.ForRange(aq, 0, num_quadrature_points, body)],
                declarations=[
                    L.ArrayDecl(contraction, values=data.contraction, const=True),
                    L.ArrayDecl(table, sizes=table_sizes),
                ],
                input=[],
                output=[],
            )
        )

    return parts
