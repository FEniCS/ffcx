# Copyright (C) 2026 Jørgen S. Dokken and Paul T. Kühner
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Code generation for the interpolations evaluated in a kernel.

The expression an interpolation interpolates is compiled into an expression
kernel of its own, evaluated at the interpolation points of the target element,
see `ffcx.analysis.analyze_ufl_objects`. Integral and expression kernels call it
before anything that reads the interpolation.

A proxy coefficient takes the place of an interpolation of an expression. Its
degrees of freedom are the interpolation matrix of the target element applied
to the values the kernel writes.
"""

from __future__ import annotations

import typing

import basix.ufl
import numpy as np
import ufl

import ffcx.codegeneration.lnodes as L

if typing.TYPE_CHECKING:
    from ffcx.analysis import ProxyCoefficient
    from ffcx.codegeneration.backend import FFCXBackend
    from ffcx.ir.representation import InterpolationIR, InterpolationKernelIR


def generate_interpolations(ir: InterpolationIR, backend: FFCXBackend) -> list[L.LNode]:
    """Generate the evaluation of the interpolations in a kernel.

    Args:
        ir: The interpolations evaluated in the kernel.
        backend: The backend of the kernel.

    Returns:
        Code to place before anything that reads a proxy coefficient.
    """
    return _proxy_coefficients(ir.proxies, backend)


def _call_kernel(
    kernel: InterpolationKernelIR, label: str, backend: FFCXBackend, side: int = 0
) -> tuple[L.Symbol, list[L.Declaration], list[L.LNode]]:
    """Evaluate an interpolated expression by calling its expression kernel.

    The coefficients and constants that the kernel reads are packed into arrays
    of their own, in the order that the kernel reads them.

    Args:
        kernel: The kernel to call.
        label: Prefix of the names of the arrays declared.
        backend: The backend of the calling kernel.
        side: The cell to evaluate on, 0 or, on an interior facet, 1 for the
            cell restricted to '-'.

    Returns:
        The array the kernel writes, ``[point][component]``, and the
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
        [
            (symbols.coefficient_offsets[c] + side * c.ufl_element().dim, c.ufl_element().dim)
            for c in kernel.coefficients
        ],
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
                    _shift(symbols.coordinate_dofs, side * _num_coordinate_dofs(kernel)),
                    _shift(symbols.entity_local_index, side),
                    _shift(symbols.quadrature_permutation, side),
                    L.Symbol("custom_data", dtype=L.DataType.SCALAR),
                ),
            )
        )
    )
    return values, declarations, statements


def _shift(pointer: L.Symbol, offset: int) -> L.LExpr:
    """The pointer advanced by ``offset`` entries, to the data of the second cell."""
    return pointer if offset == 0 else L.Add(pointer, offset)


def _num_coordinate_dofs(kernel: InterpolationKernelIR) -> int:
    """The number of entries of `coordinate_dofs` for one cell, three per node."""
    domain = ufl.domain.extract_unique_domain(kernel.interpolation)
    assert isinstance(domain, ufl.Mesh)
    element = domain.ufl_coordinate_element()
    assert isinstance(element, basix.ufl._ElementBase)
    return 3 * int(element.basix_element.dim)


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
    # On an interior facet the dofs are built for both cells, the second
    # after the first, as the coefficients are passed.
    sides = 2 if backend.access.integral_type.startswith("interior_facet") else 1
    size = max(offset + sides * proxy.ufl_element().dim for proxy, offset in offsets.items())

    parts: list[L.LNode] = [L.ArrayDecl(pw, sizes=size)]
    for i, kernel in enumerate(kernels):
        proxy = typing.cast("ProxyCoefficient", kernel.interpolation)
        for side in range(sides):
            label = f"proxy{i}" if side == 0 else f"proxy{i}_{side}"
            values, declarations, statements = _call_kernel(kernel, label, backend, side)
            offset = offsets[proxy] + side * proxy.ufl_element().dim
            statements += _interpolate(proxy, values, pw, offset, label, declarations)
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
