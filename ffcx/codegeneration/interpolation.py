# Copyright (C) 2026 Jørgen S. Dokken and Paul T. Kühner
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Code generation for the sub-kernels a kernel calls.

A sub-kernel is generated ahead of the kernels that call it, in the same file,
see `ffcx.analysis.analyze_ufl_objects`. Integral and expression kernels call
their sub-kernels before anything that reads what they build:

- A proxy coefficient takes the place of an interpolation of an expression. Its
  sub-kernel evaluates the expression at the interpolation points of the target
  element, and the degrees of freedom are the interpolation matrix of the
  target element applied to the values it writes.
- An interpolation of an expression that is only linear in an argument has its
  element table built for each cell, by contracting the values its sub-kernel
  writes for each degree of freedom of the argument with the target element
  table.
"""

from __future__ import annotations

import typing

import basix.ufl
import numpy as np
import ufl

import ffcx.codegeneration.lnodes as L
from ffcx.analysis import ProxyCoefficient

if typing.TYPE_CHECKING:
    from ffcx.codegeneration.backend import FFCXBackend
    from ffcx.codegeneration.C.external import CExternalOperator
    from ffcx.ir.elementtables import InterpolatedTableT
    from ffcx.ir.representation import SubKernelIR


def generate_sub_kernels(
    sub_kernels: list[SubKernelIR],
    backend: FFCXBackend,
    interpolated_tables: dict[str, InterpolatedTableT],
) -> list[L.LNode]:
    """Generate the calls to the sub-kernels of a kernel, and what they build.

    Args:
        sub_kernels: The sub-kernels the kernel calls.
        backend: The backend of the kernel.
        interpolated_tables: The element tables built for each cell, by name.

    Returns:
        Code to place before anything that reads a proxy coefficient or an
        element table built for each cell.
    """
    proxies: list[SubKernelIR] = []
    tables: list[SubKernelIR] = []
    for kernel in sub_kernels:
        if isinstance(kernel.target, ProxyCoefficient):
            proxies.append(kernel)
        else:
            tables.append(kernel)
    return _proxy_coefficients(proxies, backend) + _interpolated_tables(
        tables, backend, interpolated_tables
    )


def _call_kernel(
    kernel: SubKernelIR, label: str, backend: FFCXBackend, side: int = 0
) -> tuple[L.Symbol, list[L.Declaration], list[L.LNode]]:
    """Call a sub-kernel.

    The coefficients and constants that the kernel reads are packed into arrays
    of their own, in the order that the kernel reads them.

    Args:
        kernel: The kernel to call.
        label: Prefix of the names of the arrays declared.
        backend: The backend of the calling kernel.
        side: The cell to evaluate on, 0 or, on an interior facet, 1 for the
            cell restricted to '-'.

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
    # A kernel adds to the tensor it writes
    statements.append(L.ForRange(index, 0, size, [L.Assign(values[index], 0.0)]))
    statements.append(
        L.Statement(
            L.CallOp(
                kernel.name,
                (values, coefficients, constants, *_kernel_arguments(kernel, backend, side)),
            )
        )
    )
    return values, declarations, statements


def _kernel_arguments(kernel: SubKernelIR, backend: FFCXBackend, side: int) -> list[L.LExpr]:
    """The arguments a sub-kernel shares with the kernel that calls it, for one side."""
    symbols = backend.symbols
    return [
        _shift(symbols.coordinate_dofs, side * _num_coordinate_dofs(kernel, backend)),
        _shift(symbols.entity_local_index, side),
        _shift(symbols.quadrature_permutation, side),
        L.Symbol("custom_data", dtype=L.DataType.SCALAR),
    ]


def _call_external(
    kernel: SubKernelIR,
    operands: L.Symbol,
    label: str,
    backend: FFCXBackend,
    side: int,
    declarations: list[L.Declaration],
) -> tuple[L.Symbol, list[L.LNode]]:
    """Call the function of an external operator on its operands at the points.

    Args:
        kernel: The sub-kernel that evaluated the operands.
        operands: The array the sub-kernel wrote the operands to.
        label: Prefix of the names of the arrays declared.
        backend: The backend of the calling kernel.
        side: See `_call_kernel`.
        declarations: Declarations to add to.

    Returns:
        The array the function writes the values of the operator to,
        ``[point][component]``, and the statements that fill it.
    """
    proxy = kernel.target
    assert isinstance(proxy, ProxyCoefficient)
    # Checked by the C backend before the kernels are generated, see
    # `ffcx.codegeneration.C.external` (not imported at runtime: circular)
    external = typing.cast("CExternalOperator", proxy.external)
    element = proxy.ufl_element()
    # Output of the external function: the operator's values at the points,
    # [point][component]
    values = L.Symbol(f"{label}_external", dtype=L.DataType.SCALAR)
    size = kernel.shape[0] * element.reference_value_size
    declarations.append(L.ArrayDecl(values, sizes=size))
    call = L.CallOp(
        external.c_name,
        (
            values,
            operands,
            L.Symbol("NULL", dtype=L.DataType.SCALAR),
            *_kernel_arguments(kernel, backend, side),
        ),
        generated=False,
    )
    return values, [L.Statement(call)]


def _shift(pointer: L.Symbol, offset: int) -> L.LExpr:
    """The pointer advanced by ``offset`` entries, to the data of the second cell."""
    return pointer if offset == 0 else L.Add(pointer, offset)


def _num_coordinate_dofs(kernel: SubKernelIR, backend: FFCXBackend) -> int:
    """The number of entries of `coordinate_dofs` for one cell, three per node.

    The layout is that of the calling kernel, whose integration domain differs
    from the domain of the target space when that is a submesh of facets.
    """
    element = backend.access.integration_domain_coordinate_element
    if element is None:
        domain = ufl.domain.extract_unique_domain(kernel.target)
        assert isinstance(domain, ufl.Mesh)
        element = domain.ufl_coordinate_element()
    assert isinstance(element, basix.ufl._ElementBase)
    return 3 * int(element.basix_element.dim)


def _proxy_coefficients(kernels: list[SubKernelIR], backend: FFCXBackend) -> list[L.LNode]:
    """Generate the degrees of freedom of the proxy coefficients.

    They are written where they are read from, see
    `ffcx.codegeneration.symbols.FFCXBackendSymbols.coefficient_dof_access`.
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
        proxy = kernel.target
        assert isinstance(proxy, ProxyCoefficient)
        for side in range(sides):
            label = f"proxy{i}" if side == 0 else f"proxy{i}_{side}"
            values, declarations, statements = _call_kernel(kernel, label, backend, side)
            if proxy.external is not None:
                values, call = _call_external(kernel, values, label, backend, side, declarations)
                statements += call
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


def _interpolated_tables(
    kernels: list[SubKernelIR],
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
    needed = {table.interpolation for table in interpolated_tables.values()}
    parts: list[L.LNode] = []
    evaluated = {}
    for i, kernel in enumerate(kernels):
        if kernel.target not in needed:
            continue
        values, declarations, statements = _call_kernel(kernel, f"interpolation{i}", backend)
        evaluated[kernel.target] = (values, kernel.shape)
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
        values, (_, value_size, *dof_dims) = evaluated[data.interpolation]
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
