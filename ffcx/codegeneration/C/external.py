# Copyright (C) 2026 Jørgen S. Dokken and Paul T. Kühner
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""External operators evaluated by C functions.

FFCx evaluates a `ufl.ExternalOperator` in the kernel: it evaluates the
operands at the interpolation points of the element of the operator's space,
calls a function that maps those values onto the values of the operator, and
applies the interpolation matrix of the element to get its degrees of freedom.
The operator is then read in the kernel as a function of that space, built for
the cell.

FFCx provides no external operator of its own. The C backend compiles any
subclass of `ufl.ExternalOperator` that provides its C function, see
`CExternalOperator`.
"""

from __future__ import annotations

import hashlib
import typing

import numpy.typing as npt
import ufl

from ffcx.codegeneration.utils import dtype_to_c_type, dtype_to_scalar_dtype


@typing.runtime_checkable
class CExternalOperator(typing.Protocol):
    """An external operator evaluated by a C function, as the C backend requires.

    The function has the signature of a ``tabulate_tensor``, with ``T`` the
    scalar type and ``U`` the geometry type of the generated code::

        void c_name(T* restrict values, const T* restrict operands,
                    const T* restrict c, const U* restrict coordinate_dofs,
                    const int* restrict entity_local_index,
                    const uint8_t* restrict quadrature_permutation,
                    void* custom_data);

    ``operands`` holds the operands at the points, ``[point][component]``, with
    the components of the operands one after the other and each flattened in
    row-major order. ``values`` receives the values of the operator at the
    points, ``[point][component]``. ``c`` is not used and is NULL. The remaining
    arguments are those of the kernel that calls the function.

    As UFL rebuilds the operator when it processes a form, the subclass must
    carry these attributes through `_ufl_expr_reconstruct_`, and include them
    in `__eq__` and `__hash__`, or UFL merges operators evaluated by different
    functions.
    """

    c_name: str
    """The name of the function."""

    c_code: str | None
    """The definition of the function, which is placed in the generated file.

    If None, the function is declared, and must be provided when the generated
    code is linked.
    """


def signature(operator: ufl.ExternalOperator) -> str:
    """The C function of an external operator, for the signature of the generated code.

    Empty if the operator provides no C function.
    """
    if not isinstance(operator, CExternalOperator):
        return ""
    code = hashlib.sha1((operator.c_code or "").encode("utf-8")).hexdigest()
    return f"{operator.c_name};{code}"


declaration_template = """
void {name}({scalar_type}* restrict values,
            const {scalar_type}* restrict operands,
            const {scalar_type}* restrict c,
            const {geom_type}* restrict coordinate_dofs,
            const int* restrict entity_local_index,
            const uint8_t* restrict quadrature_permutation,
            void* custom_data);
"""


def generator(
    externals: list[ufl.ExternalOperator], options: dict[str, int | float | npt.DTypeLike]
) -> list[tuple[str, str]]:
    """Generate the definitions of the functions of external operators, or their declarations.

    Args:
        externals: The external operators.
        options: The options of the generated code.

    Returns:
        The code for the header, and for the source file, of each function.

    Raises:
        NotImplementedError: If an external operator does not provide a C function.
        ValueError: If external operators name the same function but give
            different code for it.
    """
    functions: dict[str, str | None] = {}
    for external in externals:
        if not isinstance(external, CExternalOperator):
            raise NotImplementedError(
                f"The C backend evaluates an external operator by a C function, which "
                f"{type(external).__name__} does not provide: see "
                f"{CExternalOperator.__module__}.{CExternalOperator.__name__}."
            )
        code = functions.setdefault(external.c_name, external.c_code)
        if code != external.c_code:
            raise ValueError(
                f"External operators name the function {external.c_name}, with different code."
            )

    scalar_type = options["scalar_type"]
    blocks = []
    for name, code in functions.items():
        if code is None:
            # Provided when the generated code is linked
            declaration = declaration_template.format(
                name=name,
                scalar_type=dtype_to_c_type(scalar_type),  # type: ignore[arg-type]
                geom_type=dtype_to_c_type(dtype_to_scalar_dtype(scalar_type)),  # type: ignore[arg-type]
            )
            blocks.append(("", f"// Function of an external operator{declaration}"))
        else:
            blocks.append(("", f"\n// Function of an external operator\n{code}\n"))
    return blocks
