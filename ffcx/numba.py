# Copyright (C) 2020-2024 Michal Habera, Chris Richardson and Garth N. Wells
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Numba specific functionality."""

from collections.abc import Callable
from typing import Any

import numba
import numpy.typing as npt


def numba_ufcx_kernel_signature(dtype: npt.DTypeLike, xdtype: npt.DTypeLike) -> Any:
    """Return a Numba C signature for the UFCx ``tabulate_tensor`` interface.

    Args:
        dtype: The scalar type for the finite element data.
        xdtype: The geometry float type.

    Returns:
        A Numba signature (``numba.core.typing.templates.Signature``).
    """
    return numba.types.void(
        numba.types.CPointer(numba.from_dtype(dtype)),
        numba.types.CPointer(numba.from_dtype(dtype)),
        numba.types.CPointer(numba.from_dtype(dtype)),
        numba.types.CPointer(numba.from_dtype(xdtype)),
        numba.types.CPointer(numba.types.intc),
        numba.types.CPointer(numba.types.uint8),
        numba.types.CPointer(numba.types.void),
    )


@numba.extending.intrinsic
def empty_void_pointer(
    typingctx: Any,
) -> tuple[Any, Callable[[Any, Any, Any, Any], Any]]:
    """Custom intrinsic to return an empty void* pointer.

    This function creates a void pointer initialized to null (0).
    This is used to pass a nullptr to the UFCx tabulate_tensor interface.

    Args:
        typingctx: The typing context.

    Returns:
        A Numba signature and a code generation function that returns a void pointer.
    """

    def codegen(context: Any, builder: Any, signature: Any, args: Any) -> Any:
        null_ptr = context.get_constant(numba.types.voidptr, 0)
        return null_ptr

    sig = numba.types.voidptr()
    return sig, codegen


@numba.extending.intrinsic
def get_void_pointer(typingctx: Any, arr: npt.NDArray) -> tuple[numba.types.RawPointer, Callable]:
    """Custom intrinsic to get a void* pointer from a NumPy array.

    This function takes a NumPy array and returns a void pointer to the array's data.
    This is used to pass custom data organised in a NumPy array
    to the UFCx tabulate_tensor interface.

    Args:
        typingctx: The typing context.
        arr: The NumPy array to get the void pointer to the first element from.
        In a multi-dimensional NumPy array, the memory is laid out in a contiguous
        block of memory, see
        https://numpy.org/doc/stable/reference/arrays.ndarray.html#internal-memory-layout-of-an-ndarray

    Returns:
        sig: A Numba signature, which specifies the numba type (here voidptr),
        codegen: A code generation function, which returns the LLVM IR to cast
        the raw data pointer to the first element of the of the contiguous block of memory
        of the NumPy array to void*.
    """
    if not isinstance(arr, numba.types.Array):
        raise TypeError("Expected a NumPy array")

    def codegen(context, builder, signature, args):
        """Generate LLVM IR code to convert a NumPy array to a void* pointer.

        This function generates the necessary LLVM IR instructions to:
        1. Allocate memory for the array on the stack.
        2. Cast the allocated memory to a void* pointer.

        Args:
            context: The LLVM context.
            builder: The LLVM IR builder.
            signature: The function signature.
            args: The input arguments (NumPy array).

        Returns:
            A void* pointer to the array's data.
        """
        [arr] = args
        raw_ptr = numba.core.cgutils.alloca_once_value(builder, arr)
        void_ptr = builder.bitcast(raw_ptr, context.get_value_type(numba.types.voidptr))
        return void_ptr

    sig = numba.types.voidptr(arr)
    return sig, codegen
