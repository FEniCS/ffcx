# Copyright (C) 2020-2024 Michal Habera, Chris Richardson and Garth N. Wells
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Utilities."""

import numpy as np
import numpy.typing as npt

try:
    import numba
except ImportError:
    numba = None


def dtype_to_c_type(dtype: npt.DTypeLike | str) -> str:
    """For a NumPy dtype, return the corresponding C type.

    Args:
        dtype: Numpy data type,

    Returns:
        Corresponding C type.
    """
    # Note: Possible aliases, e.g. numpy.longdouble, should test against char ID
    if np.dtype(dtype).char == "g":
        return "long double"
    if np.dtype(dtype) == np.intc:
        return "int"
    elif np.dtype(dtype).char == "f":
        return "float"
    elif np.dtype(dtype).char == "d":
        return "double"
    elif np.dtype(dtype) == np.complex64:
        return "float _Complex"
    elif np.dtype(dtype) == np.complex128:
        return "double _Complex"
    else:
        raise RuntimeError(f"Unknown NumPy type for: {dtype}")


def dtype_to_scalar_dtype(dtype: npt.DTypeLike | str) -> np.dtype:
    """For a NumPy dtype, return the corresponding real dtype.

    Args:
        dtype: Numpy data type

    Returns:
        ``numpy.dtype`` for the real component of ``dtype``.
    """
    if np.issubdtype(dtype, np.floating):
        return np.dtype(dtype)
    elif np.issubdtype(dtype, np.complexfloating):
        return np.dtype(dtype).type(0).real.dtype
    elif np.issubdtype(dtype, np.integer):
        return np.dtype(dtype)
    else:
        raise RuntimeError(f"Cannot get value dtype for '{dtype}'. ")


def numba_ufcx_kernel_signature(dtype: npt.DTypeLike, xdtype: npt.DTypeLike):
    """Return a Numba C signature for the UFCx ``tabulate_tensor`` interface.

    Args:
        dtype: The scalar type for the finite element data.
        xdtype: The geometry float type.

    Returns:
        A Numba signature (``numba.core.typing.templates.Signature``).

    Raises:
        ImportError: If ``numba`` cannot be imported.
    """
    try:
        import numba.types as types
        from numba import from_dtype

        return types.void(
            types.CPointer(from_dtype(dtype)),
            types.CPointer(from_dtype(dtype)),
            types.CPointer(from_dtype(dtype)),
            types.CPointer(from_dtype(xdtype)),
            types.CPointer(types.intc),
            types.CPointer(types.uint8),
            types.CPointer(types.void),
        )
    except ImportError as e:
        raise e


if numba is not None:

    @numba.extending.intrinsic
    def empty_void_pointer(typingctx):
        """Custom intrinsic to return an empty void* pointer.

        This function creates a void pointer initialized to null (0).
        This is used to pass a nullptr to the UFCx tabulate_tensor interface.

        Args:
            typingctx: The typing context.

        Returns:
            A Numba signature and a code generation function that returns a void pointer.
        """

        def codegen(context, builder, signature, args):
            null_ptr = context.get_constant(numba.types.voidptr, 0)
            return null_ptr

        sig = numba.types.voidptr()
        return sig, codegen

    @numba.extending.intrinsic
    def get_void_pointer(typingctx, arr):
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

    def create_voidptr_to_dtype_ptr_caster(target_dtype):
        """Factory that creates a Numba intrinsic casting void* to CPointer(target_dtype).

        The produced intrinsic accepts either `CPointer(void)` or `voidptr` as input
        (matching how UFCx kernels receive `custom_data`) and returns a
        `CPointer(target_dtype)` for convenient indexed access in Numba cfuncs.

        Args:
            target_dtype: A Numba scalar type (e.g. `numba.types.float64`).

        Returns:
            A Numba intrinsic function that performs the cast.
        """

        @numba.extending.intrinsic
        def voidptr_to_dtype_ptr(typingctx, src):
            # Accept void pointers in various Numba representations:
            # - CPointer(void): from UFCx cfunc signatures (shows as 'none*')
            # - voidptr: from numba.cfunc("...(voidptr)") signatures
            is_cpointer_void = (
                isinstance(src, numba.types.CPointer) and src.dtype == numba.types.void
            )
            is_voidptr = src == numba.types.voidptr

            if is_cpointer_void or is_voidptr:
                result_type = numba.types.CPointer(target_dtype)
                sig = result_type(src)

                def codegen(context, builder, signature, args):
                    [src_val] = args
                    dst_type = context.get_value_type(result_type)
                    return builder.bitcast(src_val, dst_type)

                return sig, codegen

            # Raise a clear error if the source type is not a void pointer
            msg = (
                "voidptr_to_dtype_ptr expects a void pointer (CPointer(void) or voidptr), "
                f"got {src}. Ensure you are passing a void* (e.g., custom_data from UFCx kernel signature)."
            )
            raise numba.core.errors.TypingError(msg)

        return voidptr_to_dtype_ptr

    # Pre-built casters for common types (convenience wrappers)
    voidptr_to_float64_ptr = create_voidptr_to_dtype_ptr_caster(numba.types.float64)
    voidptr_to_float32_ptr = create_voidptr_to_dtype_ptr_caster(numba.types.float32)
    voidptr_to_int32_ptr = create_voidptr_to_dtype_ptr_caster(numba.types.int32)
    voidptr_to_int64_ptr = create_voidptr_to_dtype_ptr_caster(numba.types.int64)
