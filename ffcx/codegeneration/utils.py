# Copyright (C) 2020-2024 Michal Habera, Chris Richardson and Garth N. Wells
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Utilities."""

import typing

import numpy as np
import numpy.typing as npt


def dtype_to_c_type(dtype: typing.Union[npt.DTypeLike, str]) -> str:
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


def dtype_to_scalar_dtype(dtype: typing.Union[npt.DTypeLike, str]) -> np.dtype:
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
        )
    except ImportError as e:
        raise e
