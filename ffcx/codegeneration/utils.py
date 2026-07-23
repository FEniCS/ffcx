# Copyright (C) 2020-2024 Michal Habera, Chris Richardson and Garth N. Wells
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Utilities."""

import numpy as np
import numpy.typing as npt


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
