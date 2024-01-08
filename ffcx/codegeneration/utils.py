# Copyright (C) 2020-2023 Michal Habera and Chris Richardson
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import typing

import numpy as _np
import numpy.typing as _npt


# def cdtype_to_numpy(cdtype: str) -> str:
#     """Map a C data type string NumPy datatype string."""
#     if cdtype == "double":
#         return "float64"
#     elif cdtype == "double _Complex":
#         return "complex128"
#     elif cdtype == "float":
#         return "float32"
#     elif cdtype == "float _Complex":
#         return "complex64"
#     elif cdtype == "long double":
#         return "longdouble"
#     else:
#         raise RuntimeError(f"Unknown NumPy type for: {cdtype}")


def dtype_to_c_type(dtype: typing.Union[_npt.DTypeLike, str]) -> str:
    """For a NumPy dtype, return the corresponding C type.

    Args:
        dtype: Numpy data type

    Returns:
        Corresponding C type
    """
    if _np.dtype(dtype) == _np.float32:
        return "float"
    elif _np.dtype(dtype) == _np.float64:
        return "double"
    elif _np.dtype(dtype) == _np.complex64:
        return "float _Complex"
    elif _np.dtype(dtype) == _np.complex128:
        return "double _Complex"
    elif _np.dtype(dtype) == _np.intc:
        return "int"
    else:
        raise RuntimeError(f"Unknown NumPy type for: {dtype}")


def dtype_to_scalar_dtype(dtype: typing.Union[_npt.DTypeLike, str]) -> str:
    """For a NumPy dtype, return the corresponding real dtype.

    Args:
        dtype: Numpy data type

    Returns:
        ``numpy.dtype`` for the real component of ``dtype``.
    """
    _dtype = _np.dtype(dtype)
    return _np.real(_dtype.type(0)).dtype


def scalar_to_value_type(scalar_type: str) -> str:
    """The C value type associated with a C scalar type.

    Args:
      scalar_type: A C type.

    Returns:
      The value type associated with ``scalar_type``. E.g., if
      ``scalar_type`` is ``float _Complex`` the return value is 'float'.

    """
    return scalar_type.replace(' _Complex', '')
