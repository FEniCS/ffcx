# Copyright (C) 2020-2023 Michal Habera and Chris Richardson
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import typing

import numpy as _np
import numpy.typing as _npt


def dtype_to_c_type(dtype: typing.Union[_npt.DTypeLike, str]) -> str:
    """For a NumPy dtype, return the corresponding C type.

    Args:
        dtype: Numpy data type

    Returns:
        Corresponding C type
    """
    # Note: Possible aliases, e.g. numpy.longdouble, should test against char ID
    if _np.dtype(dtype).char == 'g':
        return "long double"
    if _np.dtype(dtype) == _np.intc:
        return "int"
    elif _np.dtype(dtype).char == 'f':
        return "float"
    elif _np.dtype(dtype).char == 'd':
        return "double"
    elif _np.dtype(dtype) == _np.complex64:
        return "float _Complex"
    elif _np.dtype(dtype) == _np.complex128:
        return "double _Complex"
    else:
        raise RuntimeError(f"Unknown NumPy type for: {dtype}")


def dtype_to_scalar_dtype(dtype: typing.Union[_npt.DTypeLike, str]) -> _npt.DTypeLike:
    """For a NumPy dtype, return the corresponding real dtype.

    Args:
        dtype: Numpy data type

    Returns:
        ``numpy.dtype`` for the real component of ``dtype``.
    """
    if _np.issubdtype(dtype, _np.floating):
        return dtype
    elif _np.issubdtype(dtype, _np.complexfloating):
        return _np.dtype(dtype).type(0).real.dtype
    else:
        raise RuntimeError("Cannot get value dtype. ")
