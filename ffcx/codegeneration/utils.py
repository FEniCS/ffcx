# Copyright (C) 2020-2023 Michal Habera and Chris Richardson
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numpy.typing as _npt


def cdtype_to_numpy(cdtype: str):
    """Map a C data type string NumPy datatype string."""
    if cdtype == "double":
        return "float64"
    elif cdtype == "double _Complex":
        return "complex128"
    elif cdtype == "float":
        return "float32"
    elif cdtype == "float _Complex":
        return "complex64"
    elif cdtype == "long double":
        return "longdouble"
    else:
        raise RuntimeError(f"Unknown NumPy type for: {cdtype}")


def scalar_to_value_type(scalar_type: str) -> str:
    """The C value type associated with a C scalar type.

    Args:
      scalar_type: A C type.

    Returns:
      The value type associated with ``scalar_type``. E.g., if
      ``scalar_type`` is ``float _Complex`` the return value is 'float'.

    """
    return scalar_type.replace(' _Complex', '')


def numba_ufcx_kernel_signature(dtype: _npt.DTypeLike, xdtype: _npt.DTypeLike):
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
        from numba import from_dtype
        import numba.types as types
        return types.void(types.CPointer(from_dtype(dtype)), types.CPointer(from_dtype(dtype)),
                          types.CPointer(from_dtype(dtype)), types.CPointer(from_dtype(xdtype)),
                          types.CPointer(types.intc), types.CPointer(types.uint8))
    except ImportError as e:
        raise e
