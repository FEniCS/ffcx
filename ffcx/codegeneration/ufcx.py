# Copyright (C) 2024 Jack S. Hale
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later


def create_numba_tabulate_tensor_signature(scalar_type, real_type):
    import numba
    
    c_signature = numba.types.void(
        numba.types.CPointer(numba.typeof(scalar_type)),
        numba.types.CPointer(numba.typeof(scalar_type)),
        numba.types.CPointer(numba.typeof(scalar_type)),
        numba.types.CPointer(numba.typeof(real_type)),
        numba.types.CPointer(numba.types.int32),
        numba.types.CPointer(numba.types.int32))

    return c_signature
