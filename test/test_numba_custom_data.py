# Test that the Numba voidptr -> typed pointer caster factory works in ffcx utils
import ctypes

import numpy as np
import pytest

import ffcx.codegeneration.utils as codegen_utils

# Skip the tests if Numba is not available in the environment.
numba = pytest.importorskip("numba")

float64_ptr_caster = codegen_utils._create_voidptr_to_dtype_ptr_caster(numba.types.float64)
int32_ptr_caster = codegen_utils._create_voidptr_to_dtype_ptr_caster(numba.types.int32)


def test_numba_voidptr_struct_like_mixed_types():
    """Test reading a struct-like mixed-type buffer: float64 + int32.

    We create a NumPy structured array with fields ('scale', float64) and
    ('id', int32) with padding to align to 16 bytes. The kernel casts the
    void* to float64* and int32* and reads the corresponding offsets.
    """
    sig = codegen_utils.numba_ufcx_kernel_signature(np.float64, np.float64)

    @numba.cfunc(sig, nopython=True)
    def tabulate(b_, w_, c_, coords_, local_index, orientation, custom_data):
        b = numba.carray(b_, (1,), dtype=np.float64)
        fptr = float64_ptr_caster(custom_data)
        iptr = int32_ptr_caster(custom_data)
        scale = fptr[0]
        # int32 index for offset 8 bytes == 8/4 == 2
        id_val = iptr[2]
        b[0] = scale + id_val

    b = np.zeros(1, dtype=np.float64)
    w = np.zeros(1, dtype=np.float64)
    c = np.zeros(1, dtype=np.float64)
    coords = np.zeros(9, dtype=np.float64)
    local_index = np.array([0], dtype=np.int32)
    orientation = np.array([0], dtype=np.uint8)

    # structured dtype with C-compatible alignment
    dtype = np.dtype([("scale", np.float64), ("id", np.int32)], align=True)
    arr = np.zeros(1, dtype=dtype)
    arr["scale"][0] = 1.25
    arr["id"][0] = 5

    ptr = arr.ctypes.data

    tabulate.ctypes(
        b.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        w.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        c.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        coords.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        local_index.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        orientation.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)),
        ctypes.c_void_p(ptr),
    )

    assert b[0] == pytest.approx(6.25)
