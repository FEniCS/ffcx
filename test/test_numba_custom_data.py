# Test that the Numba voidptr -> typed pointer caster works in ffcx utils
import ctypes
import numpy as np
import pytest

numba = pytest.importorskip("numba")

from ffcx.codegeneration.utils import (
    numba_ufcx_kernel_signature,
    voidptr_to_float64_ptr,
    voidptr_to_int32_ptr,
)


def test_numba_voidptr_caster_basic():
    """Simple test: Numba cfunc reads a double from custom_data via the caster."""
    sig = numba_ufcx_kernel_signature(np.float64, np.float64)

    @numba.cfunc(sig, nopython=True)
    def tabulate(b_, w_, c_, coords_, local_index, orientation, custom_data):
        b = numba.carray(b_, (1,), dtype=np.float64)
        # Cast void* to float64*
        typed = voidptr_to_float64_ptr(custom_data)
        b[0] = typed[0]

    # Prepare arguments
    b = np.zeros(1, dtype=np.float64)
    w = np.zeros(1, dtype=np.float64)
    c = np.zeros(1, dtype=np.float64)
    coords = np.zeros(9, dtype=np.float64)
    local_index = np.array([0], dtype=np.int32)
    orientation = np.array([0], dtype=np.uint8)

    # custom_data: single double value
    val = np.array([2.5], dtype=np.float64)
    val_ptr = val.ctypes.data

    # Call the compiled cfunc via ctypes
    tabulate.ctypes(
        b.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        w.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        c.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        coords.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        local_index.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        orientation.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)),
        ctypes.c_void_p(val_ptr),
    )

    assert b[0] == pytest.approx(2.5)


def test_numba_voidptr_caster_int32():
    """Test casting void* to int32* and reading an integer value."""
    sig = numba_ufcx_kernel_signature(np.float64, np.float64)

    @numba.cfunc(sig, nopython=True)
    def tabulate(b_, w_, c_, coords_, local_index, orientation, custom_data):
        b = numba.carray(b_, (1,), dtype=np.float64)
        typed = voidptr_to_int32_ptr(custom_data)
        # Promote int32 to float64 for the output
        b[0] = typed[0]

    b = np.zeros(1, dtype=np.float64)
    w = np.zeros(1, dtype=np.float64)
    c = np.zeros(1, dtype=np.float64)
    coords = np.zeros(9, dtype=np.float64)
    local_index = np.array([0], dtype=np.int32)
    orientation = np.array([0], dtype=np.uint8)

    val = np.array([7], dtype=np.int32)
    val_ptr = val.ctypes.data

    tabulate.ctypes(
        b.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        w.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        c.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        coords.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        local_index.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        orientation.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)),
        ctypes.c_void_p(val_ptr),
    )

    assert b[0] == pytest.approx(7.0)


def test_numba_voidptr_caster_multiple_params():
    """Test reading multiple float64 parameters from custom_data."""
    sig = numba_ufcx_kernel_signature(np.float64, np.float64)

    @numba.cfunc(sig, nopython=True)
    def tabulate(b_, w_, c_, coords_, local_index, orientation, custom_data):
        b = numba.carray(b_, (1,), dtype=np.float64)
        typed = voidptr_to_float64_ptr(custom_data)
        b[0] = typed[0] + typed[1] + typed[2]

    b = np.zeros(1, dtype=np.float64)
    w = np.zeros(1, dtype=np.float64)
    c = np.zeros(1, dtype=np.float64)
    coords = np.zeros(9, dtype=np.float64)
    local_index = np.array([0], dtype=np.int32)
    orientation = np.array([0], dtype=np.uint8)

    vals = np.array([1.5, 2.0, 3.0], dtype=np.float64)
    vals_ptr = vals.ctypes.data

    tabulate.ctypes(
        b.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        w.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        c.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        coords.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        local_index.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        orientation.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)),
        ctypes.c_void_p(vals_ptr),
    )

    assert b[0] == pytest.approx(6.5)


def test_numba_voidptr_struct_like_mixed_types():
    """Test reading a struct-like mixed-type buffer: float64 + int32.

    We create a NumPy structured array with fields ('scale', float64) and
    ('id', int32) with padding to align to 16 bytes. The kernel casts the
    void* to float64* and int32* and reads the corresponding offsets.
    """
    sig = numba_ufcx_kernel_signature(np.float64, np.float64)

    @numba.cfunc(sig, nopython=True)
    def tabulate(b_, w_, c_, coords_, local_index, orientation, custom_data):
        b = numba.carray(b_, (1,), dtype=np.float64)
        fptr = voidptr_to_float64_ptr(custom_data)
        iptr = voidptr_to_int32_ptr(custom_data)
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

    # structured dtype: float64 at offset 0, int32 at offset 8, pad int32
    dtype = np.dtype([("scale", np.float64), ("id", np.int32), ("pad", np.int32)])
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
