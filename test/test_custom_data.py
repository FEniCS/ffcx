# Copyright (C) 2024 Susanne Claus
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numpy as np
import pytest


def test_tabulate_tensor_integral_add_values():
    pytest.importorskip("cffi")

    from cffi import FFI

    # Define custom tabulate tensor function in C with a struct
    # Step 1: Define the function in C and set up the CFFI builder
    ffibuilder = FFI()
    ffibuilder.set_source(
        "_cffi_kernelA",
        r"""
        typedef struct {
            size_t size;
            double* values;
        } cell_data;

        void tabulate_tensor_integral_add_values(double* restrict A,
                                                const double* restrict w,
                                                const double* restrict c,
                                                const double* restrict coordinate_dofs,
                                                const int* restrict entity_local_index,
                                                const uint8_t* restrict quadrature_permutation,
                                                void* custom_data)
        {
            // Cast the void* custom_data to cell_data*
            cell_data* custom_data_ptr = (cell_data*)custom_data;

            // Access the custom data
            size_t size = custom_data_ptr->size;
            double* values = custom_data_ptr->values;

            // Use the values in your computations
            for (size_t i = 0; i < size; i++) {
                A[0] += values[i];
            }
        }
        """,
    )
    ffibuilder.cdef(
        """
        typedef struct {
            size_t size;
            double* values;
        } cell_data;

        void tabulate_tensor_integral_add_values(double* restrict A,
                                                const double* restrict w,
                                                const double* restrict c,
                                                const double* restrict coordinate_dofs,
                                                const int* restrict entity_local_index,
                                                const uint8_t* restrict quadrature_permutation,
                                                void* custom_data);
        """
    )

    # Step 2: Compile the C code
    ffibuilder.compile(verbose=True)

    # Step 3: Import the compiled library
    from _cffi_kernelA import ffi, lib

    # Define cell data
    values = np.array([2.0, 1.0], dtype=np.float64)
    size = len(values)
    expected_result = np.array([3.0], dtype=np.float64)

    # Define the input arguments
    A = np.zeros(1, dtype=np.float64)
    w = np.array([1.0], dtype=np.float64)
    c = np.array([0.0], dtype=np.float64)
    coordinate_dofs = np.array(
        [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0], dtype=np.float64
    )
    entity_local_index = np.array([0], dtype=np.int32)
    quadrature_permutation = np.array([0], dtype=np.uint8)

    # Cast the arguments to the appropriate C types
    A_ptr = ffi.cast("double*", A.ctypes.data)
    w_ptr = ffi.cast("double*", w.ctypes.data)
    c_ptr = ffi.cast("double*", c.ctypes.data)
    coordinate_dofs_ptr = ffi.cast("double*", coordinate_dofs.ctypes.data)
    entity_local_index_ptr = ffi.cast("int*", entity_local_index.ctypes.data)
    quadrature_permutation_ptr = ffi.cast("uint8_t*", quadrature_permutation.ctypes.data)

    # Use ffi.from_buffer to create a CFFI pointer from the NumPy array
    values_ptr = ffi.cast("double*", values.ctypes.data)

    # Allocate memory for the struct
    custom_data = ffi.new("cell_data*")
    custom_data.size = size
    custom_data.values = values_ptr

    # Cast the struct to void*
    custom_data_ptr = ffi.cast("void*", custom_data)

    # Call the function
    lib.tabulate_tensor_integral_add_values(
        A_ptr,
        w_ptr,
        c_ptr,
        coordinate_dofs_ptr,
        entity_local_index_ptr,
        quadrature_permutation_ptr,
        custom_data_ptr,
    )

    # Assert the result
    np.testing.assert_allclose(A, expected_result, rtol=1e-5)
