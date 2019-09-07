import pytest
import numpy as np
import cffi

import ufl
import ffc.codegeneration.jit


def float_to_type(name):
    """Map a string name to C and NumPy types"""
    if name == "double":
        return "double", np.float64
    elif name == "double complex":
        return "double _Complex", np.complex128
    elif name == "float":
        return "float", np.float32
    elif name == "float complex":
        return "float _Complex", np.complex64
    elif name == "long double":
        return "long double", np.longdouble
    else:
        raise RuntimeError("Unknown C type for: {}".format(name))


def test_expression():
    e = ufl.VectorElement("P", "triangle", 1)

    f = ufl.Coefficient(e)

    a_mat = np.array([[1.0, 2.0], [1.0, 1.0]])
    a = ufl.as_matrix(a_mat)
    expr = ufl.dot(a, f)

    points = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    obj, module = ffc.codegeneration.jit.compile_expressions([(expr, points)])

    ffi = cffi.FFI()
    kernel = obj[0][0]

    c_type, np_type = float_to_type("double")

    A = np.zeros((2, 3), dtype=np_type)
    f_mat = np.array([[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]])
    # Coefficient storage XXXYYY
    w = np.array(f_mat.flatten(), dtype=np_type)
    c = np.array([], dtype=np_type)

    coords = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0], dtype=np.float64)
    kernel.tabulate_expression(
        ffi.cast('{type} *'.format(type=c_type), A.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), w.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), c.ctypes.data),
        ffi.cast('double *', coords.ctypes.data))

    assert np.allclose(A, np.dot(a_mat, f_mat))
