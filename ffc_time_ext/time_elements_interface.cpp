// Copyright (C) 2012 Benjamin Kehlet
//
// This file is part of FFC.
//
// FFC is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// FFC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with FFC.  If not, see <http://www.gnu.org/licenses/>.
//
// First added:  2012-08-20
// Last changed: 2013-07-11


#include <Python.h>

#include <numpy/arrayobject.h>

#include "time_elements.h"

static PyObject *compute_lobatto_interface(PyObject *dummy, PyObject *args)
{
  int degree;

  /* parse argument tuple */
  if (!PyArg_ParseTuple(args, "i", &degree))
  {
    return NULL; /* PyArg_ParseTuple has raised an exception */
  }

  npy_intp n = degree+1;

  PyArrayObject *py_array_points =  (PyArrayObject*) PyArray_SimpleNew(1, &n, NPY_DOUBLE);

  double *points = (double*)  PyArray_DATA(py_array_points);

  compute_lobatto_points(points, degree);

  // return values
  return (PyObject*) py_array_points;
}

static PyObject *compute_radau_interface(PyObject *dummy, PyObject *args)
{
  int degree;

  /* parse argument tuple */
  if (!PyArg_ParseTuple(args, "i", &degree))
  {
    return NULL; /* PyArg_ParseTuple has raised an exception */
  }

  npy_intp n = degree+1;

  PyArrayObject *py_array_points =  (PyArrayObject*) PyArray_SimpleNew(1, &n, NPY_DOUBLE);

  double *points = (double*)  PyArray_DATA(py_array_points);

  compute_radau_points(points, degree);

  // return values
  return (PyObject*) py_array_points;
}

static PyObject *compute_legendre_coeffs_interface(PyObject *dummy, PyObject *args)
{
  PyArrayObject *points_array;

  /* parse argument tuple */
  if (!PyArg_ParseTuple(args, "O!",
			&PyArray_Type, &points_array))
  {
    return NULL; /* PyArg_ParseTuple has raised an exception */
  }

  const npy_intp num_points = PyArray_DIMS(points_array)[0];
  npy_intp dims[2] = { num_points, num_points };
  PyArrayObject *py_array_coeffs =  (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_DOUBLE);

  double *coeffs = (double*)  PyArray_DATA(py_array_coeffs);

  compute_legendre_coeffs(coeffs, (double*) PyArray_DATA(points_array), num_points, num_points);

  // return values
  return (PyObject*) py_array_coeffs;
}

static char compute_lobatto_doc[] = \
  "Doc string for compute_lobatto_points";
static char compute_radau_doc[] = \
  "Doc string for compute_radau_points";
static char compute_legendre_coeffs_doc[] = \
  "Doc string for compute_legendre_coeffs";


static PyMethodDef time_elements_ext_methods[] = {
  {"compute_lobatto_points",
   compute_lobatto_interface,
   METH_VARARGS,
   compute_lobatto_doc},
  {"compute_radau_points",
   compute_radau_interface,
   METH_VARARGS,
   compute_radau_doc},
  {"compute_legendre_coeffs",
   compute_legendre_coeffs_interface,
   METH_VARARGS,
   compute_legendre_coeffs_doc},

  {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "time_elements_ext",             /* m_name */
        "This is a module..",            /* m_doc */
        -1,                              /* m_size */
        time_elements_ext_methods,       /* m_methods */
        NULL,                            /* m_reload */
        NULL,                            /* m_traverse */
        NULL,                            /* m_clear */
        NULL,                            /* m_free */
    };

    PyMODINIT_FUNC
    PyInit_time_elements_ext(void)
    {
        PyObject *m;
        m = PyModule_Create(&moduledef);
        import_array();
        if(m==NULL)
           return NULL;
        return m;
    }
#else
    PyMODINIT_FUNC
    inittime_elements_ext(void)
    {
        (void)Py_InitModule("time_elements_ext", time_elements_ext_methods);
        import_array();
    }
#endif

