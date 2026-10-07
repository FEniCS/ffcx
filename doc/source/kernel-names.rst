Generated kernel names
======================

FFCx normally identifies generated integral kernels using a hash of the UFL
form and compiler state. When compiling a Python UFL file with the ``ffcx``
command, FFCx additionally emits a stable kernel function name derived from
the Python variable containing the form. For example,

.. code-block:: python

   mass = inner(u, v) * dx

produces a function named using the output namespace, form name, integral
type, integral-group index and cell type. For an input file named ``forms.py``
on a triangular mesh, the name is

.. code-block:: c

   tabulate_tensor_forms_mass_cell_0_triangle

The generated header declares this function. The hash-named function remains
available as a forwarding symbol for compatibility.

Explicit names
--------------

Set ``ffcx_kernel_name`` integral metadata to specify the complete function
name:

.. code-block:: python

   mass = inner(u, v) * dx(metadata={"ffcx_kernel_name": "tabulate_tensor_p1_mass"})

This produces ``tabulate_tensor_p1_mass`` exactly: FFCx does not add the output
namespace or cell type to an explicit name. Explicit names must contain only
letters, digits and underscores, and must not start with a digit. Names must be
unique within a generated output file. Integrals that FFCx combines into one
kernel must use the same explicit name.

Compilation through the Python API without an object name or explicit
``ffcx_kernel_name`` continues to use only the hash-derived kernel name.
