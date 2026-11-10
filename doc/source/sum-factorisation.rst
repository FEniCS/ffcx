Sum factorisation
=================

FFCx enables sum factorisation by default wherever the elements and
quadrature rule support it. The ``sum_factorization`` code-generation
option defaults to ``True``; applications do not need to set it explicitly.

For eligible cell integrals, FFCx evaluates coefficients, the spatial
coordinate and the cell Jacobian on the tensor-product quadrature rule, then
contracts the argument tables one direction at a time. This reduces the
work for a bilinear form of degree :math:`p` in dimension :math:`d` from
:math:`O(p^{3d})` to :math:`O(p^{2d+1})` when the number of quadrature points
per direction grows proportionally to :math:`p`.

Element and quadrature requirements
-----------------------------------

Directional contractions require cell integrals on quadrilaterals or
hexahedra, with argument elements whose basis tables admit a supported
tensor-product representation and a tensor-product quadrature rule.
FFCx constructs tensor-product rules for the supported quadrature schemes
on these cells.

FFCx preserves the supplied element's degree-of-freedom (DOF) ordering.
Use Basix's ``create_tp_element`` to create an element with tensor-product
DOF ordering, then wrap it for UFL:

.. code-block:: python

   import basix
   import basix.ufl

   element = basix.ufl.wrap_element(
       basix.create_tp_element(
           basix.ElementFamily.P,
           basix.CellType.hexahedron,
           degree=4,
           lagrange_variant=basix.LagrangeVariant.gll_warped,
       )
   )

Use ``basix.CellType.quadrilateral`` for a quad mesh. Vector and tensor
spaces can be built from this scalar element using
``basix.ufl.blocked_element``.

In DOLFINx, create the function space using this element and compile the
form normally. Given a hexahedral mesh ``mesh``:

.. code-block:: python

   import ufl
   from dolfinx import fem

   V = fem.functionspace(mesh, element)
   u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
   a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
   compiled_a = fem.form(a)

With the direct FFCx JIT interface, the same default applies:

.. code-block:: python

   from ffcx.codegeneration.jit import compile_forms

   compiled_forms, module, code = compile_forms([a])

Automatic fallback
------------------

FFCx uses ordinary quadrature when the integral cannot use directional
contractions. This includes other cell types, facet integrals, argument
tables without supported tensor factors, assembly of only the diagonal
(``part="diagonal"``), and custom or vertex quadrature rules without
tensor-product information. These cases do not require
disabling the option. A coefficient or geometry table without tensor
factors can still be evaluated normally inside a supported contraction.
FFCx also uses ordinary quadrature when a kernel's arrays would not fit the
stack limit described below.

Stack usage
-----------

A sum-factorised kernel evaluates the integrand, coefficients and geometry
at every point of the quadrature rule and keeps these values in arrays on
the stack. The space they take grows with the number of quadrature points,
so cubically with the degree on hexahedra, and with the number of terms in
the integrand. For example, an elasticity form on hexahedra takes about
126 KiB of stack per kernel call at degree 4 and about 510 KiB at degree 8.

The main thread of a program usually has 8 MiB of stack, but other threads
can have much less: 512 KiB for secondary threads on macOS, and 128 KiB for
threads on musl-based systems such as Alpine Linux. A kernel that exceeds
its thread's stack crashes the program.

The ``sum_factorization_stack_limit`` option bounds this space, in bytes.
It defaults to 256 KiB. When a kernel's arrays would exceed it, FFCx splits
the terms of the integrand into batches and evaluates and contracts one
batch at a time, so that each batch's arrays fit. This costs some speed:
the integrand is evaluated once per batch, and the element tensor is
updated once per batch. If even a single term does not fit, for example
with very large quadrature rules, FFCx uses ordinary quadrature for that
integral. Set the option to ``0`` for no limit.

With the FFCx JIT:

.. code-block:: python

   compiled_forms, module, code = compile_forms(
       [a], options={"sum_factorization_stack_limit": 96 * 1024}
   )

With DOLFINx:

.. code-block:: python

   compiled_a = fem.form(
       a, form_compiler_options={"sum_factorization_stack_limit": 96 * 1024}
   )

From the command line:

.. code-block:: console

   ffcx --sum_factorization_stack_limit 98304 forms.py

Lower the limit if kernels run on threads with small stacks: about 96 KiB is
safe for 128 KiB threads. Raise it, or set it to ``0``, if kernels only run on
threads with large stacks, such as the main thread, and the fastest kernels
at high degree matter more.

Disabling sum factorisation
---------------------------

To use ordinary quadrature explicitly with the FFCx JIT:

.. code-block:: python

   compiled_forms, module, code = compile_forms(
       [a], options={"sum_factorization": False}
   )

With DOLFINx, set the corresponding form compiler option:

.. code-block:: python

   compiled_a = fem.form(
       a, form_compiler_options={"sum_factorization": False}
   )

From the command line:

.. code-block:: console

   ffcx --no-sum_factorization forms.py

This option controls kernel generation. C compiler flags are configured
separately; see :doc:`compiler-flags` for suggested optimization and FMA
settings.
