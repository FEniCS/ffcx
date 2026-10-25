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
