Compiler flags
==============

FFCx generates C kernels. The C compiler flags used to build those kernels
can have a substantial effect on performance. Set these flags through the
JIT compilation interface, or through your build system when compiling
generated C files.

Suggested options
-----------------

For GCC and Clang, start with ``-O2 -ffp-contract=fast``. Benchmark the
following configurations on representative forms and polynomial degrees:

.. list-table:: Suggested GCC and Clang configurations
   :header-rows: 1
   :widths: 25 40 35

   * - Configuration
     - Flags
     - Use
   * - Starting point
     - ``-O2 -ffp-contract=fast``
     - Optimize kernels and permit fused multiply-add instructions.
   * - Local x86-64 CPU
     - ``-O2 -march=native -ffp-contract=fast``
     - Also enable instructions supported by the compilation host.
   * - Additional optimization
     - ``-O3 -ffp-contract=fast``
     - Compare against O2; O3 is not faster for every kernel.
   * - Relaxed floating-point arithmetic
     - ``-O3 -ffast-math -ffp-contract=fast``
     - Evaluate when the application's numerical requirements allow it.

On x86-64, ``-march=native`` can also be added to either O3 configuration.
It can produce binaries that do not run on other processors. Omit it or
select an explicit CPU target supported by every machine that will run the
kernel, including machines sharing a JIT cache. Other architectures have
their own CPU-target options.

``-ffp-contract=fast`` permits the compiler to combine multiplication and
addition into a fused multiply-add (FMA) instruction when the target CPU
supports it. FMA rounds once, rather than rounding the product and sum
separately, so results can differ in their last bits. This option does not
require ``-ffast-math``. Use ``-ffp-contract=off`` when separate rounding is
required.

``-ffast-math`` permits broader changes, including reassociation of
floating-point expressions and assumptions about NaNs, infinities and signed
zeros. Check accuracy and convergence for the application before adopting
it. Compare O2 and O3 using the same CPU-target and floating-point options
to distinguish optimization-level effects from changes in arithmetic.

These switches are compiler-specific. For Microsoft Visual C++, start with
``["/O2", "/fp:precise"]``. Visual Studio 2022 also provides ``/fp:contract``
to permit FMA, subject to the selected target's instruction support. See the
`GCC optimization options <https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html>`_,
`Clang floating-point options <https://clang.llvm.org/docs/UsersManual.html#controlling-floating-point-behavior>`_
and `MSVC floating-point options <https://learn.microsoft.com/en-us/cpp/build/reference/fp-specify-floating-point-behavior>`_
for details.

Setting flags with the FFCx JIT
-------------------------------

Pass a list of C compiler arguments as ``cffi_extra_compile_args`` to
``ffcx.codegeneration.jit.compile_forms``. Given a UFL form ``a``:

.. code-block:: python

   from pathlib import Path

   from ffcx.codegeneration.jit import compile_forms

   compiler_flags = ["-O2", "-ffp-contract=fast"]
   compiled_forms, module, code = compile_forms(
       [a],
       cache_dir=Path(".cache/ffcx"),
       cffi_extra_compile_args=compiler_flags,
   )

``compile_expressions`` accepts the same keyword argument. FFCx includes
the compiler arguments in its JIT signature, so changing them produces a
different cached module. Use ``cffi_verbose=True`` to see the compiler
invocation. Select the C compiler through the build environment, for
example by setting ``CC=clang`` before starting Python.

C compiler arguments are separate from FFCx code-generation options such
as ``scalar_type`` and ``sum_factorization``, which belong in the
``options`` dictionary. :doc:`sum-factorisation` is enabled by default
where supported.

Setting flags through DOLFINx
-----------------------------

When using DOLFINx, pass the compiler arguments through ``jit_options``:

.. code-block:: python

   from dolfinx import fem

   compiled_a = fem.form(
       a,
       jit_options={
           "cffi_extra_compile_args": ["-O2", "-ffp-contract=fast"],
       },
   )

Use ``form_compiler_options`` for FFCx code-generation options. DOLFINx
also supports a ``dolfinx_jit_options.json`` file in the working directory:

.. code-block:: json

   {
       "cffi_extra_compile_args": ["-O2", "-ffp-contract=fast"]
   }

An explicitly supplied ``jit_options`` value takes precedence over the
corresponding file setting.

Compiling generated C files
---------------------------

The ``ffcx`` command generates C source and header files. Apply compiler
flags when building that source, for example in a CMake target:

.. code-block:: cmake

   target_compile_options(my_kernels PRIVATE -O2 -ffp-contract=fast)

Use the compiler options appropriate to the selected compiler and target.
Changing FFCx code-generation options does not set these C compiler flags.
