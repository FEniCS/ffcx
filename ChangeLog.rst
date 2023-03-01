Changelog
=========

0.6.0
-----
See https://github.com/FEniCS/ffcx/compare/v0.5.0...v0.6.0

0.5.0
-----
See: https://github.com/FEniCS/ffcx/compare/v0.5.0...v0.4.0

0.4.0
-----
See: https://github.com/FEniCS/ffcx/compare/v0.4.0...v0.3.0

0.3.0
-----
See: https://github.com/FEniCS/ffcx/compare/v0.3.0...v0.2.0

0.2.0
-----

- No changes

0.1.0
-----
Alpha release of ffcx

2018.2.0.dev0
-------------

- No changes

2018.1.0.dev0 (no release)
--------------------------

- Forked FFCx

2017.2.0 (2017-12-05)
---------------------

- Some fixes for ufc::eval for esoteric element combinations
- Reimplement code generation for all ufc classes with new class
  ufc::coordinate_mapping which can map between coordinates, compute
  jacobians, etc. for a coordinate mapping parameterized by a specific
  finite element.
- New functions in ufc::finite_element:
  - evaluate_reference_basis
  - evaluate_reference_basis_derivatives
  - transform_reference_basis_derivatives
  - tabulate_reference_dof_coordinates
- New functions in ufc::dofmap:
  - num_global_support_dofs
  - num_element_support_dofs
- Improved docstrings for parts of ufc.h
- FFC now accepts Q and DQ finite element families defined on quadrilaterals and hexahedrons
- Some fixes for ufc_geometry.h for quadrilateral and hexahedron cells

2017.1.0.post2 (2017-09-12)
---------------------------

- Change PyPI package name to fenics-ffc.

2017.1.0 (2017-05-09)
---------------------

- Let ffc -O parameter take an optional integer level like -O2, -O0
- Implement blockwise optimizations in uflacs code generation
- Expose uflacs optimization parameters through parameter system

2016.2.0 (2016-11-30)
---------------------

- Jit compiler now compiles elements separately from forms to avoid duplicate work
- Add parameter max_signature_length to optionally shorten signatures in the jit cache
- Move uflacs module into ffc.uflacs
- Remove installation of pkg-config and CMake files (UFC path and
  compiler flags are available from ffc module)
- Add dependency on dijitso and remove dependency on instant
- Add experimental Bitbucket pipelines
- Tidy the repo after UFC and UFLACS merge, and general spring cleanup. This
  includes removal of instructions how to merge two repos, commit hash
  c8389032268041fe94682790cb773663bdf27286.

2016.1.0 (2016-06-23)
---------------------

- Add function get_ufc_include to get path to ufc.h
- Merge UFLACS into FFC
- Generalize ufc interface to non-affine parameterized coordinates
- Add ufc::coordinate_mapping class
- Make ufc interface depend on C++11 features requiring gcc version >= 4.8
- Add function ufc_signature() to the form compiler interface
- Add function git_commit_hash()

1.6.0 (2015-07-28)
------------------

- Rename and modify a number of UFC interface functions. See docstrings in ufc.h for details.
- Bump required SWIG version to 3.0.3
- Disable dual basis (tabulate_coordinates and evaluate_dofs) for enriched
  elements until correct implementation is brought up

1.5.0 (2015-01-12)
------------------

- Remove FErari support
- Add support for new integral type custom_integral
- Support for new form compiler backend "uflacs", downloaded separately

1.4.0 (2014-06-02)
------------------

- Add support for integrals that know which coefficients they use
- Many bug fixes for facet integrals over manifolds
- Merge UFC into FFC; ChangeLog for UFC appended below
- Various updates mirroring UFL changes
- Experimental: New custom integral with user defined quadrature points

1.3.0 (2014-01-07)
------------------

- Fix bug with runtime check of SWIG version
- Move DOLFIN wrappers here from DOLFIN
- Add support for new UFL operators cell_avg and facet_avg
- Add new reference data handling system, now data is kept in an external repository
- Fix bugs with ignoring quadrature rule arguments
- Use cpp optimization by default in jit compiler

1.2.0 (2013-03-24)
------------------

- New feature: Add basic support for point integrals on vertices
- New feature: Add general support for m-dimensional cells in n-dimensional space (n >= m, n, m = 1, 2, 3)

1.1.0 (2013-01-07)
------------------

- Fix bug for Conditionals related to DG constant Coefficients. Bug #1082048.
- Fix bug for Conditionals, precedence rules for And and Or. Bug #1075149.
- Changed data structure from list to deque when pop(0) operation is needed, speeding up split_expression operation considerable
- Other minor fixes

1.0.0 (2011-12-07)
------------------

- Issue warning when form integration requires more than 100 points

1.0-rc1 (2011-11-28)
--------------------

- Fix bug with coordinates on facet integrals (intervals). Bug #888682.
- Add support for FacetArea, new geometric quantity in UFL.
- Fix bug in optimised quadrature code, AlgebraOperators demo. Bug #890859.
- Fix bug with undeclared variables in optimised quadrature code. Bug #883202.

1.0-beta2 (2011-10-11)
----------------------

- Added support for bessel functions, bessel_* (I,J,K,Y), in UFL.
- Added support for error function, erf(), new math function in UFL.
- Fix dof map 'need_entities' for Real spaces
- Improve performance for basis function computation

1.0-beta (2011-08-11)
---------------------

- Improve formatting of floats with up to one non-zero decimal place.
- Fix bug involving zeros in products and sums. Bug #804160.
- Fix bug for new conditions '&&', '||' and '!' in UFL. Bug #802560.
- Fix bug involving VectorElement with dim=1. Bug #798578.
- Fix bug with mixed element of symmetric tensor elements. Bug #745646.
- Fix bug when using geometric coordinates with one quadrature point

0.9.10 (2011-05-16)
-------------------

- Change license from GPL v3 or later to LGPL v3 or later
- Add some schemes for low-order simplices
- Request quadrature schemes by polynomial degree (not longer by number
  of points in each direction)
- Get quadrature schemes via ffc.quadrature_schemes
- Improved lock handling in JIT compiler
- Include common_cell in form signature
- Add possibility to set swig binary and swig path

0.9.9 (2011-02-23)
------------------

- Add support for generating error control forms with option -e
- Updates for UFC 2.0
- Set minimal degree to 1 in automatic degree selection for expressions
- Add command-line option -f no_ferari
- Add support for plotting of elements
- Add utility function compute_tensor_representation

0.9.4 (2010-09-01)
------------------

- Added memory cache in jit(), for preprocessed forms
- Added support for Conditional and added demo/Conditional.ufl.
- Added support for new geometric quantity Circumradius in UFL.
- Added support for new geometric quantity CellVolume in UFL.

0.9.3 (2010-07-01)
------------------

- Make global_dimension for Real return an int instead of double, bug # 592088
- Add support for facet normal in 1D.
- Expose -feliminate_zeros for quadrature optimisations to give user more
  control
- Remove return of form in compile_form
- Remove object_names argument to compile_element
- Rename ElementUnion -> EnrichedElement
- Add support for tan() and inverse trigonometric functions
- Added support for ElementUnion (i.e. span of combinations of elements)
- Added support for Bubble elements
- Added support for UFL.SpatialCoordinate.

0.9.2 (2010-02-17)
------------------

- Bug fix in removal of unused variables in Piola-mapped terms for tensor
  representation

0.9.1 (2010-02-15)
------------------

- Add back support for FErari optimizations
- Bug fixes in JIT compiler

0.9.0 (2010-02-02)
------------------

- Updates for FIAT 0.9.0
- Updates for UFC 1.4.0 (now supporting the full interface)
- Automatic selection of representation
- Change quadrature_order --> quadrature_degree
- Split compile() --> compile_form(), compile_element()
- Major cleanup and reorganization of code (flatter directories)
- Updates for changes in UFL: Argument, Coefficient, FormData

0.7.1
-----

- Handle setting quadrature degree when it is set to None in UFL form
- Added demo: HyperElasticity.ufl

0.7.0
-----

- Move contents of TODO to: https://blueprints.launchpad.net/ffc
- Support for restriction of finite elements to only consider facet dofs
- Use quadrature_order from metadata when integrating terms using tensor representation
- Use loop to reset the entries of the local element tensor
- Added new symbolic classes for quadrature optimisation (speed up compilation)
- Added demos: Biharmonic.ufl, div(grad(v)) term;
               ReactionDiffusion.ufl, tuple notation;
               MetaData.ufl, how to attach metadata to the measure;
               ElementRestriction.ufl, restriction of elements to facets
- Tabulate the coordinates of the integration points in the tabulate_tensor() function
- Change command line option '-f split_implementation' -> '-f split'
- Renaming of files and restructuring of the compiler directory
- Added option -q rule (--quadrature-rule rule) to specify which rule to use
  for integration of a given integral. (Can also bet set through the metadata
  through "quadrature_rule"). No rules have yet been implemented, so default
  is the FIAT rule.
- Remove support for old style .form files/format

0.6.2 (2009-04-07)
------------------

- Experimental support for UFL, supporting both .form and .ufl
- Moved configuration and construction of python extension module to ufc_module

0.6.1 (2009-02-18)
------------------

- Initial work on UFL transition
- Minor bug fixes
- The version of ufc and swig is included in the form signature
- Better system configuration for JIT compiled forms
- The JIT compiled python extension module use shared_ptr for all classes

0.6.0 (2009-01-05)
------------------

- Update DOLFIN output format (-l dolfin) for DOLFIN 0.9.0
- Cross-platform fixes for test scripts
- Minor bug fix for quadrature code generation (forms affected by this bug would not be able to compile
- Fix bug with output of ``*.py``.
- Permit dot product bewteen rectangular matrices (Frobenius norm)

0.5.1 (2008-10-20)
------------------

- New operator skew()
- Allow JIT compilation of elements and dof maps
- Rewrite JIT compiler to rely on Instant for caching
- Display flop count for evaluating the element tensor during compilation
- Add arguments language and representation to options dictionary
- Fix installation on Windows
- Add option -f split_implementation for separate .h and .cpp files

0.5.0 (2008-06-23)
------------------

- Remove default restriction +/- for Constant
- Make JIT optimization (-O0 / -O2) optional
- Add in-memory cache to speed up JIT compiler for repeated assembly
- Allow subdomain integrals without needing full range of integrals
- Allow simple subdomain integral specification dx(0), dx(1), ds(0) etc

0.4.5 (2008-04-30)
------------------

- Optimizations in generated quadrature code
- Change formatting of floats from %g to %e, fixes problem with too long integers
- Bug fix for order of values in interpolate_vertex_values, now according to UFC
- Speed up JIT compiler
- Add index ranges to form printing
- Throw runtime error in functions not generated
- Update DOLFIN format for new location of include files

0.4.4 (2008-02-18)
------------------

- RT, BDM, BDFM and Nedelec now working in 2D and 3D
- New element type QuadratureElement
- Add support for 1D elements
- Add experimental support for new Darcy-Stokes element
- Use FIAT transformed spaces instead of mapping in FFC
- Updates for UFC 1.1
- Implement caching of forms/modules in ~/.ffc/cache for JIT compiler
- Add script ffc-clean
- New operators lhs() and rhs()
- Bug fixes in simplify
- Bug fixes for Nedelec and BDFM
- Fix bug in mult()
- Fix bug with restrictions on exterior facet integrals
- Fix bug in grad() for vectors
- Add divergence operator for matrices

0.4.3 (2007-10-23)
------------------

- Require FIAT to use UFC reference cells
- Fix bug in form simplification
- Rename abs --> modulus to avoid conflict with builtin abs
- Fix bug in operators invert, abs, sqrt
- Fix bug in integral tabulation
- Add BDFM and Nedelec elements (nonworking)
- Fix bug in JIT compiler

0.4.2 (2007-08-31)
------------------

- Change license from GPL v2 to GPL v3 or later
- Add JIT (just-in-time) compiler
- Fix bug for constants on interior facets

0.4.1 (2007-06-22)
------------------

- Fix bug in simplification of forms
- Optimize removal of unused terms in code formattting

0.4.0 (2007-06-20)
------------------

- Move to UFC interface for code generation
- Major rewrite, restructure, cleanup
- Add support for Brezzi-Douglas-Marini (BDM) elements
- Add support for Raviart-Thomas (RT) elements
- Add support for Discontinuous Galerkin (DG) methods
- Operators jump() and avg()
- Add quadrature compilation mode (experimental)
- Simplification of forms
- Operators sqrt(), abs() and inverse
- Improved Python interface
- Add flag -f precision=n
- Generate code for basis functions and derivatives
- Use Set from set module for Python2.3 compatibility

0.3.5 (2006-12-01)
------------------

- Bug fixes
- Move from Numeric to numpy

0.3.4 (2006-10-27)
------------------

- Updates for new DOLFIN mesh library
- Add support for evaluation of functionals
- Add operator outer() for outer product of vector-valued functions
- Enable optimization of linear forms (in addition to bilinear forms)
- Remove DOLFIN SWIG format
- Fix bug in ffc -v/--version (thanks to Ola Skavhaug)
- Consolidate DOLFIN and DOLFIN SWIG formats (patch from Johan Jansson)
- Fix bug in optimized compilation (-O) for some forms ("too many values to unpack")

0.3.3 (2006-09-05)
------------------

- Fix bug in operator div()
- Add operation count (number of multiplications) with -d0
- Add hint for printing more informative error messages (flag -d1)
- Modify implementation of vertexeval()
- Add support for boundary integrals (Garth N. Wells)

0.3.2 (2006-04-01)
------------------

- Add support for FErari optimizations, new flag -O

0.3.1 (2006-03-28)
------------------

- Remove verbose output: silence means success
- Generate empty boundary integral eval() to please Intel C++ compiler
- New classes TestFunction and TrialFunction

0.3.0 (2006-03-01)
------------------

- Work on manual, document command-line and user-interfaces
- Name change: u --> U
- Add compilation of elements without form
- Add generation of FiniteElementSpec in DOLFIN formats
- Fix bugs in raw and XML formats
- Fix bug in LaTeX format
- Fix path and predefine tokens to enable import in .form file
- Report number of entries in reference tensor during compilation

0.2.5 (2005-12-28)
------------------

- Add demo Stabilization.form
- Further speedup computation of reference tensor (use ufunc Numeric.add)

0.2.4 (2005-12-05)
------------------

- Report time taken to compute reference tensor
- Restructure computation of reference tensor to use less memory.
  As a side effect, the speed has also been improved.
- Update for DOLFIN name change node --> vertex
- Update finite element interface for DOLFIN
- Check for FIAT bug in discontinuous vector Lagrange elements
- Fix signatures for vector-valued elements

0.2.3 (2005-11-28)
------------------

- New fast Numeric/BLAS based algorithm for computing reference tensor
- Bug fix: reassign indices for complete subexpressions
- Bug fix: operator Function * Integral
- Check tensor notation for completeness
- Bug fix: mixed elements with more than two function spaces
- Don't declare unused coefficients (or gcc will complain)

0.2.2 (2005-11-14)
------------------

- Add command-line argument -v / --version
- Add new operator mean() for projection onto piecewise constants
- Add support for projections
- Bug fix for higher order mixed elements: declaration of edge/face_ordering
- Generate code for sub elements of mixed elements
- Add new test form: TensorWeighteLaplacian
- Add new test form: EnergyNorm
- Fix bugs in mult() and vec() (skavhaug)
- Reset correct entries of G for interior in BLAS mode
- Only assign to entries of G that meet nonzero entries of A in BLAS mode

0.2.1 (2005-10-11)
------------------

- Only generate declarations that are needed according to format
- Check for missing options and add missing default options
- Simplify usage of FFC as Python module: from ffc import *
- Fix bug in division with constants
- Generate output for BLAS (with option -f blas)
- Add new XML output format
- Remove command-line option --license (collect in compiler options -f)
- Modify demo Mass.form to use 3:rd order Lagrange on tets
- Fix bug in dofmap() for equal order mixed elements
- Add compiler option -d debuglevel
- Fix Python Numeric bug: vdot --> dot

0.2.0 (2005-09-23)
------------------

- Generate function vertexeval() for evaluation at vertices
- Add support for arbitrary mixed elements
- Add man page
- Work on manual, chapters on form language, quickstart and installation
- Handle exceptions gracefully in command-line interface
- Use new template fenicsmanual.cls for manual
- Add new operators grad, div, rot (curl), D, rank, trace, dot, cross
- Factorize common reference tensors from terms with equal signatures
- Collect small building blocks for form algebra in common module tokens.py

0.1.9 (2005-07-05)
------------------

- Complete support for general order Lagrange elements on triangles and tetrahedra
- Compute reordering of dofs on tets correctly
- Update manual with ordering of dofs
- Break compilation into two phases: build() and write()
- Add new output format ASE (Matt Knepley)
- Improve python interface to FFC
- Remove excessive logging at compilation
- Fix bug in raw output format

0.1.8 (2005-05-17)
------------------

- Access data through map in DOLFIN format
- Experimental support for computation of coordinate maps
- Add first draft of manual
- Experimental support for computation of dof maps
- Allow specification of the number of components for vector Lagrange
- Count the number of zeros dropped
- Fix bug in handling command-line arguments
- Use module sets instead of built-in set (fix for Python 2.3)
- Handle constant indices correctly (bug reported by Garth N. Wells)

0.1.7 (2005-05-02)
------------------

- Write version number to output
- Add command-line option for choosing license
- Display usage if no input is given
- Bug fix for finding correct prefix of file name
- Automatically choose name of output file (if not supplied)
- Use FIAT tabulation mode for vector-valued elements (speedup a factor 5)
- Use FIAT tabulation mode for scalar elements (speedup a factor 1000)
- Fig bug in demo elasticity.form (change order of u and v)
- Make references to constants const in DOLFIN format
- Don't generate code for unused entries of geometry tensor
- Update formats to write numeric constants with full precision

0.1.6 (2005-03-17)
------------------

- Add support for mixing multiple different finite elements
- Add support for division with constants
- Fix index bug (reverse order of multi-indices)

0.1.5 (2005-03-14)
------------------

- Automatically choose the correct quadrature rule for precomputation
- Add test program for verification of FIAT quadrature rules
- Fix bug for derivative of sum
- Improve common interface for debugging: add indentation
- Add support for constants
- Fix bug for sums of more than one term (make copies of references in lists)
- Add '_' in naming of geometry tensor (needed for large dimensions)
- Add example elasticity.form
- Cleanup build_indices()

0.1.4-1 (2005-02-07)
--------------------

- Fix version number and remove build directory from tarball

0.1.4 (2005-02-04)
------------------

- Fix bug for systems, seems to work now
- Add common interface for debugging
- Modify DOLFIN output to initialize functions
- Create unique numbers for each function
- Use namespaces for DOLFIN output instead of class names
- Temporary implementation of dof mapping for vector-valued elements
- Make DOLFIN output format put entries into PETSc block
- Change name of coefficient data: c%d[%d] -> c[%d][%d]
- Change ordering of basis functions (one component at a time)
- Add example poissonsystem.form
- Modifications for new version of FIAT (FIAT-L)
  FIAT version 0.1 a factor 5 slower (no memoization)
  FIAT version 0.1.1 a little faster, only a factor 2 slower
- Add setup.py script

0.1.3 (2004-12-06)
------------------

- Fix bug in DOLFIN format (missing value when zero)
- Add output of reference tensor to LaTeX format
- Make raw output format print data with full precision
- Add component diagram
- Change order of declaration of basis functions
- Add new output format raw

0.1.2 (2004-11-17)
------------------

- Add command-line interface ffc
- Add support for functions (coefficients)
- Add support for constants
- Allow multiple forms (left- and right-hand side) in same file
- Add test examples: poisson.form, mass.form, navierstokes.form
- Wrap FIAT to create vector-valued finite element spaces
- Check ranks of operands
- Clean up algebra, add base class Element
- Add some documentation (class diagram)
- Add support for LaTeX output

0.1.1-1 (2004-11-10)
--------------------

- Add missing file declaration.py

0.1.1 (2004-11-10)
------------------

- Make output variable names configurable
- Clean up DOLFIN code generation
- Post-process form to create reference, geometry, and element tensors
- Experimental support for general tensor-valued elements
- Clean up and improve index reassignment
- Use string formatting for generation of output
- Change index ordering to access row-wise

0.1.0 (2004-10-22)
------------------

- First iteration of the FEniCS Form Compiler
- Change boost::shared_ptr --> std::shared_ptr

ChangeLog for UFC
=================

UFC was merged into FFC 2014-02-18. Below is the ChangeLog for
UFC at the time of the merge. From this point onward, UFC version
numbering restarts at the same version number as FFC and the rest
of FEniCS.

2.3.0 (2014-01-07)
------------------

- Use std::vector<std::vector<std::size_t> > for topology data
- Remove vertex coordinates from ufc::cell
- Improve detection of compatible Python libraries
- Add current swigversion to the JIT compiled extension module
- Remove dofmap::max_local_dimension()
- Remove cell argument from dofmap::local_dimension()

2.2.0 (2013-03-24)
------------------

- Add new class ufc::point_integral
- Use CMake to configure JIT compilation of forms
- Generate UseUFC.cmake during configuration
- Remove init_mesh(), init_cell(), init_mesh_finalize()
- Remove ufc::mesh and add a vector of num_mesh_entities to global_dimension() and tabulate_dofs().

2.1.0 (2013-01-07)
------------------

- Fix bug introduced by SWIG 2.0.5, which treated uint as Python long
- Add optimization SWIG flags, fixing bug lp:987657

2.0.5 (2011-12-07)
------------------

- Improve configuration of libboost-math

2.0.4 (2011-11-28)
------------------

- Add boost_math_tr1 to library flags when JIT compiling an
  extension module

2.0.3 (2011-10-26)
------------------

- CMake config improvements

2.0.2 (2011-08-11)
------------------

- Some tweaks of installation

2.0.1 (2011-05-16)
------------------

- Make SWIG version >= 2.0 a requirement
- Add possibility to set swig binary and swig path
- Add missing const for map_{from,to}_reference_cell

2.0.0 (2011-02-23)
------------------

- Add quadrature version of tabulate_tensor
- Add finite_element::map_{from,to}_reference_cell
- Add finite_element::{topological,geometric}_dimension
- Add dofmap::topological_dimension
- Rename num_foo_integrals --> num_foo_domains
- Rename dof_map --> dofmap
- Add finite_element::create
- Add dofmap::create

1.4.2 (2010-09-01)
------------------

- Move to CMake build system

1.4.1 (2010-07-01)
------------------

- Make functions introduced in UFC 1.1 mandatory (now pure virtual)
- Update templates to allow constructor arguments in form classes

1.4.0 (2010-02-01)
------------------

- Changed behavior of create_foo_integral (returning 0 when integral is 0)
- Bug fixes in installation

1.2.0 (2009-09-23)
------------------

- Add new function ufc::dof_map::max_local_dimension()
- Change ufc::dof_map::local_dimension() to ufc::dof_map::local_dimension(const ufc::cell c)

1.1.2 (2009-04-07)
------------------

- Added configuration and building of python extension module to ufc_utils.build_ufc_module

1.1.1 (2009-02-20)
------------------

- The extension module is now not built, if the conditions for shared_ptr are not met
- Added SCons build system
- The swig generated extension module will be compiled with shared_ptr support if boost is found on system and swig is of version 1.3.35 or higher
- The swig generated extension module is named ufc.py and expose all ufc base classes to python
- Added a swig generated extention module to ufc. UFC now depends on swig
- Changed name of the python utility module from "ufc" to "ufc_utils"

1.1.0 (2008-02-18)
------------------

- Add new function ufc::finite_element::evaluate_dofs
- Add new function ufc::finite_element::evaluate_basis_all
- Add new function ufc::finite_element::evaluate_basis_derivatives_all
- Add new function ufc::dof_map::geometric_dimension
- Add new function ufc::dof_map::num_entity_dofs
- Add new function ufc::dof_map::tabulate_entity_dofs

1.0.0 (2007-06-17)
------------------

- Release of UFC 1.0
