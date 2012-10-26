Plan for test driven development of form compiler
=================================================

* Repeat for interval, triangle, tetrahedron, quadrilateral, hexahedron:

    + Expression code generation tests (move tests from site-packages)

        - Compute expressions of x only, assuming x[] given

    + Initial geometry computations

        - Compute J[], detJ, Jinv[]

        - Compute x[] from xi[]

        - Compute xi[] from x[]

        - Compute integration scaling factor

    + Initial coefficients

        - Handle piecewise constant coefficients

    + Form structure

        - Compute functional of x assuming quadrature rule given

        - Handle piecewise constant basis functions

        - Handle non-constant basis functions

        - Handle basis function gradients (with mapping applied)

    + Secondary coefficients

        - Handle non-constant coefficients

        - Handle coefficient gradients (with mapping applied)

    + Secondary geometry computations

        - Compute facet integration scaling factor

        - Compute n[]

        - Compute volume, etc.

* Integrate with sfc/ffc!

* Saving optimization and tuning of algorithms for later when stuff works

