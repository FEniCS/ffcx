Plan for test driven development of form compiler
=================================================

* Move expression code generation tests from site-packages/test_uflacs

* Repeat for interval, triangle, tetrahedron, quadrilateral, hexahedron:

    + Basic expressions of geometry point integral (x[] input)

        - Compute expression of x[] only

        - Compute ufl_cell.J -> J[]

        - Compute ufl_cell.detJ -> detJ

        - Compute ufl_cell.Jinv -> K[]

        - Compute ufl_cell.xi -> xi[] from x[], v0, K

    + Integration

        - Assume a quadrature rule (one midpoint just for testing)

        - Generate quadrature loop, defining xi[]

        - Compute x[] from cell_xi[]

        - Compute cell_xi[] from facet_xi[]

        - Compute x[] from facet_xi[]

        - Compute cell integration scaling factor (abs(detJ)?)

        - Compute facet integration scaling factor

    + Coefficients

        - Get scalar piecewise constant coefficients

        - Get vector piecewise constant coefficients

        - Get tensor piecewise constant coefficients (with symmetries)

        - Compute non-constant coefficients from dofs and basis functions

        - Compute local coefficient gradients

        - Compute mapped global coefficient gradients

    + Secondary geometry computations

        - Compute cell barycenter

        - Compute cell volume

        - Compute cell radius

        - Compute cell normal n[]?

        - Compute facet barycenter

        - Compute facet area

        - Compute facet radius

        - Compute facet normal n[]

    + Form structure

        - Handle loop over test space

        - Handle loops over test and trial space

        - Compute piecewise constant basis functions

        - Compute non-constant basis functions

        - Compute basis function gradients (with mapping K applied)

* Integrate with ffc:

    + Add clean slate uflacs representation

    + Make uflacs representation delegate to quadrature representation for elements etc.

    + Make uflacs representation call uflacs to generate tabulate_tensor body

        - point_integral

        - cell_integral

        - exterior_facet_integral

        - interior_facet_integral

* Saving optimization and tuning of algorithms for later when stuff works

