# -*- coding: utf-8 -*-
# Copyright (C) 2009-2016 Anders Logg and Martin Sandve Alnæs
#
# This file is part of UFLACS.
#
# UFLACS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFLACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>.

# Note: Most of the code in this file is a direct translation from the old implementation in FFC


from ufl import product
from ffc.uflacs.backends.ufc.generator import ufc_generator
from ffc.uflacs.backends.ufc.utils import generate_return_new_switch, generate_return_int_switch, generate_error


class ufc_finite_element(ufc_generator):
    "Each function maps to a keyword in the template. See documentation of ufc_generator."
    def __init__(self):
        ufc_generator.__init__(self, "finite_element")

    def cell_shape(self, L, cell_shape):
        return L.Return(L.Symbol("ufc::shape::" + cell_shape))

    def topological_dimension(self, L, topological_dimension):
        return L.Return(topological_dimension)

    def geometric_dimension(self, L, geometric_dimension):
        return L.Return(geometric_dimension)

    def space_dimension(self, L, space_dimension):
        return L.Return(space_dimension)

    def value_rank(self, L, value_shape):
        return L.Return(len(value_shape))

    def value_dimension(self, L, value_shape):
        return generate_return_int_switch(L, "i", value_shape, 1)

    def value_size(self, L, value_shape):
        return L.Return(product(value_shape))

    def reference_value_rank(self, L, reference_value_shape):
        return L.Return(len(reference_value_shape))

    def reference_value_dimension(self, L, reference_value_shape):
        return generate_return_int_switch(L, "i", reference_value_shape, 1)

    def reference_value_size(self, L, reference_value_shape):
        return L.Return(product(reference_value_shape))

    def degree(self, L, degree):
        return L.Return(degree)

    def family(self, L, family):
        return L.Return(L.LiteralString(family))

    def num_sub_elements(self, L, num_sub_elements):
        return L.Return(num_sub_elements)

    def create_sub_element(self, L, ir):
        classnames = ir["create_sub_element"]
        return generate_return_new_switch(L, "i", classnames, factory=ir["jit"])

    def evaluate_basis(self, L, ir, parameters): # FIXME: Get rid of this
        from ffc.codegeneration import _evaluate_basis, indent
        return indent(_evaluate_basis(ir["evaluate_basis"]), 4)

    def evaluate_basis_all(self, L, ir, parameters):
        # FIXME: port this, then translate into reference version
        from ffc.codegeneration import _evaluate_basis_all, indent
        return indent(_evaluate_basis_all(ir["evaluate_basis"]), 4)

    def evaluate_basis_derivatives(self, L, ir, parameters): # FIXME: Get rid of this
        # FIXME: port this, then translate into reference version
        from ffc.codegeneration import _evaluate_basis_derivatives, indent
        return indent(_evaluate_basis_derivatives(ir["evaluate_basis"]), 4)

    def evaluate_basis_derivatives_all(self, L, ir, parameters):
        # FIXME: port this, then translate into reference version
        from ffc.codegeneration import _evaluate_basis_derivatives_all, indent
        return indent(_evaluate_basis_derivatives_all(ir["evaluate_basis"]), 4)

    def evaluate_dof(self, L, ir, parameters): # FIXME: Get rid of this
        # FIXME: port this, then translate into reference version
        # Codes generated together
        from ffc.evaluatedof import evaluate_dof_and_dofs
        from ffc.codegeneration import indent
        (evaluate_dof_code, evaluate_dofs_code) \
          = evaluate_dof_and_dofs(ir["evaluate_dof"])
        return indent(evaluate_dof_code, 4)

    def evaluate_dofs(self, L, ir, parameters):
        """Generate code for evaluate_dofs."""
        """
        - evaluate_dof needs to be split into invert_mapping + evaluate_dof or similar?

          f = M fhat;  nu(f) = nu(M fhat) = nuhat(M^-1 f) = sum_i w_i M^-1 f(x_i)

          // Get fixed set of points on reference element
          num_points = element->num_dof_evaluation_points();
          double X[num_points*tdim];
          element->tabulate_dof_evaluation_points(X);

          // Compute geometry in these points
          domain->compute_geometry(reference_points, num_point, X, J, detJ, K, coordinate_dofs, cell_orientation);

          // Computed by dolfin
          for ip
            fvalues[ip][:] = f.evaluate(point[ip])[:];

          // Finally: nu_j(f) = sum_component sum_ip weights[j][ip][component] fvalues[ip][component]
          element->evaluate_dofs(fdofs, fvalues, J, detJ, K)
        """
        # FIXME: port this, then translate into reference version
        # Codes generated together
        from ffc.evaluatedof import evaluate_dof_and_dofs
        from ffc.codegeneration import indent
        (evaluate_dof_code, evaluate_dofs_code) \
          = evaluate_dof_and_dofs(ir["evaluate_dof"])
        return indent(evaluate_dofs_code, 4)

    def interpolate_vertex_values(self, L, ir, parameters): # FIXME: port this
        # FIXME: port this, then translate into reference version
        from ffc.codegeneration import interpolate_vertex_values, indent
        return indent(interpolate_vertex_values(ir["interpolate_vertex_values"]), 4)

    def tabulate_dof_coordinates(self, L, ir, parameters): # FIXME: port this
        from ffc.codegeneration import _tabulate_dof_coordinates, indent
        return indent(_tabulate_dof_coordinates(ir["tabulate_dof_coordinates"]), 4)

        # TODO: For a transition period, let finite_element and dofmap depend on a class affine_<cellname>_domain?
        # TODO: Call _tabulate_dof_reference_coordinates to tabulate X[ndofs][tdim],
        # then call affine_domain::compute_physical_coordinates(x, X, coordinate_dofs)

        ir = ir["tabulate_dof_coordinates"]

        # Raise error if tabulate_dof_coordinates is ill-defined
        if not ir:
            msg = "tabulate_dof_coordinates is not defined for this element"
            return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])

        # Extract coordinates and cell dimension
        gdim = ir["gdim"]
        tdim = ir["tdim"]
        points = ir["points"]

        # Output argument
        dof_coordinates = L.FlattenedArray(L.Symbol("dof_coordinates"),
                                           dims=(len(points), gdim))

        # Input argument
        coordinate_dofs = L.Symbol("coordinate_dofs")

        # Loop indices
        i = L.Symbol("i")
        k = L.Symbol("k")
        ip = L.Symbol("ip")

        # Basis symbol
        phi = L.Symbol("phi")

        # TODO: This is used where we still assume an affine mesh. Get rid of all places that use it.
        from ffc.evaluatedof import affine_weights
        # TODO: This code assumes an affine coordinate field.
        #       Ok for now in here, this function must be removed anyway.
        # Create code for evaluating affine coordinate basis functions
        num_scalar_xdofs = tdim + 1
        cg1_basis = affine_weights(tdim)
        phi_values = [phi_comp for X in points for phi_comp in cg1_basis(X)]
        assert len(phi_values) == len(points) * num_scalar_xdofs

        code = [
            L.ArrayDecl("const double", phi,
                        (len(points) * num_scalar_xdofs,),
                        values=phi_values),
            L.ForRange(ip, 0, len(points), body=
                L.ForRange(i, 0, gdim, body=
                    L.ForRange(k, 0, num_scalar_xdofs, body=
                        L.AssignAdd(dof_coordinates[ip][i],
                                    coordinate_dofs[gdim*k + i]
                                    * phi[ip*num_scalar_xdofs + k])))),
            ]
        return code

    def tabulate_reference_dof_coordinates(self, L, ir, parameters):
        # TODO: Change signature to avoid copy? E.g.
        # virtual const std::vector<double> & tabulate_reference_dof_coordinates() const = 0;
        # See integral::enabled_coefficients for example

        # TODO: ensure points is a numpy array,
        #   get tdim from points.shape[1],
        #   place points in ir directly instead of the subdict
        ir = ir["tabulate_dof_coordinates"]

        # Raise error if tabulate_reference_dof_coordinates is ill-defined
        if not ir:
            msg = "tabulate_reference_dof_coordinates is not defined for this element"
            return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])

        # Extract coordinates and cell dimension
        tdim = ir["tdim"]
        points = ir["points"]

        # Output argument
        reference_dof_coordinates = L.Symbol("reference_dof_coordinates")

        # Reference coordinates
        dof_X = L.Symbol("dof_X")
        dof_X_values = [X[jj] for X in points for jj in range(tdim)]
        decl = L.ArrayDecl("static const double", dof_X,
                           (len(points) * tdim,), values=dof_X_values)
        copy = L.MemCopy(dof_X, reference_dof_coordinates, tdim*len(points))

        code = [decl, copy]
        return code

    def evaluate_reference_basis(self, L, ir, parameters):
        data = ir["evaluate_basis"]
        #from ffc.uflacs.backends.ufc.evaluatebasis import generate_evaluate_reference_basis
        #return generate_evaluate_reference_basis(L, data, parameters)
        return L.Comment("Missing implementation")

    def evaluate_reference_basis_derivatives(self, L, ir, parameters):
        data = ir["evaluate_basis"]
        #from ffc.uflacs.backends.ufc.evalderivs import generate_evaluate_reference_basis_derivatives
        #return generate_evaluate_reference_basis_derivatives(L, data, parameters)
        return L.Comment("Missing implementation")


"""
TODO: Document new ufc functions evaluate_reference_basis and evaluate_reference_basis_derivatives

TODO: Add support for mappings to finite_element, something like:

    /// Return true if the basis needs to be mapped between physical and reference frames
    virtual bool needs_mapping() const = 0;

    /// Map values from reference frame to physical frame
    virtual void map_from_reference_values(double * values, const double * J) const = 0;

    /// Map values from physical frame to reference frame
    virtual void map_to_reference_values(double * values, const double * J) const = 0;

    // TODO: Need mapping of derivatives as well
"""


"""
TODO: Remove unused ufc::cell from interpolate_vertex_values
"""


''' /// FUTURE SPLIT IMPLEMENTATION OF EVALUATE_BASIS:
    /// Evaluate basis function i at given point x in cell
    virtual void evaluate_basis(std::size_t i,
                                double* values,
                                const double* x,
                                const double* coordinate_dofs,
                                int cell_orientation) const;
    /// Evaluate all basis functions at given point x in cell
    virtual void evaluate_basis_all(double* values,
                                    const double* x,
                                    const double* coordinate_dofs,
                                    int cell_orientation) const
    ... and derivatives
    {
      const std::size_t gdim = 3;
      const std::size_t tdim = 2;
      const std::size_t num_points = 1;

      // domain::
      double X[num_points*tdim]; // X[i] -> X[ip*tdim + i]
      compute_reference_coordinates(X, num_points, x, coordinate_dofs, cell_orientation);

      // domain::
      double J[num_points*gdim*tdim]; // J[i,j] -> J[ip*gdim*tdim + i*tdim + j]
      compute_jacobians(J, num_points, X, coordinate_dofs, cell_orientation);

      // domain::
      double detJ[num_points]; // detJ -> detJ[ip]
      compute_jacobian_determinants(detJ, num_points, J);

      // domain::
      double K[num_points*tdim*gdim]; // K[i,j] -> K[ip*tdim*gdim + i*gdim + j]
      compute_jacobian_inverses(K, num_points, J, detJ);

      // domain:: (inverse of compute_reference_coordinates)
      //double x[num_points*gdim]; // x[i] -> x[ip*gdim + i]
      //compute_physical_coordinates(x, num_points, X, K, coordinate_dofs, cell_orientation);

      // domain:: (combining the above)
      //compute_geometry(x, J, detJ, K, num_points, X, coordinate_dofs, cell_orientation);

      // phi[ip*ndofs*rvs + idof*rvs + jcomp]
      double reference_basis_values[num_points*num_dofs*reference_value_size];
      compute_reference_basis(reference_basis_values, num_points, X);

      // phi[ip*nder*ndofs*rvs + iderivative*ndofs*rvs + idof*rvs + jcomp]
      double reference_basis_derivatives[num_points*num_derivatives*num_dofs*reference_value_size];
      compute_reference_basis_derivatives(reference_basis_derivatives, derivative_order, num_points, X);

      double physical_basis_values[num_points*num_dofs*value_size]; // phi -> phi[ip*ndofs*pvs + idof*pvs + icomp]
      compute_physical_basis[_derivatives](physical_basis_values, num_points, reference_basis_values, J, detJ, K);
    }
'''
