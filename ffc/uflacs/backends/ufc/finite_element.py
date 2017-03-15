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



"""
// TODO: UFC functions that return constant arrays
// could use something like this to reduce copying,
// returning pointers to static data:

  template<T>
  class array_view
  {
  public:
    array_view(std::size_t size, const T * data):
      size(size), data(data)
    const std::size_t size;
    const T * data;
    const T operator[](std::size_t i) const
    { return data[i]; }
  };

  array_view<int> form::original_coefficient_positions() const
  {
    static const int data = { 0, 1, 2, 5 };
    return array_view<int>{4, data};
  }

  array_view<bool> integral::enabled_coefficients() const
  {
    static const bool data = { true, true, true, false, false, true };
    return array_view<bool>{6, data};
  }
"""


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
        # TODO: Change signature to
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
        #reference_dof_coordinates = L.FlattenedArray(L.Symbol("reference_dof_coordinates"), dims=(len(points), tdim))
        reference_dof_coordinates = L.Symbol("reference_dof_coordinates")

        # Loop indices
        j = L.Symbol("j")
        ip = L.Symbol("ip")

        # Reference coordinates
        dof_X = L.Symbol("dof_X")
        dof_X_values = [X[jj] for X in points for jj in range(tdim)]
        code = [
            L.ArrayDecl("static const double", dof_X,
                        (len(points) * tdim,), values=dof_X_values),
            L.ForRange(ip, 0, len(points), body=
                L.ForRange(j, 0, tdim, body=
                    L.AssignAdd(reference_dof_coordinates[ip*tdim + j], dof_X[ip*tdim + j]))),
            ]
        return code

    def evaluate_reference_basis(self, L, ir, parameters): # FIXME: NEW implement!
        """TODO: Add this signature to finite_element:
        evaluate_reference_basis(...)
        """
        data = ir["evaluate_basis"]
        from ffc.uflacs.backends.ufc.evaluatebasis import generate_evaluate_reference_basis
        return generate_evaluate_reference_basis(L, data)

    def evaluate_reference_basis_derivatives(self, L, ir, parameters): # FIXME: NEW implement!
        """TODO: Add this signature to finite_element:
        evaluate_reference_basis_derivatives(...)
        """
        #data = ir["evaluate_reference_basis_derivatives"]
        msg = "FIXME NOT IMPLEMENTED"
        return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])
