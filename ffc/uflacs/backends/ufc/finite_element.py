# -*- coding: utf-8 -*-
# Copyright (C) 2009-2017 Anders Logg and Martin Sandve Alnæs
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


# Note: Much of the code in this file is a direct translation
# from the old implementation in FFC, although some improvements
# have been made to the generated code.


from collections import defaultdict
import numpy

from ufl import product
from ffc.uflacs.backends.ufc.generator import ufc_generator
from ffc.uflacs.backends.ufc.utils import generate_return_new_switch, generate_return_int_switch, generate_error

from ffc.uflacs.elementtables import clamp_table_small_numbers
from ffc.uflacs.backends.ufc.evaluatebasis import generate_evaluate_reference_basis
from ffc.uflacs.backends.ufc.evaluatebasis import _generate_compute_basisvalues, tabulate_coefficients
from ffc.uflacs.backends.ufc.evaluatebasisderivatives import generate_evaluate_basis_derivatives_all
from ffc.uflacs.backends.ufc.evaluatebasisderivatives import generate_evaluate_basis_derivatives
from ffc.uflacs.backends.ufc.evalderivs import generate_evaluate_reference_basis_derivatives
from ffc.uflacs.backends.ufc.evalderivs import _generate_combinations
from ffc.uflacs.backends.ufc.evaluatedof import generate_evaluate_dof, generate_evaluate_dofs, reference_to_physical_map

from ffc.uflacs.backends.ufc.jacobian import jacobian, inverse_jacobian, orientation, fiat_coordinate_mapping, _mapping_transform

index_type = "std::size_t"

def compute_basis_values(L, data, dof_data):
    basisvalues = L.Symbol("basisvalues")
    Y = L.Symbol("Y")
    element_cellname = data["cellname"]
    embedded_degree = dof_data["embedded_degree"]
    num_members = dof_data["num_expansion_members"]
    return _generate_compute_basisvalues(L, basisvalues, Y, element_cellname,
                                         embedded_degree, num_members)


def compute_values(L, data, dof_data):
    """This function computes the value of the basisfunction as the dot product
    of the coefficients and basisvalues."""

    # Initialise return code.
    code = [L.Comment("Compute value(s)")]

    # Get dof data.
    num_components = dof_data["num_components"]
    physical_offset = dof_data["physical_offset"]

    basisvalues = L.Symbol("basisvalues")
    values = L.Symbol("values")
    r = L.Symbol("r")
    lines = []
    if data["reference_value_size"] != 1:
        # Loop number of components.
        for i in range(num_components):
            coefficients = L.Symbol("coefficients%d" % i)
            lines += [L.AssignAdd(values[i + physical_offset], coefficients[r]*basisvalues[r])]
    else:
        coefficients = L.Symbol("coefficients0")
        lines = [L.AssignAdd(L.Dereference(values), coefficients[r]*basisvalues[r])]

    # Get number of members of the expansion set and generate loop.
    num_mem = dof_data["num_expansion_members"]
    code += [L.ForRange(r, 0, num_mem, index_type=index_type, body=lines)]

    code += _mapping_transform(L, data, dof_data, values, physical_offset)

    return code

def generate_element_mapping(mapping, i, num_reference_components, tdim, gdim, J, detJ, K):
    # Select transformation to apply
    if mapping == "affine":
        assert num_reference_components == 1
        num_physical_components = 1
        M_scale = 1
        M_row = [1]  # M_row[0] == 1
    elif mapping == "contravariant piola":
        assert num_reference_components == tdim
        num_physical_components = gdim
        M_scale = 1.0 / detJ
        M_row = [J[i, jj] for jj in range(tdim)]
    elif mapping == "covariant piola":
        assert num_reference_components == tdim
        num_physical_components = gdim
        M_scale = 1.0
        M_row = [K[jj, i] for jj in range(tdim)]
    elif mapping == "double covariant piola":
        assert num_reference_components == tdim**2
        num_physical_components = gdim**2
        # g_il = K_ji G_jk K_kl = K_ji K_kl G_jk
        i0 = i // tdim  # i in the line above
        i1 = i % tdim   # l ...
        M_scale = 1.0
        M_row = [K[jj,i0]*K[kk,i1] for jj in range(tdim) for kk in range(tdim)]
    elif mapping == "double contravariant piola":
        assert num_reference_components == tdim**2
        num_physical_components = gdim**2
        # g_il = (det J)^(-2) Jij G_jk Jlk = (det J)^(-2) Jij Jlk G_jk
        i0 = i // tdim  # i in the line above
        i1 = i % tdim   # l ...
        M_scale = 1.0 / (detJ*detJ)
        M_row = [J[i0,jj]*J[i1,kk] for jj in range(tdim) for kk in range(tdim)]
    else:
        error("Unknown mapping: %s" % mapping)
    return M_scale, M_row, num_physical_components


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

    def evaluate_basis(self, L, ir, parameters):
        data = ir["evaluate_basis"]

        # FIXME: does this make sense?
        if not data:
            msg = "evaluate_basis is not defined for this element"
            return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])

        # Handle unsupported elements.
        if isinstance(data, str):
            msg = "evaluate_basis: %s" % data
            return [generate_error(L, msg, parameters["convert_exceptions_to_warnings"])]

        # Get the element cell name and geometric dimension.
        element_cellname = data["cellname"]
        gdim = data["geometric_dimension"]
        tdim = data["topological_dimension"]

        # Generate run time code to evaluate an element basisfunction
        # at an arbitrary point. The value(s) of the basisfunction
        # is/are computed as in FIAT as the dot product of the
        # coefficients (computed at compile time) and basisvalues
        # which are dependent on the coordinate and thus have to be
        # computed at run time.

        # The function should work for all elements supported by FIAT,
        # but it remains untested for tensor valued elements.

        # Get code snippets for Jacobian, Inverse of Jacobian and
        # mapping of coordinates from physical element to the FIAT
        # reference element.

        cm = L.Symbol("cm")
        X = L.Symbol("X")
        code = [L.ArrayDecl("double", X, (tdim), values=0)]

        J = L.Symbol("J")
        code += [L.ArrayDecl("double", J, (gdim*tdim,))]

        detJ = L.Symbol("detJ")
        code += [L.VariableDecl("double", detJ)]

        K = L.Symbol("K")
        code += [L.ArrayDecl("double", K, (gdim*tdim,))]

        x = L.Symbol("x")
        coordinate_dofs = L.Symbol("coordinate_dofs")
        cell_orientation = L.Symbol("cell_orientation")

        no_cm_code = [  L.Call("compute_jacobian_"+element_cellname+"_"+str(gdim)+"d",
                               (J, coordinate_dofs)),
                        L.Call("compute_jacobian_inverse_"+element_cellname+"_"+str(gdim)+"d",
                               (K, detJ, J))]

        if data["needs_oriented"] and tdim != gdim:
            no_cm_code += orientation(L)

        if any((d["embedded_degree"] > 0) for d in data["dofs_data"]):
            k = L.Symbol("k")
            Y = L.Symbol("Y")
            no_cm_code += fiat_coordinate_mapping(L, element_cellname, gdim)
            if element_cellname in ('interval', 'triangle', 'tetrahedron'):
                no_cm_code += [L.Comment("Map to FFC reference coordinate"),
                               L.ForRange(k, 0, tdim, index_type=index_type, body=[L.Assign(X[k], (Y[k] + 1.0)/2.0)])]
            else:
                no_cm_code += [L.ForRange(k, 0, tdim, index_type=index_type, body=[L.Assign(X[k], Y[k])])]

        code += [L.If(cm, L.Call("cm->compute_reference_geometry",
                                (X, J, L.AddressOf(detJ), K, 1, x, coordinate_dofs, cell_orientation))),
                 L.Else(no_cm_code)]

        reference_value_size = data["reference_value_size"]
        num_dofs = len(data["dofs_data"])
        ref_values = L.Symbol("ref_values")
        code += [L.Comment("Evaluate basis on reference element"),
                 L.ArrayDecl("double", ref_values, num_dofs*reference_value_size),
                 L.Call("evaluate_reference_basis",(ref_values, 1, X))]

        physical_value_size = data["physical_value_size"]
        physical_values = L.Symbol("physical_values")
        i = L.Symbol("i")
        k = L.Symbol("k")
        values = L.Symbol("values")
        code += [L.Comment("Push forward"),
                 L.ArrayDecl("double", physical_values, num_dofs*physical_value_size),
                 L.Call("transform_reference_basis_derivatives",(physical_values, 0, 1,
                                                                 ref_values, X, J, L.AddressOf(detJ), K,
                                                                 cell_orientation)),
                 L.ForRange(k, 0, physical_value_size, index_type=index_type, body=[
                     L.Assign(values[k], physical_values[physical_value_size*i + k])
                 ])]

        return code

    def evaluate_basis_all(self, L, ir, parameters):
        data=ir["evaluate_basis"]

        # Handle unsupported elements.
        if isinstance(data, str):
            msg = "evaluate_basis_all: %s" % data
            return [generate_error(L, msg, parameters["convert_exceptions_to_warnings"])]

        physical_value_size = data["physical_value_size"]
        space_dimension = data["space_dimension"]

        x = L.Symbol("x")
        coordinate_dofs = L.Symbol("coordinate_dofs")
        cell_orientation = L.Symbol("cell_orientation")
        values = L.Symbol("values")

        # Special case where space dimension is one (constant elements).
        if space_dimension == 1:
            code = [L.Comment("Element is constant, calling evaluate_basis."),
                    L.Call("evaluate_basis",
                           (0, values, x, coordinate_dofs, cell_orientation))]
            return code

        r = L.Symbol("r")
        dof_values = L.Symbol("dof_values")
        if physical_value_size == 1:
            code = [ L.Comment("Helper variable to hold value of a single dof."),
                     L.VariableDecl("double", dof_values, 0.0),
                     L.Comment("Loop dofs and call evaluate_basis"),
                     L.ForRange(r, 0, space_dimension, index_type=index_type,
                                body=[L.Call("evaluate_basis",
                                             (r, L.AddressOf(dof_values), x,
                                              coordinate_dofs, cell_orientation)),
                                      L.Assign(values[r], dof_values)]
                               )
                   ]
        else:
            s = L.Symbol("s")
            code = [L.Comment("Helper variable to hold values of a single dof."),
                    L.ArrayDecl("double", dof_values, physical_value_size, 0.0),
                    L.Comment("Loop dofs and call evaluate_basis"),
                    L.ForRange(r, 0, space_dimension, index_type=index_type,
                               body=[L.Call("evaluate_basis",
                                             (r, dof_values, x,
                                              coordinate_dofs, cell_orientation)),
                                     L.ForRange(s, 0, physical_value_size,
                                                index_type=index_type,
                                                body=[L.Assign(values[r*physical_value_size+s], dof_values[s])])
                                    ]
                              )
                   ]

        return code

    def evaluate_basis_derivatives(self, L, ir, parameters):
        # FIXME: Get rid of this
        return generate_evaluate_basis_derivatives(L, ir["evaluate_basis"])

    def evaluate_basis_derivatives_all(self, L, ir, parameters):
        return generate_evaluate_basis_derivatives_all(L, ir["evaluate_basis"])

        """
        // Legacy version:
        evaluate_basis_derivatives_all(std::size_t n,
                                       double * values,
                                       const double * x,
                                       const double * coordinate_dofs,
                                       int cell_orientation)
        // Suggestion for new version:
        new_evaluate_basis_derivatives(double * values,
                                       std::size_t order,
                                       std::size_t num_points,
                                       const double * x,
                                       const double * coordinate_dofs,
                                       int cell_orientation,
                                       const ufc::coordinate_mapping * cm)
        """

        # TODO: This is a refactoring step to allow rewriting code
        # generation to use coordinate_mapping in one stage, before
        # making it available as an argument from dolfin in the next stage.
        affine_coordinate_mapping_classname = ir["affine_coordinate_mapping_classname"]

        # Output arguments:
        values = L.Symbol("values")

        # Input arguments:
        #order = L.Symbol("order")
        order = L.Symbol("n")
        x = L.Symbol("x")
        coordinate_dofs = L.Symbol("coordinate_dofs")
        cell_orientation = L.Symbol("cell_orientation")

        # Internal variables:
        #num_points = L.Symbol("num_points")
        num_points = 1  # Always 1 in legacy API
        reference_values = L.Symbol("reference_values")
        X = L.Symbol("X")
        J = L.Symbol("J")
        detJ = L.Symbol("detJ")
        K = L.Symbol("K")
        ip = L.Symbol("ip")

        gdim = ir["geometric_dimension"]
        tdim = ir["topological_dimension"]

        code = [
            # Create local affine coordinate mapping object
            # TODO: Get this as input instead to support non-affine
            L.VariableDecl(affine_coordinate_mapping_classname, "cm"),
            L.ForRange(ip, 0, num_points, index_type=index_type, body=[
                L.ArrayDecl("double", X, (tdim,)),
                L.ArrayDecl("double", J, (gdim*tdim,)),
                L.ArrayDecl("double", detJ, (1,)),
                L.ArrayDecl("double", K, (tdim*gdim,)),
                L.Call("cm.compute_reference_geometry",
                       (X, J, detJ, K, num_points, x, coordinate_dofs, cell_orientation)),
                L.Call("evaluate_reference_basis_derivatives",
                       (reference_values, order, num_points, X)),
                L.Call("transform_reference_basis_derivatives",
                       (values, order, num_points, reference_values, X, J, detJ, K, cell_orientation)),
            ])
        ]
        return code

    def evaluate_dof(self, L, ir, parameters):
        return generate_evaluate_dof(L, ir["evaluate_dof"])


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
        # FIXME: translate into reference version
        return generate_evaluate_dofs(L, ir["evaluate_dof"])

    def interpolate_vertex_values(self, L, ir, parameters):
        irdata = ir["interpolate_vertex_values"]

        # Raise error if interpolate_vertex_values is ill-defined
        if not irdata:
            msg = "interpolate_vertex_values is not defined for this element"
            return [generate_error(L, msg, parameters["convert_exceptions_to_warnings"])]

        # Handle unsupported elements.
        if isinstance(irdata, str):
            msg = "interpolate_vertex_values: %s" % irdata
            return [generate_error(L, msg, parameters["convert_exceptions_to_warnings"])]

        # Add code for Jacobian if necessary
        code = []
        gdim = irdata["geometric_dimension"]
        tdim = irdata["topological_dimension"]
        cell_shape = ir["cell_shape"]
        if irdata["needs_jacobian"]:
            code += jacobian(L, gdim, tdim, cell_shape)
            code += inverse_jacobian(L, gdim, tdim, cell_shape)
            if irdata["needs_oriented"] and tdim != gdim:
                code += orientation(L)

        # Compute total value dimension for (mixed) element
        total_dim = irdata["physical_value_size"]

        # Generate code for each element
        value_offset = 0
        space_offset = 0
        for data in irdata["element_data"]:
            # Add vertex interpolation for this element
            code += [L.Comment("Evaluate function and change variables")]

            # Extract vertex values for all basis functions
            vertex_values = data["basis_values"]
            value_size = data["physical_value_size"]
            space_dim = data["space_dim"]
            mapping = data["mapping"]

            J = L.Symbol("J")
            J = L.FlattenedArray(J, dims=(gdim, tdim))
            detJ = L.Symbol("detJ")
            K = L.Symbol("K")
            K = L.FlattenedArray(K, dims=(tdim, gdim))

            # Create code for each value dimension:
            for k in range(value_size):
                # Create code for each vertex x_j
                for (j, values_at_vertex) in enumerate(vertex_values):

                    if value_size == 1:
                        values_at_vertex = [values_at_vertex]

                    values = clamp_table_small_numbers(values_at_vertex)

                    # Map basis functions using appropriate mapping
                    # FIXME: sort out all non-affine mappings and make into a function
                    # components = change_of_variables(values_at_vertex, k)

                    w = []
                    if mapping == 'affine':
                        w = values[k]
                    elif mapping == 'contravariant piola':
                        for index in range(space_dim):
                            w += [sum(J[k, p]*values[p][index]
                                      for p in range(tdim))/detJ]
                    elif mapping == 'covariant piola':
                        for index in range(space_dim):
                            w += [sum(K[p, k]*values[p][index]
                                      for p in range(tdim))]
                    elif mapping == 'double covariant piola':
                        for index in range(space_dim):
                            w += [sum(K[p, k//tdim]*values[p][q][index]*K[q, k % tdim]
                                      for q in range(tdim) for p in range(tdim))]
                    elif mapping == 'double contravariant piola':
                        for index in range(space_dim):
                            w += [sum(J[k//tdim, p]*values[p][q][index]*J[k % tdim, q]
                                      for q in range(tdim) for p in range(tdim))/(detJ*detJ)]
                    else:
                        error("Unknown mapping: %s" % mapping)

                    # Contract coefficients and basis functions
                    dof_values = L.Symbol("dof_values")
                    dof_list = [dof_values[i + space_offset] for i in range(space_dim)]
                    value = sum(p*q for (p, q) in zip(dof_list, w))

                    # Assign value to correct vertex
                    index = j * total_dim + (k + value_offset)
                    v_values = L.Symbol("vertex_values")
                    code += [L.Assign(v_values[index], value)]

            # Update offsets for value- and space dimension
            value_offset += data["physical_value_size"]
            space_offset += data["space_dim"]

        return code

    def tabulate_dof_coordinates(self, L, ir, parameters):
        ir = ir["tabulate_dof_coordinates"]

        # Raise error if tabulate_dof_coordinates is ill-defined
        if not ir:
            msg = "tabulate_dof_coordinates is not defined for this element"
            return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])

        # Extract coordinates and cell dimension
        gdim = ir["gdim"]
        tdim = ir["tdim"]
        points = ir["points"]

        # Extract cellshape
        cell_shape = ir["cell_shape"]

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

        # TODO: Get rid of all places that use reference_to_physical_map, it is restricted to a basis of degree 1
        # Create code for evaluating coordinate mapping
        num_scalar_xdofs = _num_vertices(cell_shape)
        cg1_basis = reference_to_physical_map(cell_shape)
        phi_values = numpy.asarray([phi_comp for X in points for phi_comp in cg1_basis(X)])
        assert len(phi_values) == len(points) * num_scalar_xdofs

        # TODO: Use precision parameter here
        phi_values = clamp_table_small_numbers(phi_values)

        code = [
            L.Assign(
                dof_coordinates[ip][i],
                sum(phi_values[ip*num_scalar_xdofs + k] * coordinate_dofs[gdim*k + i]
                    for k in range(num_scalar_xdofs))
            )
            for ip in range(len(points))
            for i in range(gdim)
        ]

        # FIXME: This code assumes an affine coordinate field.
        #        To get around that limitation, make this function take another argument
        #            const ufc::coordinate_mapping * cm
        #        and generate code like this:
        """
        index_type X[tdim*num_dofs];
        tabulate_dof_coordinates(X);
        cm->compute_physical_coordinates(x, X, coordinate_dofs);
        """

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
        if isinstance(data, str):
            msg = "evaluate_reference_basis: %s" % data
            return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])

        return generate_evaluate_reference_basis(L, data, parameters)

    def evaluate_reference_basis_derivatives(self, L, ir, parameters):
        data = ir["evaluate_basis"]
        if isinstance(data, str):
            msg = "evaluate_reference_basis_derivatives: %s" % data
            return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])

        return generate_evaluate_reference_basis_derivatives(L, data, parameters)

    def transform_reference_basis_derivatives(self, L, ir, parameters):
        data = ir["evaluate_basis"]
        if isinstance(data, str):
            msg = "transform_reference_basis_derivatives: %s" % data
            return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])

        # Get some known dimensions
        #element_cellname = data["cellname"]
        gdim = data["geometric_dimension"]
        tdim = data["topological_dimension"]
        max_degree = data["max_degree"]
        reference_value_size = data["reference_value_size"]
        physical_value_size = data["physical_value_size"]
        num_dofs = len(data["dofs_data"])

        max_g_d = gdim**max_degree
        max_t_d = tdim**max_degree

        # Output arguments
        values_symbol = L.Symbol("values")

        # Input arguments
        order = L.Symbol("order")
        num_points = L.Symbol("num_points")  # FIXME: Currently assuming 1 point?
        reference_values = L.Symbol("reference_values")
        J = L.Symbol("J")
        detJ = L.Symbol("detJ")
        K = L.Symbol("K")

        # Internal variables
        transform = L.Symbol("transform")

        # Indices, I've tried to use these for a consistent purpose
        ip = L.Symbol("ip") # point
        i = L.Symbol("i")   # physical component
        j = L.Symbol("j")   # reference component
        k = L.Symbol("k")   # order
        r = L.Symbol("r")   # physical derivative number
        s = L.Symbol("s")   # reference derivative number
        d = L.Symbol("d")   # dof

        combinations_code = []
        if max_degree == 0:
            # Don't need combinations
            num_derivatives_t = 1  # TODO: I think this is the right thing to do to make this still work for order=0?
            num_derivatives_g = 1
        elif tdim == gdim:
            num_derivatives_t = L.Symbol("num_derivatives")
            num_derivatives_g = num_derivatives_t
            combinations_code += [
                L.VariableDecl("const " + index_type, num_derivatives_t,
                               L.Call("std::pow", (tdim, order))),
            ]

            # Add array declarations of combinations
            combinations_code_t, combinations_t = _generate_combinations(L, tdim, max_degree, order, num_derivatives_t)
            combinations_code += combinations_code_t
            combinations_g = combinations_t
        else:
            num_derivatives_t = L.Symbol("num_derivatives_t")
            num_derivatives_g = L.Symbol("num_derivatives_g")
            combinations_code += [
                L.VariableDecl("const " + index_type, num_derivatives_t,
                               L.Call("std::pow", (tdim, order))),
                L.VariableDecl("const " + index_type, num_derivatives_g,
                               L.Call("std::pow", (gdim, order))),
            ]
            # Add array declarations of combinations
            combinations_code_t, combinations_t = _generate_combinations(L, tdim, max_degree, order, num_derivatives_t, suffix="_t")
            combinations_code_g, combinations_g = _generate_combinations(L, gdim, max_degree, order, num_derivatives_g, suffix="_g")
            combinations_code += combinations_code_t
            combinations_code += combinations_code_g

        # Define expected dimensions of argument arrays
        J = L.FlattenedArray(J, dims=(num_points, gdim, tdim))
        detJ = L.FlattenedArray(detJ, dims=(num_points,))
        K = L.FlattenedArray(K, dims=(num_points, tdim, gdim))

        values = L.FlattenedArray(values_symbol,
            dims=(num_points, num_dofs, num_derivatives_g, physical_value_size))
        reference_values = L.FlattenedArray(reference_values,
            dims=(num_points, num_dofs, num_derivatives_t, reference_value_size))

        # Generate code to compute the derivative transform matrix
        transform_matrix_code = [
            # Initialize transform matrix to all 1.0
            L.ArrayDecl("double", transform, (max_g_d, max_t_d)),
            L.ForRanges(
                (r, 0, num_derivatives_g),
                (s, 0, num_derivatives_t),
                index_type=index_type,
                body=L.Assign(transform[r, s], 1.0)
            ),
            ]
        if max_degree > 0:
            transform_matrix_code += [
                # Compute transform matrix entries, each a product of K entries
                L.ForRanges(
                    (r, 0, num_derivatives_g),
                    (s, 0, num_derivatives_t),
                    (k, 0, order),
                    index_type=index_type,
                    body=L.AssignMul(transform[r, s],
                                     K[ip, combinations_t[s, k], combinations_g[r, k]])
                ),
            ]

        # Initialize values to 0, will be added to inside loops
        values_init_code = [
            L.MemZero(values_symbol, num_points * num_dofs * num_derivatives_g * physical_value_size),
            ]

        # Make offsets available in generated code
        reference_offsets = L.Symbol("reference_offsets")
        physical_offsets = L.Symbol("physical_offsets")
        dof_attributes_code = [
            L.ArrayDecl("const " + index_type, reference_offsets, (num_dofs,),
                        values=[dof_data["reference_offset"] for dof_data in data["dofs_data"]]),
            L.ArrayDecl("const " + index_type, physical_offsets, (num_dofs,),
                        values=[dof_data["physical_offset"] for dof_data in data["dofs_data"]]),
            ]

        # Build dof lists for each mapping type
        mapping_dofs = defaultdict(list)
        for idof, dof_data in enumerate(data["dofs_data"]):
            mapping_dofs[dof_data["mapping"]].append(idof)

        # Generate code for each mapping type
        d = L.Symbol("d")
        transform_apply_code = []
        for mapping in sorted(mapping_dofs):
            # Get list of dofs using this mapping
            idofs = mapping_dofs[mapping]

            # Select iteration approach over dofs
            if idofs == list(range(idofs[0], idofs[-1]+1)):
                # Contiguous
                dofrange = (d, idofs[0], idofs[-1]+1)
                idof = d
            else:
                # Stored const array of dof indices
                idofs_symbol = L.Symbol("%s_dofs" % mapping.replace(" ", "_"))
                dof_attributes_code += [
                    L.ArrayDecl("const " + index_type, idofs_symbol,
                                (len(idofs),), values=idofs),
                ]
                dofrange = (d, 0, len(idofs))
                idof = idofs_symbol[d]

            # NB! Array access to offsets, these are not Python integers
            reference_offset = reference_offsets[idof]
            physical_offset = physical_offsets[idof]

            # How many components does each basis function with this mapping have?
            # This should be uniform, i.e. there should be only one element in this set:
            num_reference_components, = set(data["dofs_data"][i]["num_components"] for i in idofs)

            M_scale, M_row, num_physical_components = generate_element_mapping(
                mapping, i,
                num_reference_components, tdim, gdim,
                J[ip], detJ[ip], K[ip]
            )

#            transform_apply_body = [
#                L.AssignAdd(values[ip, idof, r, physical_offset + k],
#                            transform[r, s] * reference_values[ip, idof, s, reference_offset + k])
#                for k in range(num_physical_components)
#            ]

            msg = "Using %s transform to map values back to the physical element." % mapping.replace("piola", "Piola")

            mapped_value = L.Symbol("mapped_value")
            transform_apply_code += [
                L.ForRanges(
                    dofrange,
                    (s, 0, num_derivatives_t),
                    (i, 0, num_physical_components),
                    index_type=index_type, body=[
                        # Unrolled application of mapping to one physical component,
                        # for affine this automatically reduces to
                        #   mapped_value = reference_values[..., reference_offset]
                        L.Comment(msg),
                        L.VariableDecl("const double", mapped_value,
                                       M_scale * sum(M_row[jj] * reference_values[ip, idof, s, reference_offset + jj]
                                                     for jj in range(num_reference_components))),
                        # Apply derivative transformation, for order=0 this reduces to
                        # values[ip,idof,0,physical_offset+i] = transform[0,0]*mapped_value
                        L.Comment("Mapping derivatives back to the physical element"),
                        L.ForRanges(
                            (r, 0, num_derivatives_g),
                            index_type=index_type, body=[
                                L.AssignAdd(values[ip, idof, r, physical_offset + i],
                                            transform[r, s] * mapped_value)
                        ])
                ])
            ]

        # Transform for each point
        point_loop_code = [
            L.ForRange(ip, 0, num_points, index_type=index_type, body=(
                transform_matrix_code
                + transform_apply_code
            ))
        ]

        # Join code
        code = (
            combinations_code
            + values_init_code
            + dof_attributes_code
            + point_loop_code
        )
        return code


def _num_vertices(cell_shape):
    """Returns number of vertices for a given cell shape."""

    num_vertices_dict = {"interval": 2, "triangle": 3, "tetrahedron": 4, "quadrilateral": 4, "hexahedron": 8}
    return num_vertices_dict[cell_shape]
