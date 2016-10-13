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
from ffc.uflacs.backends.ufc.utils import generate_return_new_switch


def affine_weights(dim): # FIXME: This is used where we still assume an affine mesh. Get rid of all places that use it.
    "Compute coefficents for mapping from reference to physical element"
    if dim == 1:
        return lambda X: (1.0 - X[0], X[0])
    elif dim == 2:
        return lambda X: (1.0 - X[0] - X[1], X[0], X[1])
    elif dim == 3:
        return lambda X: (1.0 - X[0] - X[1] - X[2], X[0], X[1], X[2])


class ufc_finite_element(ufc_generator):
    def __init__(self):
        ufc_generator.__init__(self, "finite_element")

    def topological_dimension(self, L, ir):
        tdim = ir["topological_dimension"]
        return L.Return(L.LiteralInt(tdim))

    def geometric_dimension(self, L, ir):
        gdim = ir["geometric_dimension"]
        return L.Return(L.LiteralInt(gdim))

    def degree(self, L, ir):
        degree = ir["degree"]
        return L.Return(L.LiteralInt(degree))

    def family(self, L, ir):
        family = ir["family"]
        return L.Return(L.LiteralString(family))

    def cell_shape(self, L, ir):
        name = ir["cell_shape"]
        return L.Return(L.Symbol(name))

    def space_dimension(self, L, ir):
        value = ir["space_dimension"]
        return L.Return(L.LiteralInt(value))

    def value_rank(self, L, ir):
        sh = ir["value_dimension"]
        return L.Return(L.LiteralInt(len(sh)))

    def value_size(self, L, ir):
        sh = ir["value_dimension"]
        return L.Return(L.LiteralInt(product(sh)))

    def value_dimension(self, L, ir):
        i = L.Symbol("i")
        sh = ir["value_dimension"]
        cases = [(L.LiteralInt(j), L.Return(L.LiteralInt(k))) for j, k in enumerate(sh)]
        default = L.Return(L.LiteralInt(0))
        return L.Switch(i, cases, default=default, autoscope=False, autobreak=False)

    def reference_value_rank(self, L, ir):
        sh = ir["reference_value_dimension"]
        return L.Return(L.LiteralInt(len(sh)))

    def reference_value_size(self, L, ir):
        sh = ir["reference_value_dimension"]
        return L.Return(L.LiteralInt(product(sh)))

    def reference_value_dimension(self, L, ir):
        i = L.Symbol("i")
        sh = ir["reference_value_dimension"]
        cases = [(L.LiteralInt(j), L.Return(L.LiteralInt(k))) for j, k in enumerate(sh)]
        default = L.Return(L.LiteralInt(0))
        return L.Switch(i, cases, default=default, autoscope=False, autobreak=False)

    def evaluate_reference_basis(self, L, ir): # FIXME: NEW implement!
        from ffc.uflacs.backends.ufc.evaluatebasis import generate_evaluate_reference_basis
        return generate_evaluate_reference_basis(ir["evaluate_reference_basis"])

    def evaluate_reference_basis_derivatives(self, L, ir): # FIXME: NEW implement!
        data = ir["evaluate_reference_basis_derivatives"]
        return "FIXME"

    def evaluate_basis(self, L, ir): # FIXME: Get rid of this
        data = ir["evaluate_basis"]
        return "FIXME"

    def evaluate_basis_derivatives(self, L, ir): # FIXME: Get rid of this
        # FIXME: port this, then translate into reference version
        data = ir["evaluate_basis_derivatives"]
        return "FIXME"

    def evaluate_dof(self, L, ir): # FIXME: Get rid of this
        # FIXME: port this, then translate into reference version
        data = ir["evaluate_dof"]
        return "FIXME"


    def evaluate_basis_all(self, L, ir):
        # FIXME: port this, then translate into reference version
        data = ir["evaluate_basis_all"]
        return "FIXME"

    def evaluate_basis_derivatives_all(self, L, ir):
        # FIXME: port this, then translate into reference version
        data = ir["evaluate_basis_derivatives_all"]
        return "FIXME"

    def evaluate_dofs(self, L, ir):
        """Generate code for evaluate_dofs."""
        # FIXME: port this, then translate into reference version
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

        return "FIXME" + ir["evaluate_dofs"]

    def interpolate_vertex_values(self, L, ir): # FIXME: port this
        # FIXME: port this, then translate into reference version
        return "FIXME" + ir["interpolate_vertex_values"]

    def _tabulate_dof_reference_coordinates(self, L, ir):
        """TODO: Add this signature to finite_element:
        /// Tabulate the reference coordinates of all dofs on a cell
        virtual void tabulate_dof_reference_coordinates(double * X) const = 0;
        """
        pass

    def tabulate_dof_coordinates(self, L, ir): # FIXME: port this
        # TODO: For a transition period, let finite_element and dofmap depend on a class affine_<cellname>_domain?
        # TODO: Call _tabulate_dof_reference_coordinates to tabulate X[ndofs][tdim],
        # then call affine_domain::compute_physical_coordinates(x, X, coordinate_dofs)

        ir = ir["tabulate_dof_coordinates"]

        # Raise error if tabulate_dof_coordinates is ill-defined
        if not ir:
            msg = "tabulate_dof_coordinates is not defined for this element"
            return L.Comment(msg) #L.Raise(msg) # TODO: Error handling

        # Extract coordinates and cell dimension
        gdim = ir["gdim"]
        tdim = ir["tdim"]
        points = ir["points"]

        # Output argument
        dof_coordinates = L.FlattenedArray(L.Symbol("dof_coordinates"), dims=(len(points), gdim))

        # Input argument
        coordinate_dofs = L.Symbol("coordinate_dofs")

        # Reference coordinates
        dof_reference_coordinates = L.Symbol("dof_reference_coordinates")
        dof_reference_coordinate_values = [X[j] for X in points for j in range(tdim)]

        # Loop indices
        i = L.Symbol("i")
        k = L.Symbol("k")
        ip = L.Symbol("ip")

        # Basis symbol
        phi = L.Symbol("phi")

        # TODO: This code assumes an affine coordinate field.
        #       Ok for now in here, this function must be removed anyway.
        # Create code for evaluating affine coordinate basis functions
        num_scalar_xdofs = tdim + 1
        cg1_basis = affine_weights(tdim)
        phi_values = [phi_comp for X in points for phi_comp in cg1_basis(X)]
        assert len(phi_values) == len(points) * num_scalar_xdofs

        code = [
            L.ArrayDecl("static const double", dof_reference_coordinates,
                        (len(points) * tdim,),
                        values=dof_reference_coordinate_values),
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
        return L.StatementList(code)

    def num_sub_elements(self, L, ir):
        n = ir["num_sub_elements"]
        return L.Return(L.LiteralInt(n))

    def create_sub_element(self, L, ir):
        i = L.Symbol("i")
        classnames = ir["create_sub_element"]
        return generate_return_new_switch(L, i, classnames, factory=ir["jit"])
