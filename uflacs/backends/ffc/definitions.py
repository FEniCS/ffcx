
from uflacs.backends.ffc.common import *

class FFCDefinitionsBackend(MultiFunction):
    """FFC specific cpp formatter class."""
    def __init__(self, ir, parameters):
        MultiFunction.__init__(self)

        # Store ir and parameters
        self.ir = ir
        self.parameters = parameters

        # Configure definitions behaviour
        self.physical_coordinates_known = self.ir["integral_type"] == "custom"

    def get_includes(self):
        "Return include statements to insert at top of file."
        includes = []
        return includes

    def initial(self):
        "Return code inserted at beginning of kernel."
        code = []
        return code

    def expr(self, t, mt, tabledata, access):
        error("Unhandled type {0}".format(type(t)))

    # === Generate code definitions ===

    def argument(self, t, mt, tabledata, access):
        code = []
        return code

    def coefficient(self, t, mt, tabledata, access):
        code = []
        if mt.terminal.is_cellwise_constant():
            # For a constant coefficient we reference the dofs directly, so no definition needed
            pass
        else:
            # No need to store basis function value in its own variable, just get table value directly
            uname, begin, end = tabledata
            entity = format_entity_name(self.ir["entitytype"], mt.restriction)

            iq = names.iq
            idof = names.ic

            dof = format_code(Sub(idof, begin))
            table_access = ArrayAccess(uname, (entity, iq, dof))

            dof_access = generate_coefficient_dof_access(mt.terminal, idof)

            prod = Mul(dof_access, table_access)
            body = [AssignAdd(access, prod)]

            # Loop to accumulate linear combination of dofs and tables
            code += [VariableDecl("double", access, "0.0")]
            code += [ForRange(idof, begin, end, body=body)]

        return code

    def quadrature_weight(self, e, mt, tabledata, access):
        code = []
        return code

    def spatial_coordinate(self, e, mt, tabledata, access):
        """Return definition code for the physical spatial coordinates.

        If physical coordinates are given:
          No definition needed.

        If reference coordinates are given:
          x = sum_k xdof_k xphi_k(X)

        If reference facet coordinates are given:
          x = sum_k xdof_k xphi_k(Xf)
        """
        code = []

        if self.physical_coordinates_known:
            pass
        else:
            ffc_assert(mt.terminal.domain().coordinates() is None,
                          "Assuming coefficient field symbolically inserted before this point.")
            # Reference coordinates are known, no coordinate field, so we compute
            # this component as linear combination of vertex_coordinates "dofs" and table

            cell = mt.terminal.domain().cell()
            gdim = cell.geometric_dimension()
            num_vertices = cell.topological_dimension() + 1 # FIXME: Get from cell!

            uname, begin, end = tabledata

            # access here is e.g. x0, component 0 of x

            ffc_assert(0 <= begin <= end <= num_vertices*gdim,
                          "Assuming linear element for affine simplices here.")
            entity = format_entity_name(self.ir["entitytype"], mt.restriction)
            iq = names.iq

            if 0: # Generated loop version:
                vertex = names.ic
                table_access = ArrayAccess(uname, (entity, iq, vertex))
                dof_access = generate_domain_dof_access(num_vertices, gdim, vertex,
                                                        mt.flat_component, mt.restriction)
                prod = Mul(dof_access, table_access)

                # Loop to accumulate linear combination of dofs and tables
                code += [VariableDecl("double", access, "0.0")]
                code += [ForRange(vertex, begin, end, body=[AssignAdd(access, prod)])]

            else: # Inlined version (we know this is bounded by a small number)
                dof_access = generate_domain_dofs_access(num_vertices, gdim, mt.restriction)
                prods = []
                for idof in range(begin, end):
                    table_access = ArrayAccess(uname, (entity, iq, Sub(idof, begin)))
                    prods += [Mul(dof_access[idof], table_access)]

                # Inlined loop to accumulate linear combination of dofs and tables
                code += [VariableDecl("const double", access, Sum(prods))]

        return code

    def cell_coordinate(self, e, mt, tabledata, access):
        """Return definition code for the reference spatial coordinates.

        If reference coordinates are given:
          No definition needed.

        If physical coordinates are given and domain is affine:
          X = K*(x-x0)
        This is inserted symbolically.

        If physical coordinates are given and domain is non- affine:
          Not currently supported.
        """
        code = []
        return code

    def jacobian(self, e, mt, tabledata, access):
        """Return definition code for the Jacobian of x(X).

        J = sum_k xdof_k grad_X xphi_k(X)
        """
        code = []

        if self.physical_coordinates_known:
            pass
        else:
            ffc_assert(mt.terminal.domain().coordinates() is None,
                          "Assuming coefficient field symbolically inserted before this point.")
            # Reference coordinates are known, no coordinate field, so we compute
            # this component as linear combination of vertex_coordinates "dofs" and table

            cell = mt.terminal.domain().cell()
            gdim = cell.geometric_dimension()
            num_vertices = cell.topological_dimension() + 1 # FIXME: Get from cell!

            uname, begin, end = tabledata

            # access here is e.g. J_0, component 0 of J

            ffc_assert(0 <= (end-begin) <= num_vertices,
                          "Assuming linear element for affine simplices here.")
            entity = format_entity_name(self.ir["entitytype"], mt.restriction)
            vertex = names.ic
            iq = 0

            if 0: # Generated loop version:
                table_access = ArrayAccess(uname, iq, (entity, vertex))
                dof_access = generate_domain_dof_access(num_vertices, gdim, vertex,
                                                        mt.flat_component, mt.restriction)
                prod = Mul(dof_access, table_access)

                # Loop to accumulate linear combination of dofs and tables
                code += [VariableDecl("double", access, "0.0")]
                code += [ForRange(vertex, begin, end, body=[AssignAdd(access, prod)])]

            else: # Inlined version:
                prods = []
                dof_access = generate_domain_dofs_access(num_vertices, gdim, mt.restriction)
                for idof in range(begin, end):
                    table_access = ArrayAccess(uname, (entity, 0, Sub(idof, begin)))
                    prods += [Mul(dof_access[idof], table_access)]

                # Inlined loop to accumulate linear combination of dofs and tables
                code += [VariableDecl("const double", access, Sum(prods))]

        return code

    def cell_facet_jacobian(self, e, mt, tabledata, access):
        # Constant table defined in ufc_geometry.h
        code = []
        return code

    def cell_edge_vectors(self, e, mt, tabledata, access):
        # Constant table defined in ufc_geometry.h
        code = []
        return code

    def facet_edge_vectors(self, e, mt, tabledata):
        # Constant table defined in ufc_geometry.h
        code = []
        return code

    def cell_orientation(self, e, mt, tabledata, access):
        # Computed or constant table defined in ufc_geometry.h
        code = []
        return code

    def facet_orientation(self, e, mt, tabledata, access):
        # Constant table defined in ufc_geometry.h
        code = []
        return code

    def facet_normal(self, e, mt, tabledata, access):
        error("Expecting {0} to be replaced with lower level types in symbolic preprocessing.".format(type(e)))

    def cell_normal(self, e, mt, tabledata, access):
        error("Expecting {0} to be replaced with lower level types in symbolic preprocessing.".format(type(e)))

    def jacobian_inverse(self, e, mt, tabledata, access):
        error("Expecting {0} to be replaced with lower level types in symbolic preprocessing.".format(type(e)))

    def jacobian_determinant(self, e, mt, tabledata, access):
        error("Expecting {0} to be replaced with lower level types in symbolic preprocessing.".format(type(e)))

    def facet_jacobian(self, e, mt, tabledata, access):
        error("Expecting {0} to be replaced with lower level types in symbolic preprocessing.".format(type(e)))

    def facet_jacobian_inverse(self, e, mt, tabledata, access):
        error("Expecting {0} to be replaced with lower level types in symbolic preprocessing.".format(type(e)))

    def facet_jacobian_determinant(self, e, mt, tabledata, access):
        error("Expecting {0} to be replaced with lower level types in symbolic preprocessing.".format(type(e)))
