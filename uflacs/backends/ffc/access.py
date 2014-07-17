
from uflacs.backends.ffc.common import *

class FFCAccessBackend(MultiFunction):
    """FFC specific cpp formatter class."""
    def __init__(self, ir, parameters):
        MultiFunction.__init__(self)

        # Store ir and parameters
        self.ir = ir
        self.parameters = parameters

        # Configure definitions behaviour
        self.physical_coordinates_known = self.ir["integral_type"] == "quadrature"

        self._current_num_points, = iterkeys(self.ir["uflacs"]["expr_ir"]) # TODO: Hack! Assuming single quadrature rule here.

    def precision_float(self, f):
        # Use ufl utility to control float formatting precision, same as ffc quadrature mode uses
        return ufl.constantvalue.format_float(f)

    def get_includes(self):
        "Return include statements to insert at top of file."
        includes = []
        return includes

    # === Access to names of quantities not among the symbolic UFL types ===

    def weights_array_name(self, num_points):
        return "{0}{1}".format(names.weights, num_points)

    def points_array_name(self, num_points):
        return "{0}{1}".format(names.points, num_points)

    def quadrature_loop_index(self):
        return names.iq

    def argument_loop_index(self, iarg):
        return "{name}{num}".format(name=names.ia, num=iarg)

    def element_tensor_name(self):
        return names.A

    def element_tensor_entry(self, index):
        return ArrayAccess(names.A, index)

    # === Multifunction handlers for all modified terminal types, basic C++ types are covered by base class ===

    def expr(self, e, mt, tabledata):
        error("Missing handler for type {0}.".format(e._uflclass.__name__))

    def argument(self, e, mt, tabledata):

        # Expecting only local derivatives and values here
        assert not mt.global_derivatives
        #assert mt.global_component is None

        # No need to store basis function value in its own variable, just get table value directly
        uname, begin, end = tabledata
        entity = format_entity_name(self.ir["entitytype"], mt.restriction)

        iq = self.quadrature_loop_index()
        idof = self.argument_loop_index(mt.terminal.number())

        dof = format_code(Sub(idof, begin))
        access = ArrayAccess(uname, (entity, iq, dof))
        return access

    def coefficient(self, e, mt, tabledata):
        t = mt.terminal
        if t.is_cellwise_constant():
            access = self._constant_coefficient(e, mt, tabledata)
        else:
            access = self._varying_coefficient(e, mt, tabledata)
        return access

    def _constant_coefficient(self, e, mt, tabledata):
        # Map component to flat index
        vi2si, si2vi = build_component_numbering(mt.terminal.shape(), mt.terminal.element().symmetry())
        num_flat_components = len(si2vi)
        ffc_assert(mt.flat_component == vi2si[mt.component], "Incompatible component flattening!")

        # Offset index if on second cell in interior facet integral
        if mt.restriction == "-": # TODO: Get the notion that '-' is the second cell from a central definition?
            idof = mt.flat_component + len(si2vi)
        else:
            idof = mt.flat_component

        # Return direct reference to dof array
        return generate_coefficient_dof_access(mt.terminal, idof)

    def _varying_coefficient(self, e, mt, tabledata):
        # Format base coefficient (derivative) name
        basename = "{name}{count}".format(name=names.w, count=mt.terminal.count())
        access = format_mt_name(basename, mt)
        return access

    def quadrature_weight(self, e, mt, tabledata): # TODO: let num_points be arg here
        # TODO: Need num_points to identify weights array name?
        #        Or maybe place it in tabledata:
        #name, dummy, num_points = tabledata
        num_points = self._current_num_points # TODO: Fix this hack!
        access = ArrayAccess(self.weights_array_name(num_points), self.quadrature_loop_index())
        return access

    def spatial_coordinate(self, e, mt, tabledata):
        ffc_assert(not mt.global_derivatives, "Not expecting derivatives of SpatialCoordinates.")
        ffc_assert(not mt.local_derivatives, "Not expecting derivatives of SpatialCoordinates.")
        #ffc_assert(not mt.restriction, "Not expecting restriction of SpatialCoordinates.")
        ffc_assert(not mt.averaged, "Not expecting average of SpatialCoordinates.")

        if self.physical_coordinates_known:
            # In a context where the physical coordinates are available in existing variables.
            # TODO: Need num_points to identify points array name?
            #        Or maybe place it in tabledata:
            #name, dummy, num_points = tabledata
            num_points = self._current_num_points # TODO: Fix this hack!
            access = ArrayAccess(self.points_array_name(num_points),
                                 (self.quadrature_loop_index(), mt.flat_component))

        elif mt.terminal.domain().coordinates() is not None:
            # No special variable should exist in this case.
            error("Expecting spatial coordinate to be symbolically rewritten.")

        else:
            # In a context where physical coordinates are computed by code generated by us.
            access = format_mt_name(names.x, mt)

        return access

    def cell_coordinate(self, e, mt, tabledata):
        ffc_assert(not mt.global_derivatives, "Not expecting derivatives of CellCoordinates.")
        ffc_assert(not mt.local_derivatives, "Not expecting derivatives of CellCoordinates.")
        ffc_assert(not mt.averaged, "Not expecting average of CellCoordinates.")

        assert not mt.restriction # FIXME: Not used!

        if self.physical_coordinates_known:
            # No special variable should exist in this case.
            error("Expecting reference coordinate to be symbolically rewritten.")

        else:
            # TODO: Need num_points to identify points array name?
            #        Or maybe place it in tabledata:
            #name, dummy, num_points = tabledata
            num_points = self._current_num_points # TODO: Fix this hack!
            access = ArrayAccess(self.points_array_name(num_points),
                                 (self.quadrature_loop_index(), mt.flat_component))

        return access

    def jacobian(self, e, mt, tabledata):
        ffc_assert(not mt.global_derivatives, "Not expecting derivatives of Jacobian.")
        ffc_assert(not mt.local_derivatives, "Not expecting derivatives of Jacobian.")
        ffc_assert(not mt.averaged, "Not expecting average of Jacobian.")

        access = format_mt_name(names.J, mt)
        return access

    def cell_facet_jacobian(self, e, mt, tabledata):
        cellname = mt.terminal.domain().cell().cellname()
        if cellname in ("triangle", "tetrahedron"):
            tablename = "{0}_reference_facet_jacobian".format(cellname)
            facet = format_entity_name("facet", mt.restriction)
            access = ArrayAccess(tablename, (facet, mt.component[0], mt.component[1]))
        elif cellname == "interval":
            error("The reference facet jacobian doesn't make sense for interval cell.")
        else:
            error("Unhandled cell types {0}.".format(cellname))
        return access

    def cell_orientation(self, e, mt, tabledata):
        # TODO: error if not in manifold case
        access = "cell_orientation"
        return access

    def facet_orientation(self, e, mt, tabledata):
        cellname = mt.terminal.domain().cell().cellname()
        if cellname in ("triangle", "tetrahedron"):
            tablename = "{0}_facet_orientations".format(cellname)
            facet = format_entity_name("facet", mt.restriction)
            access = ArrayAccess(tablename, (facet,))
        elif cellname == "interval":
            error("The facet orientation doesn't make sense for interval cell.")
        else:
            error("Unhandled cell types {0}.".format(cellname))
        return access

    def facet_normal(self, e, mt, tabledata):
        error("Expecting {0} to be replaced with lower level types in symbolic preprocessing.".format(type(e)))

    def cell_normal(self, e, mt, tabledata):
        error("Expecting {0} to be replaced with lower level types in symbolic preprocessing.".format(type(e)))

    def jacobian_inverse(self, e, mt, tabledata):
        error("Expecting {0} to be replaced with lower level types in symbolic preprocessing.".format(type(e)))

    def jacobian_determinant(self, e, mt, tabledata):
        error("Expecting {0} to be replaced with lower level types in symbolic preprocessing.".format(type(e)))

    def facet_jacobian(self, e, mt, tabledata):
        error("Expecting {0} to be replaced with lower level types in symbolic preprocessing.".format(type(e)))

    def facet_jacobian_inverse(self, e, mt, tabledata):
        error("Expecting {0} to be replaced with lower level types in symbolic preprocessing.".format(type(e)))

    def facet_jacobian_determinant(self, e, mt, tabledata):
        error("Expecting {0} to be replaced with lower level types in symbolic preprocessing.".format(type(e)))
