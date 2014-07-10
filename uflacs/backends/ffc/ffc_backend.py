"""
FIXME:
- Fix elementwise bugs in graph_rebuild
- Fix FacetNormal (sign?)
- Restrict Jacobian etc. in ufl change_to_reference
- Various restriction issues, check all quantities

DONE - Inject CG1 tables for x and J at some point and pass them here through tabledata
DONE - Check the code for J carefully
HACKED- Pass num_points to weight, x, X
DONE - Handle cell restriction in generate_domain_dof_access
DONE - Make reference_facet_jacobi tables and use below
"""

from six import iterkeys
from six.moves import xrange as range
import ufl
from ufl.permutation import build_component_numbering
from ufl.classes import Argument, Coefficient, GeometricQuantity
from ufl.algorithms import MultiFunction

from ffc.log import ffc_assert, warning, error

from uflacs.codeutils.cpp_statement_formatting_rules import CppStatementFormattingRules
langfmt = CppStatementFormattingRules()

from uflacs.codeutils.format_code import (format_code,
                                          ForRange,
                                          VariableDecl,
                                          ArrayAccess,
                                          AssignAdd, Assign,
                                          Add, Sub, Mul,
                                          Sum,)

class Names: # TODO: This is not used much anymore, integrate in backend class
    def __init__(self):
        # Topology argument names
        self.vertex = "vertex"
        self.facet = "facet"

        # Geometry names
        self.vertex_coordinates = "vertex_coordinates"
        self.xi = "xi"
        self.x = "x"
        self.J = "J"
        self.K = "K"
        self.detJ = "detJ"
        self.det = "det"

        # Quadrature rule
        self.points = "points"
        self.weights = "weights"

        # Quadrature temps
        self.qw = "qw"
        self.D = "D"

        # (Base)name for intermediate registers
        self.s = "s"

        # Element tensor
        self.A = "A"

        # Coefficient dofs array
        self.w = "w"

        # Basenames for function components
        self.wbase = "w"
        self.vbase = "v"
        self.dwbase = "dw"
        self.dvbase = "dv"

        # Loop indices
        self.iq = "iq"   # Quadrature loop
        self.ic = "ic"   # Coefficient accumulation loop
        self.ia = "ia"   # Argument dof loop
        self.ild = "ild" # Local derivative accumulation loop

        # Rules, make functions?
        self.restriction_postfix = { "+": "_0", "-": "_1", None: "" } # TODO: Use this wherever we need it?

names = Names()


def format_entity_name(entitytype, r):
    if entitytype == "cell":
        entity = "0" #None # TODO: Keep 3D tables and use entity 0 for cells or make tables 2D and use None?
    elif entitytype == "facet":
        entity = names.facet + names.restriction_postfix[r]
    elif entitytype == "vertex":
        entity = names.vertex
    return entity

def format_mt_der(mt):
    # Expecting only local derivatives here
    assert not mt.global_derivatives
    # Add derivatives to name
    if mt.local_derivatives:
        der = "_d{0}".format(''.join(map(str, mt.local_derivatives)))
    else:
        der = ""
    return der

def format_mt_comp(mt):
    # Add flattened component to name (TODO: this should be the local component?)
    if mt.component:
        comp = "_c{0}".format(mt.flat_component)
    else:
        comp = ""
    return comp

def format_mt_avg(mt):
    # Add averaged state to name
    if mt.averaged:
        avg = "_a{0}".format(mt.averaged)
    else:
        avg = ""
    return avg

def format_mt_res(mt):
    return names.restriction_postfix[mt.restriction].replace("_", "_r")

def format_mt_name(basename, mt):
    access = "{basename}{avg}{res}{der}{comp}".format(basename=basename,
                                                      avg=format_mt_avg(mt),
                                                      res=format_mt_res(mt),
                                                      der=format_mt_der(mt),
                                                      comp=format_mt_comp(mt))
    return access

def ufc_restriction_postfix(restriction):
    # TODO: Get restriction postfix from somewhere central
    if restriction == "+":
        res = "_0"
    elif restriction == "-":
        res = "_1"
    else:
        res = ""
    return res

#from uflacs.backends.ffc.ffc_statement_formatter import format_element_table_access
#from ufl.utils.derivativetuples import derivative_listing_to_counts
#def generate_element_table_access(mt):
#    # FIXME: See  format_element_table_access  get_element_table_data
#    #entity = format_entity_name(self.ir["entitytype"], mt.restriction)
#    #return ArrayAccess(uname, (entity, names.iq, dof_number))
#    return "FE[0]" # FIXME

#def generate_geometry_table_access(mt):
#    return "FJ[0]" # FIXME

def generate_coefficient_dof_access(coefficient, dof_number):
    # TODO: Add domain_number = self.ir["domain_numbering"][coefficient.domain().domain_key()]
    # TODO: Flatten dofs array and use CRS lookup table.
    # TODO: Apply integral specific renumbering.
    return ArrayAccess(names.w, (coefficient.count(), dof_number))

def generate_domain_dof_access(num_vertices, gdim, vertex, component, restriction):
    # TODO: Add domain number as argument here, and {domain_offset} to array indexing:
    #domain_offset = self.ir["domain_offsets"][domain_number]
    vc = names.vertex_coordinates + names.restriction_postfix[restriction]
    return ArrayAccess(vc, Add(Mul(gdim, vertex), component))

def generate_domain_dofs_access(num_vertices, gdim, restriction):
    # TODO: Add domain number as argument here, and {domain_offset} to array indexing:
    # FIXME: Handle restriction here
    #domain_offset = self.ir["domain_offsets"][domain_number]
    return [generate_domain_dof_access(num_vertices, gdim, vertex, component, restriction)
            for component in range(gdim)
            for vertex in range(num_vertices)]

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
        ffc_assert(not mt.restriction, "Not expecting restriction of SpatialCoordinates.")
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



class FFCDefinitionsBackend(MultiFunction):
    """FFC specific cpp formatter class."""
    def __init__(self, ir, parameters):
        MultiFunction.__init__(self)

        # Store ir and parameters
        self.ir = ir
        self.parameters = parameters

        # Configure definitions behaviour
        self.physical_coordinates_known = self.ir["integral_type"] == "quadrature"

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

            gdim = mt.terminal.domain().cell().geometric_dimension()
            num_vertices = mt.terminal.domain().cell().topological_dimension() + 1 # TODO: Get from cellname?

            uname, begin, end = tabledata

            # access here is e.g. x0, component 0 of x

            ffc_assert(0 <= begin <= end <= num_vertices*gdim,
                          "Assuming linear element for affine simplices here.")
            entity = format_entity_name(self.ir["entitytype"], mt.restriction)
            vertex = names.ic
            iq = names.iq

            if 0: # Generated loop version:
                table_access = ArrayAccess(uname, (entity, iq, vertex))
                dof_access = generate_domain_dof_access(num_vertices, gdim, vertex,
                                                        mt.flat_component, mt.restriction)
                prod = Mul(dof_access, table_access)

                # Loop to accumulate linear combination of dofs and tables
                code += [VariableDecl("double", access, "0.0")]
                code += [ForRange(vertex, begin, end, body=[AssignAdd(access, prod)])]

            else: # Inlined version:
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

    def jacobian(self, e, mt, tabledata, access): # TODO: Handle jacobian and restricted in change_to_reference
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

            gdim = mt.terminal.domain().cell().geometric_dimension()
            num_vertices = mt.terminal.domain().cell().topological_dimension() + 1 # TODO: Get from cellname?

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
        # Currently the table is inserted in self.initial()
        # TODO: Define table in ufc_geometry? Or insert among regular tables?
        code = []
        return code

    def cell_orientation(self, e, mt, tabledata, access):
        code = []
        return code

    def facet_orientation(self, e, mt, tabledata, access):
        code = []
        return code

    def facet_normal(self, e, mt, tabledata, access):
        code = []
        return code

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
