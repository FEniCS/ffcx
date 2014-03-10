"""
FIXME:
- Make reference_facet_jacobi tables and use below
- Fix FacetNormal (sign?)
- Restrict Jacobian etc. in ufl change_to_reference

DONE - Inject CG1 tables for x and J at some point and pass them here through tabledata
DONE - Check the code for J carefully
HACKED- Pass num_points to weight, x, X
DONE - Handle cell restriction in generate_domain_dof_access
"""

from ufl.common import component_to_index
from ufl.permutation import build_component_numbering
from ufl.classes import Argument, Coefficient, GeometricQuantity
from ufl.algorithms import MultiFunction

from uflacs.utils.log import uflacs_assert, warning, error

from uflacs.geometry.default_names import names
from uflacs.codeutils.cpp_statement_formatting_rules import CppStatementFormattingRules
langfmt = CppStatementFormattingRules()

from uflacs.codeutils.format_code import (format_code,
                                          ForRange,
                                          VariableDecl,
                                          ArrayAccess,
                                          AssignAdd, Assign,
                                          Add, Sub, Mul,
                                          Sum,)


# TODO: The organization of code utilities is a bit messy...
#from uflacs.backends.ffc.ffc_statement_formatter import format_element_table_access
#from uflacs.elementtables.table_utils import derivative_listing_to_counts
from uflacs.elementtables.table_utils import flatten_component


# TODO: Insert in ufc_geometry or ufc_cells header
ufc_reference_facet_jacobian_tables = """
/// --- Jacobians of reference facet cell to reference cell coordinate mappings by UFC conventions ---

static const double triangle_reference_facet_jacobian[3][2][1] = {
  { {-1.0}, { 1.0} },
  { { 0.0}, { 1.0} },
  { { 1.0}, { 0.0} },
  };

static const double tetrahedron_reference_facet_jacobian[4][3][2] = {
  { {-1.0, -1.0}, {1.0, 0.0}, {0.0, 1.0} },
  { { 0.0,  0.0}, {1.0, 0.0}, {0.0, 1.0} },
  { { 1.0,  0.0}, {0.0, 0.0}, {0.0, 1.0} },
  { { 1.0,  0.0}, {0.0, 1.0}, {0.0, 0.0} },
  };
"""


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
        self.ir = ir
        self.parameters = parameters

        self.physical_coordinates_known = self.ir["domain_type"] == "quadrature"

        self._current_num_points, = self.ir["uflacs"]["expr_ir"].keys() # TODO: Hack!

    def precision_float(self, f):
        return "%g" % f # TODO: Control float formatting precision here

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
        iq = names.iq
        idof = "{0}{1}".format(names.ia, mt.terminal.number())
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
        uflacs_assert(mt.flat_component == vi2si[mt.component], "Incompatible component flattening!")

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

    def quadrature_weight(self, e, mt, tabledata):
        # TODO: Need num_points to identify weights array name?
        #        Or maybe place it in tabledata:
        #name, dummy, num_points = tabledata
        num_points = self._current_num_points # TODO: Fix this hack!
        weights = "{0}{1}".format(names.weights, num_points)
        access = ArrayAccess(weights, names.iq)
        return access

    def spatial_coordinate(self, e, mt, tabledata):
        uflacs_assert(not mt.global_derivatives, "")
        uflacs_assert(not mt.local_derivatives, "")
        uflacs_assert(not mt.restriction, "")
        uflacs_assert(not mt.averaged, "")

        if self.physical_coordinates_known:
            # TODO: Need num_points to identify points array name?
            #        Or maybe place it in tabledata:
            #name, dummy, num_points = tabledata
            num_points = self._current_num_points # TODO: Fix this hack!
            points = "{0}{1}".format(names.points, num_points)
            access = ArrayAccess(points, (names.iq, mt.flat_component))
        elif mt.terminal.domain().coordinates() is not None:
            error("Expecting spatial coordinate to be symbolically rewritten.")
        else:
            access = "{basename}{comp}".format(basename=names.x,
                                               comp=format_mt_comp(mt))
        return access

    def reference_coordinate(self, e, mt, tabledata):
        uflacs_assert(not mt.global_derivatives, "")
        uflacs_assert(not mt.local_derivatives, "")
        uflacs_assert(not mt.restriction, "")
        uflacs_assert(not mt.averaged, "")

        if self.physical_coordinates_known:
            error("Expecting reference coordinate to be symbolically rewritten.")
        else:
            # TODO: Need num_points to identify points array name?
            #        Or maybe place it in tabledata:
            #name, dummy, num_points = tabledata
            num_points = self._current_num_points # TODO: Fix this hack!
            points = "{0}{1}".format(names.points, num_points)
            access = ArrayAccess(points, (names.iq, mt.flat_component))

        return access

    def jacobian(self, e, mt, tabledata):
        uflacs_assert(not mt.global_derivatives, "")
        uflacs_assert(not mt.local_derivatives, "")
        uflacs_assert(not mt.restriction, "")
        uflacs_assert(not mt.averaged, "")

        access = "{basename}{comp}".format(basename=names.J,
                                           comp=format_mt_comp(mt))
        return access

    def reference_facet_jacobian(self, e, mt, tabledata):
        cellname = mt.terminal.domain().cell().cellname()
        if cellname in ("triangle", "tetrahedron"):
            tablename = "{0}_reference_facet_jacobian".format(cellname)
            facet = names.facet + names.restriction_postfix[mt.restriction]
            access = ArrayAccess(tablename, (facet, mt.component[0], mt.component[1]))
        elif cellname == "interval":
            error("The reference facet jacobian doesn't make sense for interval cell.")
        else:
            error("Unhandled cell types {0}.".format(cellname))
        return access

    def facet_normal(self, e, mt, tabledata): # FIXME
        access = "nFIXME" # FIXME: Computed from tables?
        return access

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
        self.ir = ir
        self.parameters = parameters

        self.physical_coordinates_known = self.ir["domain_type"] == "quadrature"

    def initial(self): # TODO: Need something like this?
        code = []
        code += [ufc_reference_facet_jacobian_tables]
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
            uflacs_assert(mt.terminal.domain().coordinates() is None,
                          "Assuming coefficient field symbolically inserted before this point.")
            # Reference coordinates are known, no coordinate field, so we compute
            # this component as linear combination of vertex_coordinates "dofs" and table

            gdim = mt.terminal.domain().cell().geometric_dimension()
            num_vertices = mt.terminal.domain().cell().topological_dimension() + 1 # TODO: Get from cellname?

            uname, begin, end = tabledata

            # access here is e.g. x0, component 0 of x

            uflacs_assert(0 <= begin <= end <= num_vertices*gdim,
                          "Assuming linear element for affine simplices here.")
            entity = format_entity_name(self.ir["entitytype"], mt.restriction)
            vertex = names.ic

            if 0: # Generated loop version:
                table_access = ArrayAccess(uname, (entity, names.iq, vertex))
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
                    table_access = ArrayAccess(uname, (entity, names.iq, Sub(idof, begin)))
                    prods += [Mul(dof_access[idof], table_access)]

                # Inlined loop to accumulate linear combination of dofs and tables
                code += [VariableDecl("const double", access, Sum(prods))]

        return code

    def reference_coordinate(self, e, mt, tabledata, access):
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
            uflacs_assert(mt.terminal.domain().coordinates() is None,
                          "Assuming coefficient field symbolically inserted before this point.")
            # Reference coordinates are known, no coordinate field, so we compute
            # this component as linear combination of vertex_coordinates "dofs" and table

            gdim = mt.terminal.domain().cell().geometric_dimension()
            num_vertices = mt.terminal.domain().cell().topological_dimension() + 1 # TODO: Get from cellname?

            uname, begin, end = tabledata

            # access here is e.g. J_0, component 0 of J

            uflacs_assert(0 <= (end-begin) <= num_vertices,
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

    def reference_facet_jacobian(self, e, mt, tabledata, access):
        # Currently the table is inserted in self.initial()
        # TODO: Define table in ufc_geometry? Or insert among regular tables?
        code = []
        return code

    def facet_normal(self, e, mt, tabledata, access): # FIXME
        "TODO: Insert symbolically?"
        code = []
        code += [VariableDecl("double", access, "1.0 /* FIXME */")] # FIXME
        return code

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
