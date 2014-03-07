
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
                                          Add, Sub, Mul,)


# TODO: The organization of code utilities is a bit messy...
#from uflacs.backends.ffc.ffc_statement_formatter import format_element_table_access
#from uflacs.elementtables.table_utils import derivative_listing_to_counts
from uflacs.elementtables.table_utils import flatten_component



def format_entity_name(entitytype, r):
    if entitytype == "cell":
        entity = "0" #None # FIXME: Keep 3D tables and use entity 0 for cells or make tables 2D and use None?
    elif entitytype == "facet":
        entity = names.facet + names.restriction_postfix[r]
    elif entitytype == "vertex":
        entity = names.vertex
    return entity



def generate_element_table_access(mt):
    # FIXME: See  format_element_table_access  get_element_table_data
    return "FE[0]" # FIXME

def generate_geometry_table_access(mt):
    return "FJ[0]" # FIXME

def generate_coefficient_dof_access(coefficient, dof_number):
    # TODO: Add domain_number = self.ir["domain_numbering"][coefficient.domain().domain_key()]
    # TODO: Flatten dofs array and use CRS lookup table.
    # TODO: Apply integral specific renumbering.
    return ArrayAccess(names.w, (coefficient.count(), dof_number))

def generate_domain_dof_access(num_vertices, gdim, vertex, component):
    # TODO: Add domain number as argument here, and {domain_offset} to array indexing:
    #domain_offset = self.ir["domain_offsets"][domain_number]
    name = names.vertex_coordinates
    return "{name}[{gdim}*{vertex} + {component}]".format(**locals())


class FFCAccessBackend(MultiFunction):
    """FFC specific cpp formatter class."""
    def __init__(self, ir, parameters):
        MultiFunction.__init__(self)
        self.ir = ir
        self.parameters = parameters

    def precision_float(self, f):
        return "%g" % f # TODO: Control float formatting precision here

    # === Multifunction handlers for all modified terminal types, basic C++ types are covered by base class ===

    def expr(self, e, mt, tabledata):
        error("Missing handler for type {0}.".format(str(e)))

    def argument(self, e, mt, tabledata):

        # Expecting only local derivatives and values here
        assert not (mt.global_derivatives and any(mt.global_derivatives)) # FIXME: Clarify format
        #assert mt.global_component is None

        # No need to store basis function value in its own variable, just get table value directly
        uname, begin, end = tabledata
        entity = format_entity_name(self.ir["entitytype"], mt.restriction) # FIXME: Clarify mt.restriction format
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
        flat_component = vi2si[mt.component]

        # Offset index if on second cell in interior facet integral
        if mt.restriction == "-":
            flat_component += len(si2vi)

        # Return direct reference to dof array
        return generate_coefficient_dof_access(mt.terminal, flat_component)

    def _varying_coefficient(self, e, mt, tabledata):
        t = mt.terminal

        # Expecting only local derivatives and values here
        assert not (mt.global_derivatives and any(mt.global_derivatives))
        #assert mt.global_component is None

        # Add derivatives to name
        if mt.local_derivatives and any(mt.local_derivatives):
            derivative_str = ''.join(map(str, mt.local_derivatives))
            der = "_d{0}".format(derivative_str)
        else:
            der = ""

        # Add flattened component to name (TODO: this should be the local component)
        if t.shape():
            flat_component = flatten_component(mt.component, t.shape(), t.element().symmetry())
            comp = "_c{0}".format(flat_component)
        else:
            comp = ""

        # Add restriction to name
        res = names.restriction_postfix[mt.restriction]

        # Format base coefficient (derivative) name
        access = "{name}{count}{der}{comp}{res}".format(name=names.w, count=t.count(), der=der, comp=comp, res=res)
        return access


    def _old_geometric_quantity(self, o, mt, tabledata):
        "Generic rendering of variable names for all piecewise constant geometric quantities."
        uflacs_assert(not mt.global_derivatives and not mt.local_derivatives,
                      "Compiler should be able to simplify derivatives of geometry.")

        # Simply using the UFL str to define the name in the generated code, ensures consistency
        name = str(o)
        if mt.restriction:
            res = names.restriction_postfix[mt.restriction]
            name = name + res

        # Indexing if there is a shape
        sh = o.shape()
        if sh:
            code = langfmt.array_access(name, component_to_index(mt.component, sh))
        else:
            code = name

        return code

    def spatial_coordinate(self, e, mt, tabledata):
        # FIXME: See define_coord_vars
        access = "x" # FIXME
        return access

    def local_coordinate(self, e, mt, tabledata):
        # FIXME: See define_coord_vars
        access = "X" # FIXME
        return access

    def jacobian(self, e, mt, tabledata):
        # FIXME: See _define_piecewise_geometry  define_piecewise_geometry
        access = "J" # FIXME
        return access

    def jacobian_inverse(self, e, mt, tabledata):
        access = "K" # FIXME
        return access

    def jacobian_determinant(self, e, mt, tabledata):
        access = "detJ" # FIXME
        return access

    def facet_jacobian(self, e, mt, tabledata):
        access = "J" # FIXME
        return access

    def facet_jacobian_inverse(self, e, mt, tabledata):
        access = "K" # FIXME
        return access

    def facet_jacobian_determinant(self, e, mt, tabledata):
        access = "detJ" # FIXME
        return access

    def facet_normal(self, e, mt, tabledata):
        access = "n" # FIXME
        return access



class FFCDefinitionsBackend(MultiFunction):
    """FFC specific cpp formatter class."""
    def __init__(self, ir, parameters):
        MultiFunction.__init__(self)
        self.ir = ir
        self.parameters = parameters

    def expr(self, t, mt, name, tabledata, access):
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
            entity = format_entity_name(self.ir["entitytype"], mt.restriction) # FIXME: Clarify mt.restriction format
            iq = names.iq
            idof = names.ic

            dof = format_code(Sub(idof, begin))
            table_access = ArrayAccess(uname, (entity, iq, dof))

            dof_access = generate_coefficient_dof_access(mt.terminal, idof)

            prod = Mul(dof_access, table_access)
            body = [AssignAdd(access, prod)]

            # Loop to accumulate linear combination of dofs and tables
            code += [VariableDecl("const double", access, "0.0")]
            code += [ForRange(idof, begin, end, body=body)]

        return code

    # FIXME: See _define_piecewise_geometry  define_piecewise_geometry

    def spatial_coordinate(self, e, mt, tabledata, access):
        # FIXME: See define_coord_vars
        code += [VariableDecl("double", access, "1.0")] # FIXME
        return code

    def local_coordinate(self, e, mt, tabledata, access):
        # FIXME: See define_coord_vars
        code += [VariableDecl("double", access, "1.0")] # FIXME
        return code

    def jacobian(self, e, mt, tabledata, access):
        code += [VariableDecl("double", access, "1.0")] # FIXME
        return code

    def jacobian_inverse(self, e, mt, tabledata, access):
        code += [VariableDecl("double", access, "1.0")] # FIXME
        return code

    def jacobian_determinant(self, e, mt, tabledata, access):
        code += [VariableDecl("double", access, "1.0")] # FIXME
        return code

    def facet_jacobian(self, e, mt, tabledata, access):
        code += [VariableDecl("double", access, "1.0")] # FIXME
        return code

    def facet_jacobian_inverse(self, e, mt, tabledata, access):
        code += [VariableDecl("double", access, "1.0")] # FIXME
        return code

    def facet_jacobian_determinant(self, e, mt, tabledata, access):
        code += [VariableDecl("double", access, "1.0")] # FIXME
        return code

    def facet_normal(self, e, mt, tabledata, access):
        code += [VariableDecl("double", access, "1.0")] # FIXME
        return code
