
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
