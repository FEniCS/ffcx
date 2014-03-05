
from ufl.common import component_to_index
from ufl.permutation import build_component_numbering
from ufl.algorithms import MultiFunction

from uflacs.utils.log import uflacs_assert, warning, error

# TODO: The organization of code utilities is a bit messy...
from uflacs.codeutils.cpp_expr_formatting_rules import CppFormattingRules
from uflacs.geometry.default_names import names
from uflacs.backends.ffc.ffc_statement_formatter import langfmt
from uflacs.backends.ffc.ffc_statement_formatter import (format_element_table_access, format_entity_name)
from uflacs.elementtables.table_utils import derivative_listing_to_counts, flatten_component


class FFCLanguageFormatter(MultiFunction, CppFormattingRules):
    """FFC specific cpp formatter class."""
    def __init__(self, ir):
        MultiFunction.__init__(self)
        CppFormattingRules.__init__(self)

    def modified_terminal(self, o, mt):
        return str(o) # FIXME

    # Core terminals
    coefficient = modified_terminal
    argument = modified_terminal
    geometric_quantity = modified_terminal

    # Modifiers
    grad = modified_terminal
    local_grad = modified_terminal
    restriction = modified_terminal
    indexed = modified_terminal
    facet_avg = modified_terminal
    cell_avg = modified_terminal

class Backend(object):
    def generate_modified_terminal_definition(self, mt):
        if isinstance(mt.terminal, Argument):
            return self.generate_argument_definition(mt)
        if isinstance(mt.terminal, Coefficient):
            return self.generate_coefficient_definition(mt)
        if isinstance(mt.terminal, GeometricQuantity):
            return self.generate_geometric_quantity_definition(mt)
        error("Unknown terminal type {0}.".format(type(mt.terminal)))

    def generate_modified_terminal_access(self, mt):
        if isinstance(mt.terminal, Argument):
            return self.generate_argument_access(mt)
        if isinstance(mt.terminal, Coefficient):
            return self.generate_coefficient_access(mt)
        if isinstance(mt.terminal, GeometricQuantity):
            return self.generate_geometric_quantity_access(mt)
        error("Unknown terminal type {0}.".format(type(mt.terminal)))


def format_entity_name(entitytype, r):
    if entitytype == "cell":
        entity = "0" #None # FIXME: Keep 3D tables and use entity 0 for cells or make tables 2D and use None?
    elif entitytype == "facet":
        entity = names.facet + names.restriction_postfix[r]
    elif entitytype == "vertex":
        entity = names.vertex
    return entity


class FFCBackend(Backend):

    def generate_element_table_access(self, mt):
        return "FE[0]"

    def generate_geometry_table_access(self, mt):
        return "FJ[0]"

    def generate_coefficient_dof_access(self, coefficient_number, dof_number):
        # TODO: Add domain number
        return langfmt.array_access(names.w, coefficient_number, dof_number)

    def generate_domain_dof_access(self, num_vertices, gdim, vertex, component):
        # TODO: Add domain number as argument here, and {domain_offset} to array indexing:
        #domain_offset = self.ir["domain_offsets"][domain_number]
        name = names.vertex_coordinates
        return "{name}[{gdim}*{vertex} + {component}]".format(**locals())


    def generate_argument_definition(self, mt):
        code = []
        return code

    def generate_argument_access(self, mt):
        # Expecting only local derivatives and values here
        assert not (mt.global_derivatives and any(mt.global_derivatives))
        #assert mt.global_component is None

        t = mt.terminal
        element = t.element()

        # Pick entity index variable name, following ufc argument names
        entity = format_entity_name(self._entitytype, mt.restriction) # FIXME

        idof = "{0}{1}".format(names.ia, t.number())

        # TODO: Local component
        flat_component = flatten_component(mt.component, element.value_shape(), element.symmetry())

        # No need to store basis function value in its own variable, just get table value directly
        access = format_element_table_access(self._ir, self._entitytype, self._num_points,
                                             element, flat_component, (), entity, idof, True)

        return access


    def generate_coefficient_access(self, mt):
        t = mt.terminal
        if t.is_cellwise_constant():
            access = self.generate_constant_coefficient_access(mt)
        else:
            access = self.generate_varying_coefficient_access(mt)
        return access

    def generate_constant_coefficient_access(self, mt):
        # Map component to flat index
        vi2si, si2vi = build_component_numbering(mt.terminal.shape(), mt.terminal.element().symmetry())
        flat_component = vi2si[mt.component]

        # Offset index if on second cell in interior facet integral
        if mt.restriction == "-":
            flat_component += len(si2vi)

        # Return direct reference to dof array
        return self.generate_coefficient_dof_access(mt.terminal.count(), flat_component)

    def generate_varying_coefficient_access(self, mt):
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
        access = "%s%d%s%s%s" % (names.w, t.number(), der, comp, res)
        return access

    def generate_coefficient_definition(self, mt):
        code = []
        if not t.is_cellwise_constant():
            name = generate_varying_coefficient_access(mt)
            value = "1.0" # FIXME: Linear combination of dofs and tables
            code += ["const double {name} = {value};".format(name, value)]
        return code


    def generate_geometric_quantity_definition(self, mt):
        code = []
        code += ["double geometry = 1.0;"] # FIXME
        return code

    def generate_geometric_quantity_access(self, mt):
        access = "geometry" # FIXME
        return access


def generate_piecewise_partition():

    for o in modified_terminals:
        mt = analyse_modified_terminals2(o)

        definition = backend.generate_modified_terminal_definition(mt)
        access = backend.generate_modified_terminal_access(mt)

        code += definition
        langfmt.cache[o] = access
