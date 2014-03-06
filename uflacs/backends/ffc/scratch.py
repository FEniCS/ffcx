
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


def format_entity_name(entitytype, r):
    if entitytype == "cell":
        entity = "0" #None # FIXME: Keep 3D tables and use entity 0 for cells or make tables 2D and use None?
    elif entitytype == "facet":
        entity = names.facet + names.restriction_postfix[r]
    elif entitytype == "vertex":
        entity = names.vertex
    return entity


class FFCLanguageFormatter(MultiFunction, CppFormattingRules):
    """FFC specific cpp formatter class."""
    def __init__(self, ir):
        MultiFunction.__init__(self)
        CppFormattingRules.__init__(self)

    def generate_element_table_access(self, mt):
        # FIXME: See  format_element_table_access  get_element_table_data
        return "FE[0]" # FIXME

    def generate_geometry_table_access(self, mt):
        return "FJ[0]" # FIXME

    def generate_coefficient_dof_access(self, coefficient_number, dof_number):
        # TODO: Add domain number
        return langfmt.array_access(names.w, coefficient_number, dof_number)

    def generate_domain_dof_access(self, num_vertices, gdim, vertex, component):
        # TODO: Add domain number as argument here, and {domain_offset} to array indexing:
        #domain_offset = self.ir["domain_offsets"][domain_number]
        name = names.vertex_coordinates
        return "{name}[{gdim}*{vertex} + {component}]".format(**locals())

    # === Multifunction handlers for all modified terminal types, basic C++ types are covered by base class ===

    def modified_terminal(self, o):
        "Analyse underlying modified terminal structure and pass on to handler of terminal type."
        mt = analyse_modified_terminal2(o)
        return self(mt.terminal, mt)

    # Terminal modifiers TODO: Any more? Add local_value/pullback
    grad = modified_terminal
    local_grad = modified_terminal
    restriction = modified_terminal
    indexed = modified_terminal
    facet_avg = modified_terminal
    cell_avg = modified_terminal

    def argument(self, e, mt):
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
                                             element, flat_component, (), entity, idof, True) # FIXME

        return access

    def coefficient(self, e, mt):
        t = mt.terminal
        if t.is_cellwise_constant():
            access = self._constant_coefficient(e, mt)
        else:
            access = self._varying_coefficient(e, mt)
        return access

    def _constant_coefficient(self, e, mt):
        # Map component to flat index
        vi2si, si2vi = build_component_numbering(mt.terminal.shape(), mt.terminal.element().symmetry())
        flat_component = vi2si[mt.component]

        # Offset index if on second cell in interior facet integral
        if mt.restriction == "-":
            flat_component += len(si2vi)

        # Return direct reference to dof array
        return self.generate_coefficient_dof_access(mt.terminal.count(), flat_component)

    def _varying_coefficient(self, e, mt):
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


    def _old_geometric_quantity(self, o, mt):
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

    def spatial_coordinate(self, e, mt):
        # FIXME: See define_coord_vars
        access = "x" # FIXME
        return access

    def local_coordinate(self, e, mt):
        # FIXME: See define_coord_vars
        access = "X" # FIXME
        return access

    def jacobian(self, e, mt):
        # FIXME: See _define_piecewise_geometry  define_piecewise_geometry
        access = "J" # FIXME
        return access

    def jacobian_inverse(self, e, mt):
        access = "K" # FIXME
        return access

    def jacobian_determinant(self, e, mt):
        access = "detJ" # FIXME
        return access

    def facet_jacobian(self, e, mt):
        access = "J" # FIXME
        return access

    def facet_jacobian_inverse(self, e, mt):
        access = "K" # FIXME
        return access

    def facet_jacobian_determinant(self, e, mt):
        access = "detJ" # FIXME
        return access

    def facet_normal(self, e, mt):
        access = "n" # FIXME
        return access

    # === Generate code definitions ===

    def generate_modified_terminal_definition(self, mt):
        if isinstance(mt.terminal, Argument):
            return self.generate_argument_definition(mt)
        if isinstance(mt.terminal, Coefficient):
            return self.generate_coefficient_definition(mt)
        if isinstance(mt.terminal, GeometricQuantity):
            return self.generate_geometric_quantity_definition(mt)
        error("Unknown terminal type {0}.".format(type(mt.terminal)))

    def generate_argument_definition(self, mt):
        code = []
        return code

    def generate_coefficient_definition(self, mt):
        # FIXME: See ffc_statement_formatter
        #  define_coord_dependent_coefficients
        #  _define_coord_dependent_coefficient
        #
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
