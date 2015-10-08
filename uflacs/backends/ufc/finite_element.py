
from ufl import product
from uflacs.backends.ufc.generator import ufc_generator
from uflacs.backends.ufc.templates import finite_element_header, finite_element_implementation
from uflacs.backends.ufc.utils import generate_return_new_switch

class ufc_finite_element(ufc_generator):
    def __init__(self):
        ufc_generator.__init__(self, finite_element_header, finite_element_implementation)

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

    def evaluate_basis(self, L, ir): # FIXME: port this
        return "FIXME" + ir["evaluate_basis"]

    def evaluate_basis_derivatives(self, L, ir): # FIXME: port this
        return "FIXME" + ir["evaluate_basis_derivatives"]

    def evaluate_basis_all(self, L, ir): # FIXME: port this
        return "FIXME" + ir["evaluate_basis_all"]

    def evaluate_basis_derivatives_all(self, L, ir): # FIXME: port this
        return "FIXME" + ir["evaluate_basis_derivatives_all"]

    def evaluate_dof(self, L, ir): # FIXME: port this
        return "FIXME" + ir["evaluate_dof"]

    def evaluate_dofs(self, L, ir): # FIXME: port this
        return "FIXME" + ir["evaluate_dofs"]

    def interpolate_vertex_values(self, L, ir): # FIXME: port this
        return "FIXME" + ir["interpolate_vertex_values"]

    def tabulate_dof_coordinates(self, L, ir): # FIXME: port this
        coords = ir["tabulate_dof_coordinates"]
        code = "FIXME" + str(coords)
        return code

    def num_sub_elements(self, L, ir):
        n = ir["num_sub_elements"]
        return L.Return(L.LiteralInt(n))

    def create_sub_element(self, L, ir):
        i = L.Symbol("i")
        classnames = ir["create_sub_element"] # FIXME: ffc provides element ids, not classname
        return generate_return_new_switch(L, i, classnames)
