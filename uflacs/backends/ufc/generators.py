
import inspect
import re
from string import Formatter

from ufl import product
from ffc.log import error, warning
from ffc.backends.ufc import *

from uflacs.language.format_lines import format_indented_lines
from uflacs.backends.ufc.templates import *

#__all__ = (["ufc_form", "ufc_dofmap", "ufc_finite_element", "ufc_integral"]
#           + ["ufc_%s_integral" % integral_type for integral_type in integral_types])

# These are all the integral types directly supported in ufc.
# TODO: Get these from somewhere for more automation.
ufc_integral_types = ("cell", "exterior_facet", "interior_facet", "vertex", "custom")

# These are the method names in ufc::form that are specialized for each integral type
integral_name_templates = (
    "max_%s_subdomain_id",
    "has_%s_integrals",
    "create_%s_integral",
    "create_default_%s_integral",
    )

# TODO: Move to language utils
def generate_return_new_switch(L, i, classnames):
    if classnames:
        cases = []
        for j, classname in enumerate(classnames):
            if classname:
                cases.append((j, L.Return(L.New(classname))))
        code = [L.Switch(i, cases, autobreak=False, autoscope=False)]
    else:
        code = []
    code.append(L.Return(L.Null()))
    return L.StatementList(code)


class ufc_generator(object):
    """Common functionality for code generators producing ufc classes.

    The generate function is the driver for generating code for a class.
    It automatically extracts template keywords and inserts the results
    from calls to self.<keyword>(language, ir), or the value of ir[keyword]
    if there is no self.<keyword>.
    """
    def __init__(self, header_template, implementation_template):
        self._header_template = header_template
        self._implementation_template = implementation_template

        r = re.compile(r"%\(([a-zA-Z0-9_]*)\)")
        self._header_keywords = set(r.findall(self._header_template))
        self._implementation_keywords = set(r.findall(self._implementation_template))

        self._keywords = sorted(self._header_keywords | self._implementation_keywords)

    def generate_snippets(self, L, ir):
        # Generate code snippets for each keyword found in templates
        snippets = {}
        for kw in self._keywords:
            # Try self.<keyword>(L, ir) if available, otherwise use ir[keyword]
            if hasattr(self, kw):
                method = getattr(self, kw)
                value = method(L, ir)
                if isinstance(value, L.CStatement):
                    value = L.Indented(value.cs_format())
                    value = format_indented_lines(value)
            else:
                value = ir.get(kw)
                if value is None:
                    error("Missing template keyword '%s' in ir for '%s'." % (kw, self.__class__.__name__))
            snippets[kw] = value

        # Error checking (can detect some bugs early when changing the interface)
        valueonly = {"classname"}
        attrs = set(name for name in dir(self) if not name.startswith("_"))
        base_attrs = set(name for name in dir(ufc_generator) if not name.startswith("_"))
        base_attrs.add("generate_snippets")
        base_attrs.add("generate")
        unused = attrs - set(self._keywords) - base_attrs
        missing = set(self._keywords) - attrs - valueonly
        if unused:
            warning("*** Unused generator functions:\n%s" % ('\n'.join(map(str,sorted(unused))),))
        if missing:
            warning("*** Missing generator functions:\n%s" % ('\n'.join(map(str,sorted(missing))),))

        return snippets

    def generate(self, L, ir):
        # Return composition of templates with generated snippets
        snippets = self.generate_snippets(L, ir)
        h = self._header_template % snippets
        cpp = self._implementation_template % snippets
        return h, cpp

    def members(self, L, ir):
        "Override in classes that need members."
        assert not ir.get("")
        return "members"

    def constructor(self, L, ir):
        "Override in classes that need constructor."
        assert not ir.get("constructor")
        return ""

    def constructor_arguments(self, L, ir):
        "Override in classes that need constructor."
        assert not ir.get("constructor_arguments")
        return ""

    def initializer_list(self, L, ir):
        "Override in classes that need constructor."
        assert not ir.get("")
        return ""

    def destructor(self, L, ir):
        "Override in classes that need destructor."
        assert not ir.get("destructor")
        return ""

    def signature(self, L, ir):
        "Default implementation of returning signature string fetched from ir."
        value = ir["signature"]
        return L.Return(L.LiteralString(value))

    def topological_dimension(self, L, ir):
        "Default implementation of returning topological dimension fetched from ir."
        value = ir["topological_dimension"]
        return L.Return(L.LiteralInt(value))

    def geometric_dimension(self, L, ir):
        "Default implementation of returning geometric dimension fetched from ir."
        value = ir["geometric_dimension"]
        return L.Return(L.LiteralInt(value))

    def create(self, L, ir):
        "Default implementation of creating a new object of the same type."
        classname = ir["classname"]
        return L.Return(L.New(classname))

def add_ufc_form_integral_methods(cls):
    """This function generates methods on the class it decorates, for each integral type.

    This allows implementing e.g. create_###_integrals once in the decorated class,
    while
    """
    # The name "foo" is chosen for familiarity for ffc developers
    impl_type = "foo"

    for template in integral_name_templates:
        implname = "_" + (template % (impl_type,))
        impl = getattr(cls, implname)
        for integral_type in ufc_integral_types:
            declname = template % (integral_type,)

            # Binding variables explicitly because Python closures don't
            # capture the value of integral_type for each iteration here
            def _delegate(self, L, ir, integral_type=integral_type, declname=declname, impl=impl):
                return impl(self, L, ir, integral_type, declname)
            _delegate.__doc__ = impl.__doc__ % {"declname": declname, "integral_type": integral_type}

            setattr(cls, declname, _delegate)
    return cls

@add_ufc_form_integral_methods
class ufc_form(ufc_generator):
    def __init__(self):
        ufc_generator.__init__(self, form_header, form_implementation)

    def num_coefficients(self, L, ir):
        value = ir["num_coefficients"]
        return L.Return(L.LiteralInt(value))

    def rank(self, L, ir):
        value = ir["rank"]
        return L.Return(L.LiteralInt(value))

    # TODO: missing 's' in ufc signature:
    def original_coefficient_position(self, L, ir): # FIXME: port this
        # Input args
        i = L.Symbol("i")
        positions = ir["original_coefficient_positions"]
        code = "FIXME"
        return code

    def create_coordinate_finite_element(self, L, ir):
        classname = ir["create_coordinate_finite_element"] # FIXME: ffc provides element id, not classname
        return L.Return(L.New(classname))
        # TODO: Use factory functions instead, here and in all create_* functions:
        #classname = ir["coordinate_finite_element_classname"] # Not in FFC
        #factoryname = make_factory_function_name(classname)
        #return L.Return(L.Call(factoryname))

    def create_coordinate_dofmap(self, L, ir):
        classname = ir["create_coordinate_dofmap"] # FIXME: ffc provides element id, not classname
        return L.Return(L.New(classname))

    def create_finite_element(self, L, ir):
        i = L.Symbol("i")
        classnames = ir["create_finite_element"] # FIXME: ffc provides element id, not classname
        return generate_return_new_switch(L, i, classnames)

    def create_dofmap(self, L, ir):
        i = L.Symbol("i")
        classnames = ir["create_dofmap"] # FIXME: ffc provides element id, not classname
        return generate_return_new_switch(L, i, classnames)

    def _max_foo_subdomain_id(self, L, ir, integral_type, declname):
        "Return implementation of ufc::form::%(declname)s()."
        value = ir[declname]
        return L.Return(L.LiteralInt(value))

    def _has_foo_integrals(self, L, ir, integral_type, declname):
        "Return implementation of ufc::form::%(declname)s()."
        value = ir[declname]
        return L.Return(L.LiteralBool(value))

    def _create_foo_integral(self, L, ir, integral_type, declname):
        "Return implementation of ufc::form::%(declname)s()."
        subdomain_id = L.Symbol("subdomain_id")
        classnames = ir[declname] # FIXME: ffc provides element id, not classname
        return generate_return_new_switch(L, subdomain_id, classnames)

    def _create_default_foo_integral(self, L, ir, integral_type, declname):
        "Return implementation of ufc::form::%(declname)s()."
        classname = ir[declname] # FIXME: ffc provides element id, not classname
        if classname:
            return L.Return(L.New(classname))
        else:
            return L.Return(L.Null())


class ufc_dofmap(ufc_generator):
    def __init__(self):
        ufc_generator.__init__(self, dofmap_header, dofmap_implementation)

    def needs_mesh_entities(self, L, ir):
        d = L.Symbol("d")
        nme = ir["needs_mesh_entities"]
        cases = [(L.LiteralInt(dim), L.Return(L.LiteralBool(need)))
                 for dim, need in enumerate(nme)]
        default = L.Return(L.LiteralBool(False))
        return L.Switch(d, cases, default=default, autoscope=False, autobreak=False)

    def global_dimension(self, L, ir): # FIXME: port this
        value = ir["global_dimension"] # FIXME: This is not an int
        code = "FIXME"
        return code

    def num_element_dofs(self, L, ir):
        value = ir["num_element_dofs"]
        return L.Return(L.LiteralInt(value))

    def num_facet_dofs(self, L, ir):
        value = ir["num_facet_dofs"]
        return L.Return(L.LiteralInt(value))

    def num_entity_dofs(self, L, ir):
        d = L.Symbol("d")
        values = ir["num_entity_dofs"]
        cases = [(i, L.Return(L.LiteralInt(value))) for i, value in enumerate(values)]
        default = L.Return(L.LiteralInt(0))
        return L.Switch(d, cases, default=default)

    def tabulate_dofs(self, L, ir): # FIXME: port this
        code = "FIXME"
        return code

    def tabulate_facet_dofs(self, L, ir): # FIXME: port this
        code = "FIXME"
        return code

    def tabulate_entity_dofs(self, L, ir): # FIXME: port this
        code = "FIXME"
        return code

    def num_sub_dofmaps(self, L, ir):
        value = ir["num_sub_dofmaps"]
        return L.Return(L.LiteralInt(value))

    def create_sub_dofmap(self, L, ir):
        i = L.Symbol("i")
        classnames = ir["create_sub_dofmap"] # FIXME: ffc provides element ids, not classname
        return generate_return_new_switch(L, i, classnames)


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


class ufc_integral(ufc_generator):
    def __init__(self, integral_type):
        assert integral_type in ufc_integral_types
        integral_header = eval("%s_integral_header" % integral_type)
        integral_implementation = eval("%s_integral_implementation" % integral_type)
        ufc_generator.__init__(self, integral_header, integral_implementation)

    def enabled_coefficients(self, L, ir):
        enabled_coefficients = ir["enabled_coefficients"]
        initializer_list = ", ".join("true" if enabled else "false"
                                     for enabled in enabled_coefficients)
        code = L.StatementList([
            # Cheating a bit with verbatim:
            L.VerbatimStatement("static const std::vector<bool> enabled({%s});" % initializer_list),
            L.Return(L.Symbol("enabled")),
            ])
        return code

    def tabulate_tensor(self, L, ir):
        # FIXME: This is where the current ffc code generation goes
        tt = ir["tabulate_tensor"]
        code = "code generated from %s" % tt
        return code

class ufc_cell_integral(ufc_integral):
    def __init__(self):
        ufc_integral.__init__(self, "cell")

class ufc_exterior_facet_integral(ufc_integral):
    def __init__(self):
        ufc_integral.__init__(self, "exterior_facet")

class ufc_interior_facet_integral(ufc_integral):
    def __init__(self):
        ufc_integral.__init__(self, "interior_facet")

class ufc_custom_integral(ufc_integral):
    def __init__(self):
        ufc_integral.__init__(self, "custom")

    def num_cells(self, L, ir):
        value = ir["num_cells"]
        return L.Return(L.LiteralInt(value))

class ufc_vertex_integral(ufc_integral):
    def __init__(self):
        ufc_integral.__init__(self, "vertex")
