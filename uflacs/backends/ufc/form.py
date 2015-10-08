
from uflacs.backends.ufc.generator import ufc_generator, integral_name_templates, ufc_integral_types
from uflacs.backends.ufc.templates import form_header, form_implementation
from uflacs.backends.ufc.utils import generate_return_new_switch

def add_ufc_form_integral_methods(cls):
    """This function generates methods on the class it decorates, for each integral type.

    This allows implementing e.g. create_###_integrals once in the decorated class,
    while
    """
    # The dummy name "foo" is chosen for familiarity for ffc developers
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

    def original_coefficient_position(self, L, ir):
        i = L.Symbol("i")

        positions = ir["original_coefficient_position"]

        position = L.Symbol("position")

        # Throwing a lot into the 'typename' string here but no plans for building a full C++ type system
        typename = "static const std::vector<std::size_t>"
        initializer_list = L.VerbatimExpr("{" + ", ".join(str(i) for i in positions) + "}")
        code = L.StatementList([
            L.VariableDecl(typename, position, value=initializer_list),
            L.Return(position[i]),
            ])
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
        """Return implementation of ufc::form::%(declname)s().

        This is extended to a number of functions by add_ufc_form_integral_methods.
        """
        value = ir[declname]
        return L.Return(L.LiteralInt(value))

    def _has_foo_integrals(self, L, ir, integral_type, declname):
        """Return implementation of ufc::form::%(declname)s().

        This is extended to a number of functions by add_ufc_form_integral_methods.
        """
        value = ir[declname]
        return L.Return(L.LiteralBool(value))

    def _create_foo_integral(self, L, ir, integral_type, declname):
        """Return implementation of ufc::form::%(declname)s().

        This is extended to a number of functions by add_ufc_form_integral_methods.
        """
        subdomain_id = L.Symbol("subdomain_id")
        classnames = ir[declname] # FIXME: ffc provides element id, not classname
        return generate_return_new_switch(L, subdomain_id, classnames)

    def _create_default_foo_integral(self, L, ir, integral_type, declname):
        """Return implementation of ufc::form::%(declname)s().

        This is extended to a number of functions by add_ufc_form_integral_methods.
        """
        classname = ir[declname] # FIXME: ffc provides element id, not classname
        if classname:
            return L.Return(L.New(classname))
        else:
            return L.Return(L.Null())
