
from uflacs.backends.ufc.generator import ufc_generator, integral_name_templates, ufc_integral_types
import uflacs.backends.ufc.templates

class ufc_integral(ufc_generator):
    def __init__(self, integral_type):
        assert integral_type in ufc_integral_types
        integral_header = getattr(uflacs.backends.ufc.templates, "%s_integral_header" % integral_type)
        integral_implementation = getattr(uflacs.backends.ufc.templates, "%s_integral_implementation" % integral_type)
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
        # FIXME: This is where the current ffc backend code generation should be injected
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
