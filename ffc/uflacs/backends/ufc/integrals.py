# -*- coding: utf-8 -*-

from ffc.uflacs.backends.ufc.generator import ufc_generator, ufc_integral_types

class ufc_integral(ufc_generator):
    "Each function maps to a keyword in the template. See documentation of ufc_generator."
    def __init__(self, integral_type):
        assert integral_type in ufc_integral_types
        ufc_generator.__init__(self, integral_type + "_integral")

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


    def tabulate_tensor_comment(self, L, ir):
        "Generate comment for tabulate_tensor."

        r = ir["representation"]
        integrals_metadata = ir["integrals_metadata"]
        integral_metadata = ir["integral_metadata"]

        lines  = [
            "This function was generated using '%s' representation" % r,
            "with the following integrals metadata:",
            ""]
        lines += [("  " + l) for l in dstr(integrals_metadata).split("\n")]
        lines += [""]
        for i, metadata in enumerate(integral_metadata):
            lines += [
                "",
                "and the following integral %d metadata:" % i,
                "",
                ]
            lines += [("  " + l) for l in dstr(metadata).split("\n")][:-1]
            lines += [""]

        return [L.Comment(line) for line in lines]


class ufc_cell_integral(ufc_integral):
    def __init__(self):
        ufc_integral.__init__(self, "cell")

class ufc_exterior_facet_integral(ufc_integral):
    def __init__(self):
        ufc_integral.__init__(self, "exterior_facet")

class ufc_interior_facet_integral(ufc_integral):
    def __init__(self):
        ufc_integral.__init__(self, "interior_facet")

class ufc_vertex_integral(ufc_integral):
    def __init__(self):
        ufc_integral.__init__(self, "vertex")

class ufc_custom_integral(ufc_integral):
    def __init__(self):
        ufc_integral.__init__(self, "custom")

    def num_cells(self, L, ir):
        value = ir["num_cells"]
        return L.Return(L.LiteralInt(value))

class ufc_interface_integral(ufc_integral):
    def __init__(self):
        ufc_integral.__init__(self, "interface")

class ufc_overlap_integral(ufc_integral):
    def __init__(self):
        ufc_integral.__init__(self, "overlap")

class ufc_cutcell_integral(ufc_integral):
    def __init__(self):
        ufc_integral.__init__(self, "cutcell")
