
from ufl.common import component_to_index
from ufl.permutation import build_component_numbering
from ufl.algorithms import MultiFunction

from uflacs.utils.log import uflacs_assert, warning, error

from uflacs.codeutils.cpp_format import CppFormatterRulesCollection

class ToyCppLanguageFormatter(MultiFunction, CppFormatterRulesCollection):
    """Example cpp formatter class, used for the test cases.
    Override the same functions for your particular target."""
    def __init__(self, dependency_handler, ir):
        MultiFunction.__init__(self)
        CppFormatterRulesCollection.__init__(self)

        # An object used to track who depends on what
        self._dependency_handler = dependency_handler
        self._coefficient_mapping = self._dependency_handler.form_argument_mapping

    def geometric_quantity(self, o, component=(), derivatives=(), restriction=None):
        "Generic rendering of variable names for all piecewise constant geometric quantities."
        uflacs_assert(not derivatives,
                      "Compiler should be able to simplify derivatives of geometry.")

        # Simply using the UFL str to define the name in the generated code, ensures consistency
        name = str(o)
        if restriction:
            name = name + restriction

        # Indexing if there is a shape
        sh = o.shape()
        if sh:
            uflacs_assert(component, "Missing component for nonscalar %r." % o)
            code = "%s[%d]" % (name, component_to_index(component, sh))
        else:
            uflacs_assert(component == (), "Component specified for scalar %r." % o)
            code = name

        # Make a record of dependency
        self._dependency_handler.require(o, component, derivatives, restriction, code)

        return code

    def facet_area(self, o, component=(), derivatives=(), restriction=None):
        uflacs_assert(restriction is None, "Assuming facet_area is not restricted.")
        return self.geometric_quantity(o, component, derivatives, restriction)

    def _piecewise_constant_coefficient(self, o, component, derivatives, restriction):
        o = self._coefficient_mapping.get(o, o)

        uflacs_assert(not derivatives,
                      "Not expecting derivatives of constant coefficients!")

        vi2si, si2vi = build_component_numbering(o.shape(), o.element().symmetry())
        comp = vi2si[component]

        if restriction == "-":
            comp += len(si2vi)

        return "w[%d][%d]" % (o.count(), comp)

    def coefficient(self, o, component=(), derivatives=(), restriction=None):
        o = self._coefficient_mapping.get(o, o)

        count = o.count()
        uflacs_assert(count >= 0,
            "Expecting positive count, provide a renumbered form argument mapping.")
        if o.is_cellwise_constant():
            return self._piecewise_constant_coefficient(o, component, derivatives, restriction)
        else:
            return self.form_argument(o, component, derivatives, restriction, "w%d" % count)

    def argument(self, o, component=(), derivatives=(), restriction=None):
        return self.form_argument(o, component, derivatives, restriction, "v%d" % o.number())

    def form_argument(self, o, component, derivatives, restriction, base_name):
        rcode = {None:"", "+": "_p", "-": "_m"}[restriction]

        def indstring(indices):
            # FIXME: Indexing flat or nested in C++?
            if 1:
                return "".join("[%s]" % i for i in indices)
            else:
                return ("[%s]" % (", ".join(map(str,indices))))

        if derivatives:
            dcodepre  = "d%d_" % (len(derivatives),)
            dcodepost = indstring(derivatives)
        else:
            dcodepre, dcodepost = "", ""

        if component:
            icode = indstring(component)
        else:
            icode = ""

        code = dcodepre + base_name + rcode + icode + dcodepost

        self._dependency_handler.require(o, component, derivatives, restriction, code)
        return code
