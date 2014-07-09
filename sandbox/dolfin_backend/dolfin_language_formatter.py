
from six.moves import map
from ufl.common import product
from ufl.common import component_to_index
from ufl.permutation import build_component_numbering
from ufl.algorithms import MultiFunction

from ffc.log import ffc_assert

from uflacs.codeutils.cpp_format import CppFormattingRules

class DolfinExpressionLanguageFormatter(MultiFunction, CppFormattingRules):
    def __init__(self, dependency_handler, ir):
        MultiFunction.__init__(self)
        CppFormattingRules.__init__(self)

        # An object used to track who depends on what
        self._dependency_handler = dependency_handler

    def geometric_quantity(self, o, component=(), derivatives=(), restriction=None):
        "Generic rendering of variable names for all piecewise constant geometric quantities."
        ffc_assert(not derivatives,
                      "Compiler should be able to simplify derivatives of geometry.")

        # Simply using the UFL str to define the name in the generated code, ensures consistency
        name = str(o)
        if restriction:
            name = name + restriction

        # Indexing if there is a shape
        sh = o.shape()
        if sh:
            ffc_assert(component, "Missing component for nonscalar %r." % o)
            code = "%s[%d]" % (name, component_to_index(component, sh))
        else:
            ffc_assert(component == (), "Component specified for scalar %r." % o)
            code = name

        # Make a record of dependency
        self._dependency_handler.require(o, component, derivatives, restriction, code)

        return code

    def facet_area(self, o, component=(), derivatives=(), restriction=None):
        ffc_assert(restriction is None, "Assuming facet_area is not restricted.")
        return self.geometric_quantity(o, component, derivatives, restriction)

    def coefficient(self, o, component=(), derivatives=(), restriction=None):
        basename = "w%d" % o.count()
        basename = self._dependency_handler.coefficient_names[o]
        #o = self._dependency_handler.form_argument_mapping.get(o,o)#[o]
        #print self._dependency_handler.form_argument_mapping
        return self.form_argument(o, component, derivatives, restriction, basename)

    def argument(self, o, component=(), derivatives=(), restriction=None):
        raise NotImplementedError("Not expecting an Argument in dolfin expression formatter.")

    def form_argument(self, o, component, derivatives, restriction, base_name):
        # FIXME: This probably needs some work, revisit when implementing
        #        actual generation of declarations of these quantities

        # FIXME: I think we need to flatten and combine component and derivatives.

        def indstring(indices):
            # FIXME: Indexing flat or nested in C++?
            if 1:
                return "".join("[%s]" % i for i in indices)
            else:
                return ("[%s]" % (", ".join(map(str, indices))))

        rcode = {None:"", "+": "_p", "-": "_m"}[restriction]

        if derivatives:
            # FIXME: Format using derivative counts
            dcodepre  = 'd%d_' % (len(derivatives),)
            dcodepost = indstring(derivatives)
        else:
            dcodepre = "v_"
            dcodepost = ""

        icode = indstring(component) if component else "[0]"

        code = dcodepre + base_name + rcode + icode + dcodepost

        self._dependency_handler.require(o, component, derivatives, restriction, code)
        return code
