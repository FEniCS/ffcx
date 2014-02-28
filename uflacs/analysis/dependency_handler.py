
from ufl.classes import Coefficient, Argument, GeometricQuantity
from ufl.sorting import sorted_expr
from uflacs.analysis.modified_terminals import analyse_modified_terminal, analyse_modified_terminal2
from uflacs.utils.log import uflacs_assert, error

class DependencyHandler2(object):
    """Class used to collect dependencies during early compilation stages."""
    def __init__(self, modified_terminals, form_argument_mapping=None):
        # FIXME: A bit confused now, which objects are relabeled and which are not?
        # Mapping from original to relabeled form argument objects
        self.form_argument_mapping = form_argument_mapping or {}

        # FIXME: Start by building a new structure here, then use that instead of the below structures
        # FIXME: Split up in more data structures, typically need to iterate over geometry,
        #        coefficients, arguments separately in other places
        # FIXME: Change from tuple to struct of (t, c, d, r), then add avg
        # FIXME: Sort modified terminals in topological ordering w.r.t. their computational dependencies

        # ... New data structures:
        self.constant_geometry_data = []
        self.varying_geometry_data = []
        self.constant_coefficient_data = []
        self.varying_coefficient_data = []
        self.argument_data = []

        for o in sorted_expr(modified_terminals):
            # FIXME: Apply form_argument_mapping prior to getting here, then remove it here
            mt = analyse_modified_terminal2(o, form_argument_mapping)
            t = mt.terminal

            if isinstance(t, GeometricQuantity):
                if t.is_piecewise_constant():
                    self.constant_geometry_data.append(mt)
                else:
                    self.varying_geometry_data.append(mt)
            if isinstance(t, Coefficient):
                if t.is_piecewise_constant():
                    self.constant_coefficient_data.append(mt)
                else:
                    self.varying_coefficient_data.append(mt)
            elif isinstance(t, Argument):
                self.argument_data.append(mt)
            else:
                error("Unkonwn terminal type %s." % str(type(t)))

        # TODO: Sorting?
        #self.constant_geometry = sorted(self.constant_geometry, key=...)
        #self.varying_geometry = sorted(self.varying_geometry, key=...)

        self.constant_coefficient_data = sorted(self.constant_coefficient_data,
                                                key=lambda mt: (mt.terminal.count(),) + mt.key)
        self.varying_coefficient_data = sorted(self.varying_coefficient_data,
                                               key=lambda mt: (mt.terminal.count(),) + mt.key)
        self.argument_data = sorted(self.argument_data,
                                    key=lambda mt: (mt.terminal.number(),) + mt.key) # TODO: Part?

        # ...

class DependencyHandler(object):
    """Class used to collect dependencies during early compilation stages."""
    def __init__(self, modified_terminals, form_argument_mapping=None):
        # FIXME: A bit confused now, which objects are relabeled and which are not?
        # Mapping from original to relabeled form argument objects
        self.form_argument_mapping = form_argument_mapping or {}

        # Analyse modified terminals and store data about them in a canonical ordering
        self.terminal_data = [analyse_modified_terminal(o, form_argument_mapping)
                              for o in sorted_expr(modified_terminals)]

        # Extract referenced functions without modifiers and duplicates and sort by count
        self.mapped_coefficients = sorted(set(td[0] for td in self.terminal_data
                                              if isinstance(td[0], Coefficient)),
                                              key=lambda x: x.count())
        self.mapped_arguments = sorted(set(td[0] for td in self.terminal_data
                                           if isinstance(td[0], Argument)),
                                           key=lambda x: x.number())

        # TODO: This is here to make the dolfin compiler work, improve somehow?
        self.coefficient_names = dict((c,"w%d" % i) for i,c in enumerate(self.mapped_coefficients))


        # TODO: Add averaging state to (c,d,r)
        # A mapping { expr: {(c,d,r): code} } used to record visited dependencies
        self.required = {}

    def require(self, o, component, derivatives, restriction, code):
        "Helper function for remembering modified terminal dependencies."

        # FIXME: Do we get mapped functions here?
        # FIXME: Check that we have this entry in terminal_data!
        # FIXME: Record which entries we do not have in terminal_data!

        reqdata = self.required.get(o)
        if reqdata is None:
            reqdata = {}
            self.required[o] = reqdata

        avg = None # FIXME: Take as input
        key = (tuple(component), tuple(derivatives), restriction)#, avg)
        oldcode = reqdata.get(key)
        uflacs_assert((not oldcode) or (oldcode == code),
                      "Generated different code for same expression.")
        reqdata[key] = code

        return code
