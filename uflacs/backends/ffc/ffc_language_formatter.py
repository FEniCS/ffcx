
class FFCLanguageFormatter(MultiFunction, CppFormattingRules):
    """FFC specific cpp formatter class."""
    def __init__(self, ir):
        self._entitytype = ir["entitytype"]
        self._gdim = ir["cell"].geometric_dimension()

        self._using_names = set()
        self._includes = set(("#include <cstring>",
                              "#include <cmath>"))

    def add_using(self, name):
        self._using_names.add(name)

    def add_include(self, name):
        self._includes.add(name)

    def get_using_statements(self):
        return ["using %s;" % name for name in sorted(self._using_names)]

    def get_includes(self):
        return sorted(self._includes)

    def facet_area(self, o, component=(), derivatives=(), restriction=None):
        uflacs_assert(restriction is None, "Assuming facet_area is not restricted.")
        return self.geometric_quantity(o, component, derivatives, restriction)
