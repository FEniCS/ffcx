
# TODO: Replace with cgen module?
# TODO: Should some of these be implemented in the format_code design instead? Or should that design call upon this class for formatting?

class CppStatementFormattingRules(object):

    # TODO: This is a duplicate of product rule in cpp_format, placed here during refactoring.
    def product(self, a, b):
        return "%s * %s" % (a, b)

    # TODO: This is a duplicate of sum rule in cpp_format, placed here during refactoring.
    def sum(self, a, b):
        return "%s + %s" % (a, b)

    # TODO: ?
    def function_return(self, value):
        return "return %s;" % (value,)

    def for_loop(self, typename, iname, ibegin, iend):
        return "for (%s %s = %s; %s < %s; ++%s)" % (typename, iname, ibegin, iname, iend, iname)

    def comment(self, line):
        return "// %s" % (line,)

    def assign(self, lvalue, rvalue):
        return "%s = %s;" % (lvalue, rvalue)

    def iadd(self, lvalue, rvalue):
        return "%s += %s;" % (lvalue, rvalue)

    def isub(self, lvalue, rvalue):
        return "%s -= %s;" % (lvalue, rvalue)

    def imul(self, lvalue, rvalue):
        return "%s *= %s;" % (lvalue, rvalue)

    def idiv(self, lvalue, rvalue):
        return "%s /= %s;" % (lvalue, rvalue)

    def precision_float(self, f):
        return "%g" % f  # TODO: Control float formatting precision here

    def precision_floats(self, fs):
        return ', '.join(self.precision_float(f) for f in fs)

    def array_access(self, name, *keys):
        br = ''.join("[%s]" % k for k in keys)
        return "%s%s" % (name, br)

    def array_decl(self, typename, name, dims, initializer=None):
        dims = (dims,) if isinstance(dims, (int, str)) else dims
        init = " = %s" % (initializer,) if initializer else ""
        br = ''.join("[%s]" % d for d in dims)
        return "%s %s%s%s;" % (typename, name, br, init)

    def var_decl(self, typename, name, initializer=None):
        if initializer:
            init = " = %s" % (initializer,)
        return "%s %s%s;" % (typename, name, init)

    def double_array_decl(self, name, dims, initializer=None):
        return self.array_decl("double", name, dims, initializer)

    def double_decl(self, name, initializer=None):
        return self.var_decl("double", name, initializer)
