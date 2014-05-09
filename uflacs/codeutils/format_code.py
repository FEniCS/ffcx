"""
Tools for stitching together code snippets.
"""


def strip_trailing_whitespace(s):
    return '\n'.join(l.rstrip() for l in s.split('\n'))

def format_float(x):
    eps = 1e-12 # FIXME: Configurable threshold
    if abs(x) < eps:
        return "0.0"

    precision = 12 # FIXME: Configurable precision
    fmt = "%%.%de" % precision
    return fmt % x

def indent(text, level, indentchar='    '):
    if level == 0:
        return text
    ind = indentchar*level
    return '\n'.join(ind + line for line in text.split('\n'))

def build_separated_list(values, sep):
    "Make a list with sep inserted between each value in values."
    items = []
    for v in values[:-1]:
        items.append((v, sep))
    items.append(values[-1])
    return items

def build_initializer_list(values, begin="{ ", sep=", ", end=" }"):
    "Build a value initializer list."
    return [begin] + build_separated_list(values, sep) + [end]

def build_recursive_initializer_list(values, sizes):
    r = len(sizes)
    assert r > 0
    assert len(values) == sizes[0]

    if r == 1:
        initializer_list = tuple(build_initializer_list(values))

    elif r == 2:
        assert len(values[0]) == sizes[1]
        inner = []
        for i0 in range(sizes[0]):
            inner.append( tuple(build_initializer_list(values[i0])) )
        initializer_list = ["{", Indented([build_separated_list(inner, ","), "}"])]

    elif r == 3:
        assert len(values[0]) == sizes[1]
        assert len(values[0][0]) == sizes[2]
        outer = []
        for i0 in range(sizes[0]):
            inner = []
            for i1 in range(sizes[1]):
                inner.append( tuple(build_initializer_list(values[i0][i1])) )
            outer.append(Indented(["{", build_separated_list(inner, ","), "}"]))
        initializer_list = ["{", Indented([build_separated_list(outer, ","), "}"])]

    else:
        error("TODO: Make recursive implementation of initializer_list formatting.")

    return initializer_list



class Code(object):
    pass

class Indented(Code):
    def __init__(self, code):
        self.code = code

    def format(self, level, indentchar, keywords):
        return format_code(self.code, level+1, indentchar, keywords)

class WithKeywords(Code): # TODO: Do we need this? Can simplify quite a bit by removing.
    def __init__(self, code, keywords):
        self.code = code
        self.keywords = keywords

    def format(self, level, indentchar, keywords):
        if keywords: # TODO: Merge with self.keywords instead
            raise RuntimeError("Doubly defined keywords not implemented.")
        fmt_keywords = {}
        for k,v in self.keywords.iteritems():
            fmt_keywords[k] = format_code(v, 0, indentchar, keywords)
        return format_code(self.code, level, indentchar, fmt_keywords)

class Block(Code):
    def __init__(self, body, start='{', end='}'):
        self.start = start
        self.body = body
        self.end = end

    def format(self, level, indentchar, keywords):
        code = [self.start, Indented(self.body), self.end]
        return format_code(code, level, indentchar, keywords)

class TemplateArgumentList(Code):

    singlelineseparators = ('<',', ','>')
    multilineseparators = ('<\n',',\n','\n>')

    def __init__(self, args, multiline=True):
        self.args = args
        self.multiline = multiline

    def format(self, level, indentchar, keywords):
        if self.multiline:
            container = Indented
            start, sep, end = self.multilineseparators
        else:
            container = tuple
            start, sep, end = self.singlelineseparators
            # Add space to avoid >> template issue
            last = self.args[-1]
            if isinstance(last, TemplateArgumentList) or (
                isinstance(last, Type) and last.template_arguments):
                end = ' ' + end
        code = [sep.join(format_code(arg, keywords=keywords) for arg in self.args)]
        code = (start, container(code), end)
        return format_code(code, level, indentchar, keywords)

class Type(Code):
    def __init__(self, name, template_arguments=None, multiline=False):
        self.name = name
        self.template_arguments = template_arguments
        self.multiline = multiline

    def format(self, level, indentchar, keywords):
        code = self.name
        if self.template_arguments:
            code = code, TemplateArgumentList(self.template_arguments, self.multiline)
        return format_code(code, level, indentchar, keywords)

class TypeDef(Code):
    def __init__(self, type_, typedef):
        self.type_ = type_
        self.typedef = typedef

    def format(self, level, indentchar, keywords):
        code = ('typedef ', self.type_, " %s;" % self.typedef)
        return format_code(code, level, indentchar, keywords)

class Namespace(Code):
    def __init__(self, name, body):
        self.name = name
        self.body = body

    def format(self, level, indentchar, keywords):
        code = ['namespace %s' % self.name, Block(self.body)]
        return format_code(code, level, indentchar, keywords)

class VariableDecl(Code):
    def __init__(self, typename, name, value=None):
        self.typename = typename
        self.name = name
        self.value = value

    def format(self, level, indentchar, keywords):
        sep = " "
        code = (self.typename, sep, self.name)
        if self.value is not None:
            code += (" = ", self.value)
        code += (";",)
        return format_code(code, level, indentchar, keywords)

class ArrayDecl(Code):
    def __init__(self, typename, name, sizes, values=None):
        self.typename = typename
        self.name = name
        self.sizes = (sizes,) if isinstance(sizes, int) else tuple(sizes)
        self.values = values

    def format(self, level, indentchar, keywords):
        sep = " "
        brackets = tuple("[%d]" % n for n in self.sizes)
        if self.values is None:
            valuescode = ""
        else:
            initializer_list = build_recursive_initializer_list(self.values, self.sizes)
            valuescode = (" = ", initializer_list)
        code = (self.typename, sep, self.name, brackets, valuescode, ";")
        return format_code(code, level, indentchar, keywords)

class ArrayAccess(Code):
    def __init__(self, arraydecl, indices):
        if isinstance(arraydecl, ArrayDecl):
            self.arrayname = arraydecl.name
        else:
            self.arrayname = arraydecl

        if isinstance(indices, (list,tuple)):
            self.indices = indices
        else:
            self.indices = (indices,)

        # Early error checking of array dimensions
        if any(isinstance(i, int) and i < 0 for i in self.indices):
            raise ValueError("Index value < 0.")

        # Additional checks possible if we get an ArrayDecl instead of just a name
        if isinstance(arraydecl, ArrayDecl):
            if len(self.indices) != len(arraydecl.sizes):
                raise ValueError("Invalid number of indices.")
            if any((isinstance(i, int) and isinstance(d, int) and i >= d)
                   for i,d in zip(self.indices, arraydecl.sizes)):
                raise ValueError("Index value >= array dimension.")

    def format(self, level, indentchar, keywords):
        brackets = tuple(("[", n, "]") for n in self.indices)
        code = (self.arrayname, brackets)
        return format_code(code, level, indentchar, keywords)

class WhileLoop(Code):
    def __init__(self, check, body=None):
        self.check = check
        self.body = body

    def format(self, level, indentchar, keywords):
        code = ("while (", self.check, ")")
        if self.body is not None:
            code = [code, Block(self.body)]
        return format_code(code, level, indentchar, keywords)

class ForLoop(Code):
    def __init__(self, init, check, increment, body=None):
        self.init = init
        self.check = check
        self.increment = increment
        self.body = body

    def format(self, level, indentchar, keywords):
        code = ("for (", self.init, "; ", self.check, "; ", self.increment, ")")
        if self.body is not None:
            code = [code, Block(self.body)]
        return format_code(code, level, indentchar, keywords)

class ForRange(Code):
    def __init__(self, name, lower, upper, body=None):
        self.name = name
        self.lower = lower
        self.upper = upper
        self.body = body

    def format(self, level, indentchar, keywords):
        init = ("int ", self.name, " = ", self.lower)
        check = (self.name, " < ", self.upper)
        increment = ("++", self.name)
        code = ForLoop(init, check, increment, body=self.body)
        return format_code(code, level, indentchar, keywords)

class Class(Code):
    def __init__(self, name, superclass=None, public_body=None,
                 protected_body=None, private_body=None,
                 template_arguments=None, template_multiline=False):
        self.name = name
        self.superclass = superclass
        self.public_body = public_body
        self.protected_body = protected_body
        self.private_body = private_body
        self.template_arguments = template_arguments
        self.template_multiline = template_multiline

    def format(self, level, indentchar, keywords):
        code = []
        if self.template_arguments:
            code += [('template', TemplateArgumentList(self.template_arguments,
                                                       self.template_multiline))]
        if self.superclass:
            code += ['class %s: public %s' % (self.name, self.superclass)]
        else:
            code += ['class %s' % self.name]
        code += ['{']
        if self.public_body:
            code += ['public:', Indented(self.public_body)]
        if self.protected_body:
            code += ['protected:', Indented(self.protected_body)]
        if self.private_body:
            code += ['private:', Indented(self.private_body)]
        code += ['};']
        return format_code(code, level, indentchar, keywords)

class Comment(Code):
    def __init__(self, comment):
        self.comment = comment

    def format(self, level, indentchar, keywords):
        code = ("// ", self.comment)
        return format_code(code, level, indentchar, keywords)

class Return(Code):
    def __init__(self, value):
        self.value

    def format(self, level, indentchar, keywords):
        code = ("return ", self.value, ";")
        return format_code(code, level, indentchar, keywords)


class AssignBase(Code):
    def __init__(self, lhs, rhs):
        self.lhs = lhs
        self.rhs = rhs

    def format(self, level, indentchar, keywords):
        code = (self.lhs, type(self).op, self.rhs, ";")
        return format_code(code, level, indentchar, keywords)

class Assign(AssignBase):
    op = " = "

class AssignAdd(AssignBase):
    op = " += "

class AssignSub(AssignBase):
    op = " -= "

class AssignMul(AssignBase):
    op = " *= "

class AssignDiv(AssignBase):
    op = " /= "


class UnOp(Code):
    def __init__(self, arg):
        self.arg = arg

    def format(self, level, indentchar, keywords):
        # TODO: Handle precedence at this level instead of in the ExprFormatter stuff?
        code = (type(self).op, self.arg)
        return format_code(code, level, indentchar, keywords)

class BinOp(Code):
    def __init__(self, lhs, rhs):
        self.lhs = lhs
        self.rhs = rhs

    def format(self, level, indentchar, keywords):
        # TODO: Handle precedence at this level instead of in the ExprFormatter stuff?
        code = (self.lhs, type(self).op, self.rhs)
        return format_code(code, level, indentchar, keywords)

class NOp(Code):
    def __init__(self, ops):
        self.ops = ops

    def format(self, level, indentchar, keywords):
        # TODO: Handle precedence at this level instead of in the ExprFormatter stuff?
        code = []
        for op in self.ops:
            code.append(op)
            code.append(type(self).op)
        code = tuple(code[:-1])
        return format_code(code, level, indentchar, keywords)


class Add(BinOp):
    op = " + "

class Sub(BinOp):
    op = " - "

class Mul(BinOp):
    op = " * "

class Div(BinOp):
    op = " / "


class Negative(UnOp):
    op = "-"


class Sum(NOp):
    op = " + "

class Product(NOp):
    op = " * "


def format_code(code, level=0, indentchar='    ', keywords=None):
    """Format code by stitching together snippets. The code can
    be built recursively using the following types:

    - str: Just a string, keywords can be provided to replace %%(name)s.

    - tuple: Concatenate items in tuple with no chars in between.

    - list: Concatenate items in tuple with newline in between.

    - Indented: Indent the code within this object one level.

    - WithKeywords: Assign keywords to code within this object.

    - Block: Wrap code in {} and indent it.

    - Namespace: Wrap code in a namespace.

    - Class: Piece together a class definition.

    - TemplateArgumentList: Format a template argument list, one line or one per type.

    - Type: Format a typename with or without template arguments.

    - TypeDef: Format a typedef for a type.

    - VariableDecl: Declaration of a variable.

    - ArrayDecl: Declaration of an array.

    - ArrayAccess: Access element of array.

    - Return: Return statement with value.

    - Assign: = statement.

    - AssignAdd: += statement.

    - AssignSub: -= statement.

    - AssignMul: *= statement.

    - AssignDiv: /= statement.

    See the respective classes for usage.
    """
    if isinstance(code, str):
        if keywords:
            code = code % keywords # TODO: Either change to new formatting or just remove keywords
        if level:
            return indent(code, level, indentchar)
        else:
            return code

    if isinstance(code, list):
        return "\n".join(format_code(item, level, indentchar, keywords) for item in code)

    if isinstance(code, tuple):
        joined = "".join(format_code(item, 0, indentchar, keywords) for item in code)
        return format_code(joined, level, indentchar, keywords)

    if isinstance(code, Code):
        return code.format(level, indentchar, keywords)

    if isinstance(code, int):
        return indent(str(code), level, indentchar)

    if isinstance(code, float):
        return indent(format_float(code), level, indentchar)

    raise RuntimeError("Unexpected type %s:\n%s" % (type(code),str(code)))

