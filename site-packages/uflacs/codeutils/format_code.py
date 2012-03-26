"Tools for stitching together code snippets."

def indent(text, level, indentchar='    '):
    if level == 0:
        return text
    ind = indentchar*level
    return '\n'.join(ind + line for line in text.split('\n'))

class Code:
    pass

class Indented(Code):
    def __init__(self, code):
        self.code = code

    def format(self, level, indentchar, keywords):
        return format_code(self.code, level+1, indentchar, keywords)

class WithKeywords(Code):
    def __init__(self, code, keywords):
        self.code = code
        self.keywords = keywords

    def format(self, level, indentchar, keywords):
        return format_code(self.code, level, indentchar, self.keywords)

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
        start, sep, end = self.multilineseparators if self.multiline else self.singlelineseparators
        container = Indented if self.multiline else tuple
        if not self.multiline:
            last = self.args[-1]
            if isinstance(last, TemplateArgumentList) or (isinstance(last, Type) and last.template_arguments):
                end = ' ' + end
        code = [sep.join(format_code(arg, keywords=keywords) for arg in self.args)]
        code = (start, container(code), end)
        return format_code(code, level, indentchar, keywords)

class Type(Code):

    def __init__(self, name, template_arguments=None, ml=False):
        self.name = name
        self.template_arguments = template_arguments
        self.ml = ml

    def format(self, level, indentchar, keywords):
        code = self.name
        if self.template_arguments:
            code = code, TemplateArgumentList(self.template_arguments, self.ml)
        return format_code(code, level, indentchar, keywords)

class TypeDef(Code):
    def __init__(self, type_, typedef):
        self.type = type_
        self.typedef = typedef

    def format(self, level, indentchar, keywords):
        code = ('typedef ',self.type," %s;" % self.typedef)
        return format_code(code, level, indentchar, keywords)

class Namespace(Code):
    def __init__(self, name, body):
        self.name = name
        self.body = body

    def format(self, level, indentchar, keywords):
        code = ['namespace %s' % self.name, Block(self.body)]
        return format_code(code, level, indentchar, keywords)

class Class(Code):
    def __init__(self, name, superclass=None, public_body=None, protected_body=None, private_body=None):
        self.name = name
        self.superclass = superclass
        self.public_body = public_body
        self.protected_body = protected_body
        self.private_body = private_body

    def format(self, level, indentchar, keywords):
        if self.superclass:
            code = ['class %s: public %s' % (self.name, self.superclass)]
        else:
            code = ['class %s' % self.name]
        code += ['{']
        if self.public_body:
            code += ['public:', Indented(self.public_body)]
        if self.protected_body:
            code += ['protected:', Indented(self.protected_body)]
        if self.private_body:
            code += ['private:', Indented(self.private_body)]
        code += ['};']
        return format_code(code, level, indentchar, keywords)

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

    - Class: Piece together a class.

    See the respective classes for usage.
    """
    if isinstance(code, str):
        if keywords:
            code = code % keywords
        return indent(code, level, indentchar)

    if isinstance(code, list):
        return "\n".join(format_code(item, level, indentchar, keywords) for item in code)

    if isinstance(code, tuple):
        joined = "".join(format_code(item, 0, indentchar) for item in code)
        return format_code(joined, level, indentchar, keywords)

    if isinstance(code, Code):
        return code.format(level, indentchar, keywords)

