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

passes = 0
fails = 0
def assertEqual(a, b):
    global passes
    global fails
    if a == b:
        passes += 1
    else:
        fails += 1
        print "Not equal:\n%s\n%s" % (a, b)

def test_format_code():
    # Reproduce a string
    assertEqual(format_code("string"), "string")

    # Insert keywords in string using <2.6 syntax
    assertEqual(format_code("%(string)s %(other)s",
                            keywords=dict(string='Hello', other='World')),
                "Hello World")

    # Basic strings with indentation
    assertEqual(format_code("string", 1), "    string")
    assertEqual(format_code("string", 2, indentchar='xy'), "xyxystring")

    # Multiline strings with indentation
    assertEqual(format_code("fee\nfie\nfoe", level=1, indentchar='q'),
                "qfee\nqfie\nqfoe")

    # One item lists are the same as strings
    assertEqual(format_code(["hello"]), "hello")
    assertEqual(format_code(["hello"], 1), "    hello")
    assertEqual(format_code(["hello"], 2, indentchar='xy'), "xyxyhello")

    # Tuples are concatenated directly
    assertEqual(format_code(("hello", "world")),
                "helloworld")
    assertEqual(format_code(("hello", "world"), 1),
                "    helloworld")
    assertEqual(format_code(("hello", "world"), 2, indentchar='xy'),
                "xyxyhelloworld")

    # Lists are joined with newlines
    assertEqual(format_code(["hello", "world"]),
                "hello\nworld")
    assertEqual(format_code(["hello", "world"], 1),
                "    hello\n    world")
    assertEqual(format_code(["hello", "world"], 2, indentchar='xy'),
                "xyxyhello\nxyxyworld")

    # Strings and lists can be put in Indented containers
    assertEqual(format_code(Indented("hei")), "    hei")
    assertEqual(format_code(Indented(["hei", Indented("verden")]), indentchar='z'),
                "zhei\nzzverden")

    # A Block is an indented body with brackets before and after
    assertEqual(format_code(["{", Indented("fee\nfie\nfoe"), "}"], indentchar='\t'),
                "{\n\tfee\n\tfie\n\tfoe\n}")
    assertEqual(format_code(Block("fee\nfie\nfoe"), indentchar='\t'),
                "{\n\tfee\n\tfie\n\tfoe\n}")
    # A Namespace is a 'namespace foo' line before a block
    assertEqual(format_code(Namespace("bar", "fee\nfie\nfoe"), indentchar='\t'),
                "namespace bar\n{\n\tfee\n\tfie\n\tfoe\n}")

    # Making a for loop
    assertEqual(format_code(["for (iq...)", Block("foo(iq);")]),
                "for (iq...)\n{\n    foo(iq);\n}")
    # Making a do loop
    assertEqual(format_code(["iq = 0;", "do", (Block("foo(iq);"), " while (iq < nq);")]),
                "iq = 0;\ndo\n{\n    foo(iq);\n} while (iq < nq);")
    # Making a for loop with keywords
    assertEqual(format_code(["for (%(i)s=0; %(i)s<n; ++%(i)s)", Block("foo(%(i)s);")], keywords={'i': 'k'}),
                "for (k=0; k<n; ++k)\n{\n    foo(k);\n}")
    assertEqual(format_code(WithKeywords(["for (%(i)s=0; %(i)s<n; ++%(i)s)", Block("foo(%(i)s);")],
                                         keywords={'i': 'k'})),
                "for (k=0; k<n; ++k)\n{\n    foo(k);\n}")

    # Formatting the same code with different keywords using WithKeywords
    tmp = ['%(a)s', '%(b)s']
    assertEqual(format_code([WithKeywords(tmp, keywords={'a': 'a', 'b':'b'}),
                             WithKeywords(tmp, keywords={'a': 'x', 'b':'y'})]),
                "a\nb\nx\ny")

    # Making a class declaration
    assertEqual(format_code(Class('Car')), 'class Car\n{\n};')
    assertEqual(format_code(Class('Car', superclass='Vehicle')),
                'class Car: public Vehicle\n{\n};')
    assertEqual(format_code(Class('Car', public_body='void eval()\n{\n}')),
                'class Car\n{\npublic:\n    void eval()\n    {\n    }\n};')
    assertEqual(format_code(Class('Car', protected_body='void eval()\n{\n}')),
                'class Car\n{\nprotected:\n    void eval()\n    {\n    }\n};')
    assertEqual(format_code(Class('Car', private_body='void eval()\n{\n}')),
                'class Car\n{\nprivate:\n    void eval()\n    {\n    }\n};')

    print "%d passes, %d fails" % (passes, fails)
