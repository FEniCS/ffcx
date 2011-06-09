
def indent(text, level, indentchar='    '):
    if level == 0:
        return text
    ind = indentchar*level
    return '\n'.join(ind + line for line in text.split('\n'))

class Indented:
    def __init__(self, code):
        self.code = code

class Block:
    def __init__(self, body, start='{', end='}'):
        self.start = start
        self.body = body
        self.end = end

class Namespace:
    def __init__(self, name, body):
        self.name = name
        self.body = body

def format_code(code, level=0, indentchar='    '):
    if isinstance(code, str):
        return indent(code, level, indentchar)

    if isinstance(code, list):
        return "\n".join(format_code(item, level, indentchar) for item in code)

    if isinstance(code, tuple):
        joined = "".join(format_code(item, 0, indentchar) for item in code)
        return format_code(joined, level, indentchar)

    if isinstance(code, Indented):
        return format_code(code.code, level+1, indentchar)

    if isinstance(code, Block):
        return format_code([code.start, Indented(code.body), code.end],
                      level, indentchar)

    if isinstance(code, Namespace):
        return format_code(['namespace %s' % code.name, Block(code.body)],
                      level, indentchar)

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
    # Basic strings with indentation
    assertEqual(format_code("string"), "string")
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

    print "%d passes, %d fails" % (passes, fails)

def tabulate_tensor_prototype():

    def build_body01():
        body = []
        body += ["decl"]
        body += ["for ()"]
        decl = ["decl01;"]
        body += [Block(decl)]
        return body

    def build_loop():
        preloop = ["preloop"]
        postloop = ["postloop"]
        loopheader = "for ()"
        body = ["body;"]
        return [preloop, loopheader, Block(body), postloop]

    def build_body1():
        body = []
        body += ["decl"]
        body += ["for (i1)"]
        decl = ["decl1;"]
        body += [Block(decl)]
        return body

    def build_body0():
        body = []
        body += ["decl"]
        body += ["for (i1)"]
        decl = [build_body1()]
        body += [Block(decl)]
        return body

    def build_bodyx():
        body = []
        body += ["compute x from iq;"]
        body += ["compute quantities depending on x"]
        body += build_i0loop()
        body += build_i1loop()
        return body

    def build_i0loop():
        header = "for (int i0 = 0; i0 < n0; ++i0)"
        decl = build_body0()
        return [header, Block(decl)]

    def build_i1loop():
        header = "for (int i1 = 0; i1 < n1; ++i1)"
        decl = build_body1()
        return [header, Block(decl)]

    def build_body():
        body = []
        body += ["compute quantities depending on cell and coeffs but not x;"]
        body += ["for (int iq = 0; iq < nq; ++iq)"]
        decl = build_bodyx()
        body += [Block(decl)]
        return body

    print format_code(build_body())

if __name__ == '__main__':
    test_format_code()
    tabulate_tensor_prototype()
