#!/usr/bin/env python

"""
"""

# These are thin wrappers on top of unittest.TestCase and unittest.main
from ufltestcase import UflTestCase, main
from uflacs.codeutils.format_code import *

class CodeUtilsTest(UflTestCase):

    def test_format_code_basics(self):
        # Reproduce a string
        self.assertEqual(format_code("string"), "string")

        # Insert keywords in string using <2.6 syntax
        self.assertEqual(format_code("%(string)s %(other)s",
                                keywords=dict(string='Hello', other='World')),
                    "Hello World")

        # Basic strings with indentation
        self.assertEqual(format_code("string", 1), "    string")
        self.assertEqual(format_code("string", 2, indentchar='xy'), "xyxystring")

        # Multiline strings with indentation
        self.assertEqual(format_code("fee\nfie\nfoe", level=1, indentchar='q'),
                    "qfee\nqfie\nqfoe")

        # One item lists are the same as strings
        self.assertEqual(format_code(["hello"]), "hello")
        self.assertEqual(format_code(["hello"], 1), "    hello")
        self.assertEqual(format_code(["hello"], 2, indentchar='xy'), "xyxyhello")

        # Tuples are concatenated directly
        self.assertEqual(format_code(("hello", "world")),
                    "helloworld")
        self.assertEqual(format_code(("hello", "world"), 1),
                    "    helloworld")
        self.assertEqual(format_code(("hello", "world"), 2, indentchar='xy'),
                    "xyxyhelloworld")

        # Lists are joined with newlines
        self.assertEqual(format_code(["hello", "world"]),
                    "hello\nworld")
        self.assertEqual(format_code(["hello", "world"], 1),
                    "    hello\n    world")
        self.assertEqual(format_code(["hello", "world"], 2, indentchar='xy'),
                    "xyxyhello\nxyxyworld")

    def test_format_code_blocks(self):
        # Strings and lists can be put in Indented containers
        self.assertEqual(format_code(Indented("hei")), "    hei")
        self.assertEqual(format_code(Indented(["hei", Indented("verden")]), indentchar='z'),
                    "zhei\nzzverden")

        # A Block is an indented body with brackets before and after
        self.assertEqual(format_code(["{", Indented("fee\nfie\nfoe"), "}"], indentchar='\t'),
                    "{\n\tfee\n\tfie\n\tfoe\n}")
        self.assertEqual(format_code(Block("fee\nfie\nfoe"), indentchar='\t'),
                    "{\n\tfee\n\tfie\n\tfoe\n}")
        # A Namespace is a 'namespace foo' line before a block
        self.assertEqual(format_code(Namespace("bar", "fee\nfie\nfoe"), indentchar='\t'),
                    "namespace bar\n{\n\tfee\n\tfie\n\tfoe\n}")

        # Making a for loop
        self.assertEqual(format_code(["for (iq...)", Block("foo(iq);")]),
                    "for (iq...)\n{\n    foo(iq);\n}")
        # Making a do loop
        self.assertEqual(format_code(["iq = 0;", "do", (Block("foo(iq);"), " while (iq < nq);")]),
                    "iq = 0;\ndo\n{\n    foo(iq);\n} while (iq < nq);")

    def test_format_code_keywords(self):
        # Making a for loop with keywords
        self.assertEqual(format_code(["for (%(i)s=0; %(i)s<n; ++%(i)s)", Block("foo(%(i)s);")], keywords={'i': 'k'}),
                    "for (k=0; k<n; ++k)\n{\n    foo(k);\n}")
        self.assertEqual(format_code(WithKeywords(["for (%(i)s=0; %(i)s<n; ++%(i)s)", Block("foo(%(i)s);")],
                                             keywords={'i': 'k'})),
                    "for (k=0; k<n; ++k)\n{\n    foo(k);\n}")

        # Formatting the same code with different keywords using WithKeywords
        tmp = ['%(a)s', '%(b)s']
        self.assertEqual(format_code([WithKeywords(tmp, keywords={'a': 'a', 'b':'b'}),
                                 WithKeywords(tmp, keywords={'a': 'x', 'b':'y'})]),
                    "a\nb\nx\ny")

    def test_format_code_class(self):
        # Making a class declaration
        self.assertEqual(format_code(Class('Car')), 'class Car\n{\n};')
        self.assertEqual(format_code(Class('Car', superclass='Vehicle')),
                    'class Car: public Vehicle\n{\n};')
        self.assertEqual(format_code(Class('Car', public_body='void eval()\n{\n}')),
                    'class Car\n{\npublic:\n    void eval()\n    {\n    }\n};')
        self.assertEqual(format_code(Class('Car', protected_body='void eval()\n{\n}')),
                    'class Car\n{\nprotected:\n    void eval()\n    {\n    }\n};')
        self.assertEqual(format_code(Class('Car', private_body='void eval()\n{\n}')),
                    'class Car\n{\nprivate:\n    void eval()\n    {\n    }\n};')

    def test_format_code_templates(self):
        def t(args, mlcode, slcode):
            self.assertEqual(format_code(TemplateArgumentList(args, False), indentchar=' '), slcode)
            self.assertEqual(format_code(TemplateArgumentList(args, True), indentchar=' '), mlcode)
        t(('A',), '<\n A\n>', '<A>')
        t(('A', 'B'), '<\n A,\n B\n>', '<A, B>')

    def xtest_format_code_template_class(self):
        code = []
        code.append(('template ', TemplateArgumentList(('typename T', 'typename R'), False)))
        code.append(Class('MyClass', public_body='void hello(int world) {}'))
        #print format_code(code)


if __name__ == "__main__":
    main()
