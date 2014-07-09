#!/usr/bin/env python
"""
Tests of generic code formatting utilities,
which focus on the overall structure of the code,
e.g. indentation, curly brace blocks, control flow
structures like loops, and class and function definitions.
Some of this is C++ specific, some is more generic.
"""

import pytest
from uflacs.codeutils.format_code import *

def test_format_basics():
    # Reproduce a string
    assert format_code("string") == "string"

    # Insert keywords in string using <2.6 syntax
    code = format_code("%(string)s %(other)s", keywords=dict(string='Hello', other='World'))
    assert code == "Hello World"

    # Basic strings with indentation
    assert format_code("string", 1) == "    string"
    assert format_code("string", 2, indentchar='xy') == "xyxystring"

    # Multiline strings with indentation
    assert format_code("fee\nfie\nfoe", level=1, indentchar='q') == "qfee\nqfie\nqfoe"

    # One item lists are the same as strings
    assert format_code(["hello"]) == "hello"
    assert format_code(["hello"], 1) == "    hello"
    assert format_code(["hello"], 2, indentchar='xy') == "xyxyhello"

    # Tuples are concatenated directly
    assert format_code(("hello", "world")) == "helloworld"
    assert format_code(("hello", "world"), 1) == "    helloworld"
    code = format_code(("hello", "world"), 2, indentchar='xy')
    assert code == "xyxyhelloworld"

    # Lists are joined with newlines
    assert format_code(["hello", "world"]) == "hello\nworld"
    assert format_code(["hello", "world"], 1) == "    hello\n    world"
    assert format_code(["hello", "world"], 2, indentchar='xy') == "xyxyhello\nxyxyworld"

def test_format_blocks():
    # Strings and lists can be put in Indented containers
    assert format_code(Indented("hei")) == "    hei"
    code = format_code(Indented(["hei", Indented("verden")]), indentchar='z')
    assert code == "zhei\nzzverden"

    # A Block is an indented body with brackets before and after
    code = format_code(["{", Indented("fee\nfie\nfoe"), "}"], indentchar='\t')
    assert code == "{\n\tfee\n\tfie\n\tfoe\n}"
    code = format_code(Block("fee\nfie\nfoe"), indentchar='\t')
    assert code == "{\n\tfee\n\tfie\n\tfoe\n}"
    # A Namespace is a 'namespace foo' line before a block
    code = format_code(Namespace("bar", "fee\nfie\nfoe"), indentchar='\t')
    assert code == "namespace bar\n{\n\tfee\n\tfie\n\tfoe\n}"

    # Making a for loop
    code = format_code(["for (iq...)", Block("foo(iq);")])
    assert code == "for (iq...)\n{\n    foo(iq);\n}"
    # Making a do loop
    code = format_code(["iq = 0;", "do", (Block("foo(iq);"), " while (iq < nq);")])
    assert code == "iq = 0;\ndo\n{\n    foo(iq);\n} while (iq < nq);"

def test_format_keywords():
    # Making a for loop with keywords
    code = format_code(["for (%(i)s=0; %(i)s<n; ++%(i)s)", Block("foo(%(i)s);")], keywords={'i': 'k'})
    assert code == "for (k=0; k<n; ++k)\n{\n    foo(k);\n}"
    format_code(WithKeywords(["for (%(i)s=0; %(i)s<n; ++%(i)s)",
                                        Block("foo(%(i)s);")], keywords={'i': 'k'}))
    assert code == "for (k=0; k<n; ++k)\n{\n    foo(k);\n}"

    # Formatting the same code with different keywords using WithKeywords
    tmp = ['%(a)s', '%(b)s']
    code = format_code([WithKeywords(tmp, keywords={'a': 'a', 'b':'b'}),
                                WithKeywords(tmp, keywords={'a': 'x', 'b':'y'})])
    assert code == "a\nb\nx\ny"

def test_format_class():
    # Making a class declaration
    assert format_code(Class('Car')) == 'class Car\n{\n};'
    code = format_code(Class('Car', superclass='Vehicle'))
    assert code == 'class Car: public Vehicle\n{\n};'
    code = format_code(Class('Car', public_body='void eval()\n{\n}'))
    assert code == 'class Car\n{\npublic:\n    void eval()\n    {\n    }\n};'
    code = format_code(Class('Car', protected_body='void eval()\n{\n}'))
    assert code == 'class Car\n{\nprotected:\n    void eval()\n    {\n    }\n};'
    code = format_code(Class('Car', private_body='void eval()\n{\n}'))
    assert code == 'class Car\n{\nprivate:\n    void eval()\n    {\n    }\n};'

def test_format_template_argument_list():
    def t(args, mlcode, slcode):
        code = format_code(TemplateArgumentList(args, False), indentchar=' ')
        assert code == slcode
        code = format_code(TemplateArgumentList(args, True), indentchar=' ')
        assert code == mlcode
    t(('A',), '<\n A\n>', '<A>')
    t(('A', 'B'), '<\n A,\n B\n>', '<A, B>')

def test_format_templated_type():
    code = format_code(Type('Foo'))
    assert code == 'Foo'
    code = format_code(Type('Foo', ('int', '123')))
    assert code == 'Foo<int, 123>'
    code = format_code(Type('Foo', ('int', Type('Bar', ('123', 'float')))))
    assert code == 'Foo<int, Bar<123, float> >'

def test_format_typedef():
    assert format_code(TypeDef('int', 'myint')) == 'typedef int myint;'
    code = format_code(TypeDef(Type('Foo', ('int', Type('Bar', ('123', 'float')))), 'Thing'))
    assert code == 'typedef Foo<int, Bar<123, float> > Thing;'

def test_format_template_class():
    expected = "template<typename T, typename R>\nclass MyClass\n{\npublic:\n    void hello(int world) {}\n};"

    code = [('template', TemplateArgumentList(('typename T', 'typename R'), False)),
            Class('MyClass', public_body='void hello(int world) {}')]
    assert format_code(code) == expected

    code = Class('MyClass', public_body='void hello(int world) {}',
                 template_arguments=('typename T', 'typename R'))
    assert format_code(code) == expected

def test_format_variable_decl():
    code = VariableDecl("double", "foo")
    expected = "double foo;"
    assert format_code(code) == expected

def test_format_array_decl():
    expected = "double foo[3];"
    code = ArrayDecl("double", "foo", 3)
    assert format_code(code) == expected
    code = ArrayDecl("double", "foo", (3,))
    assert format_code(code) == expected
    decl = code

    expected = "foo[1]"
    code = ArrayAccess(decl, 1)
    assert format_code(code) == expected
    code = ArrayAccess(decl, (1,))
    assert format_code(code) == expected

    expected = "foo[0]"
    code = ArrayAccess(decl, 0)
    assert format_code(code) == expected
    with pytest.raises(ValueError):
        ArrayAccess(decl, 3)

    code = ArrayDecl("double", "foo", (3, 4))
    expected = "double foo[3][4];"
    assert format_code(code) == expected
    decl = code

    expected = "foo[0][1]"
    code = ArrayAccess(decl, (0, 1))
    assert format_code(code) == expected
    expected = "foo[2][3]"
    code = ArrayAccess(decl, (2, 3))
    assert format_code(code) == expected

    with pytest.raises(ValueError):
        ArrayAccess(decl, 0)
    with pytest.raises(ValueError):
        ArrayAccess(decl, (0, 4))
    with pytest.raises(ValueError):
        ArrayAccess(decl, (3, 0))
    with pytest.raises(ValueError):
        ArrayAccess(decl, (3, -1))

def test_format_array_def():
    expected = "double foo[3] = { 1.0, 2.0, 3.0 };"
    code = ArrayDecl("double", "foo", 3, ["1.0", "2.0", "3.0"])
    assert format_code(code) == expected
    code = ArrayDecl("double", "foo", (3,), ["1.0", "2.0", "3.0"])
    assert format_code(code) == expected

    expected = "double foo[2][3] = {\n    { 1.0, 2.0, 3.0 },\n    { 6.0, 5.0, 4.0 }\n    };"
    code = ArrayDecl("double", "foo", (2, 3), [["1.0", "2.0", "3.0"], ["6.0", "5.0", "4.0"]])
    assert format_code(code) == expected

    expected = "double foo[2][2][3] = {\n        {\n        { 1.0, 2.0, 3.0 },\n        { 6.0, 5.0, 4.0 }\n        },\n        {\n        { 1.0, 2.0, 3.0 },\n        { 6.0, 5.0, 4.0 }\n        }\n    };"
    code = ArrayDecl("double", "foo", (2, 2, 3), [[["1.0", "2.0", "3.0"], ["6.0", "5.0", "4.0"]],
                                                [["1.0", "2.0", "3.0"], ["6.0", "5.0", "4.0"]]])
    if 0:
        print()
        print(expected)
        print()
        print(format_code(code))
        print()
    assert format_code(code) == expected

def test_class_with_arrays():
    adecl = ArrayDecl('double', 'phi', (3, 5))
    aacc = ArrayAccess(adecl, (2, 1))
    expected = "Car.phi[2][1]"
    code = ("Car", ".", aacc)
    assert format_code(code) == expected

    decls = [VariableDecl('double', 'foo'),
             VariableDecl('int', 'bar'),
             adecl,
             ]
    classcode = Class('Car', public_body=decls)

    expected_decls = [format_code(decl) for decl in decls]
    expected = 'class Car\n{\npublic:\n    %s\n};' % '\n    '.join(expected_decls)
    actual = format_code(classcode)
    assert actual == expected

def test_class_array_access():
    oname = 'pvars[iq]'
    vname = 'phi'
    shape = (3, 4)
    component = ('i0', 2)
    decl = ArrayDecl('double', vname, shape)
    dcode = format_code(decl)
    acode = format_code((oname, '.', ArrayAccess(decl, component)))
    assert dcode == 'double phi[3][4];'
    assert acode == 'pvars[iq].phi[i0][2]'

def test_while_loop():
    code = WhileLoop("--k < 3")
    actual = format_code(code)
    expected = "while (--k < 3)"
    assert actual == expected

    code = WhileLoop("--k < 3", body=["ting;", "tang;"])
    actual = format_code(code)
    expected = "while (--k < 3)\n{\n    ting;\n    tang;\n}"
    assert actual == expected

def test_for_loop():
    code = ForLoop("int i = 0", "i < 3", "++i")
    actual = format_code(code)
    expected = "for (int i = 0; i < 3; ++i)"
    assert actual == expected

    code = ForLoop("int i = 0", "i < 3", "++i", body=["ting;", "tang;"])
    actual = format_code(code)
    expected = "for (int i = 0; i < 3; ++i)\n{\n    ting;\n    tang;\n}"
    assert actual == expected

def test_for_range():
    code = ForRange("i", "0", "3")
    actual = format_code(code)
    expected = "for (int i = 0; i < 3; ++i)"
    assert actual == expected

    code = ForRange("i", "0", "3", body=["ting;", "tang;"])
    actual = format_code(code)
    expected = "for (int i = 0; i < 3; ++i)\n{\n    ting;\n    tang;\n}"
    assert actual == expected
