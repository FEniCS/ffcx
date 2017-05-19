# -*- coding: utf-8 -*-
"""
Tests of generic code formatting utilities,
which focus on the overall structure of the code,
e.g. indentation, curly brace blocks, control flow
structures like loops, and class and function definitions.
Some of this is C++ specific, some is more generic.
"""

import pytest
from ffc.uflacs.language.cnodes import *


def test_format_basics():
    # Reproduce a string
    assert format_indented_lines("string") == "string"

    # Basic strings with indentation
    assert format_indented_lines("string", 1) == "    string"
    assert format_indented_lines("string", 2) == "        string"

    # Multiline strings with indentation
    assert format_indented_lines("fee\nfie\nfoe", level=1) == "    fee\n    fie\n    foe"

    # One item lists are the same as strings
    assert format_indented_lines(["hello"]) == "hello"
    assert format_indented_lines(["hello"], 1) == "    hello"
    assert format_indented_lines(["hello"], 2) == "        hello"

    # Lists are joined with newlines
    assert format_indented_lines(["hello", "world"]) == "hello\nworld"
    assert format_indented_lines(["hello", "world"], 1) == "    hello\n    world"
    assert format_indented_lines(["hello", "world"], 2) == "        hello\n        world"

    # Tuples are joined with newlines
    assert format_indented_lines(("hello", "world")) == "hello\nworld"
    assert format_indented_lines(("hello", "world"), 1) == "    hello\n    world"
    assert format_indented_lines(("hello", "world"), 2) == "        hello\n        world"
    code = format_indented_lines(("hello", "world"), 2)
    assert code == "        hello\n        world"


def test_format_indented_lines():
    # Strings and lists can be put in Indented containers
    assert format_indented_lines(Indented("hei")) == "    hei"
    code = format_indented_lines(Indented(["hei", Indented("verden")]))
    assert code == "    hei\n        verden"
    code = format_indented_lines(["{", Indented("fee\nfie\nfoe"), "}"])
    assert code == "{\n    fee\n    fie\n    foe\n}"


def xtest_format_blocks():
    # A Scope is an indented body with brackets before and after
    code = format_code(Scope("fee\nfie\nfoe"))
    assert code == "{\n    fee\n    fie\n    foe\n}"
    # A Namespace is a 'namespace foo' line before a block
    code = format_code(Namespace("bar", "fee\nfie\nfoe"))
    assert code == "namespace bar\n{\n    fee\n    fie\n    foe\n}"

    # Making a for loop
    code = format_code(["for (iq...)", Scope("foo(iq);")])
    assert code == "for (iq...)\n{\n    foo(iq);\n}"
    # Making a do loop
    code = format_code(["iq = 0;", "do", (Scope("foo(iq);"), " while (iq < nq);")])
    assert code == "iq = 0;\ndo\n{\n    foo(iq);\n} while (iq < nq);"


def xtest_format_class():
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


def xtest_format_template_argument_list():
    def t(args, mlcode, slcode):
        code = format_code(TemplateArgumentList(args, False))
        assert code == slcode
        code = format_code(TemplateArgumentList(args, True))
        assert code == mlcode
    t(('A',), '<\n    A\n>', '<A>')
    t(('A', 'B'), '<\n    A,\n    B\n>', '<A, B>')


def xtest_format_templated_type():
    code = format_code(Type('Foo'))
    assert code == 'Foo'
    code = format_code(Type('Foo', ('int', '123')))
    assert code == 'Foo<int, 123>'
    code = format_code(Type('Foo', ('int', Type('Bar', ('123', 'float')))))
    assert code == 'Foo<int, Bar<123, float> >'


def xtest_format_typedef():
    assert format_code(TypeDef('int', 'myint')) == 'typedef int myint;'
    code = format_code(TypeDef(Type('Foo', ('int', Type('Bar', ('123', 'float')))), 'Thing'))
    assert code == 'typedef Foo<int, Bar<123, float> > Thing;'


def xtest_format_template_class():
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
    assert str(code) == expected


def test_literal_cexpr_value_conversion():
    assert bool(LiteralBool(True)) is True
    assert bool(LiteralBool(False)) is False
    assert int(LiteralInt(2)) == 2
    assert float(LiteralFloat(2.0)) == 2.0
    # assert complex(LiteralFloat(2.0+4.0j)) == 2.0+4.0j


def test_format_array_decl():
    expected = "double foo[3];"
    code = ArrayDecl("double", "foo", 3)
    assert str(code) == expected
    code = ArrayDecl("double", "foo", (3,))
    assert str(code) == expected
    decl = code

    expected = "foo[1]"
    code = ArrayAccess(decl, 1)
    assert str(code) == expected
    code = ArrayAccess(decl, (1,))
    assert str(code) == expected

    expected = "foo[0]"
    code = ArrayAccess(decl, 0)
    assert str(code) == expected
    with pytest.raises(ValueError):
        ArrayAccess(decl, 3)

    code = ArrayDecl("double", "foo", (3, 4))
    expected = "double foo[3][4];"
    assert str(code) == expected
    decl = code

    expected = "foo[0][1]"
    code = ArrayAccess(decl, (0, 1))
    assert str(code) == expected
    expected = "foo[2][3]"
    code = ArrayAccess(decl, (2, 3))
    assert str(code) == expected

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
    code = ArrayDecl("double", "foo", 3, [1.0, 2.0, 3.0])
    assert str(code) == expected

    expected = """double foo[2][3] =
    { { 1.0, 2.0, 3.0 },
      { 6.0, 5.0, 4.0 } };"""
    code = ArrayDecl("double", "foo", (2, 3), [[1.0, 2.0, 3.0], [6.0, 5.0, 4.0]])
    assert str(code) == expected

    expected = """double foo[2][2][3] =
    { { { 1.0, 2.0, 3.0 },
        { 6.0, 5.0, 4.0 } },
      { { 1.0, 2.0, 3.0 },
        { 6.0, 5.0, 4.0 } } };"""
    code = ArrayDecl("double", "foo", (2, 2, 3), [[[1.0, 2.0, 3.0], [6.0, 5.0, 4.0]],
                                                  [[1.0, 2.0, 3.0], [6.0, 5.0, 4.0]]])
    assert str(code) == expected


def test_format_array_def_zero():
    expected = "double foo[3] = {};"
    code = ArrayDecl("double", "foo", 3, [0.0, 0.0, 0.0])
    assert str(code) == expected
    code = ArrayDecl("double", "foo", (3,), [0.0, 0.0, 0.0])
    assert str(code) == expected

    expected = "double foo[2][3] = {};"
    code = ArrayDecl("double", "foo", (2, 3), [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    assert str(code) == expected


def xtest_class_with_arrays():
    adecl = ArrayDecl('double', 'phi', (3, 5))
    aacc = ArrayAccess(adecl, (2, 1))
    expected = "phi[2][1]"
    assert str(aacc) == expected

    decls = [VariableDecl('double', 'foo'),
             VariableDecl('int', 'bar'),
             adecl,
             ]
    classcode = Class('Car', public_body=decls)

    expected_decls = [str(decl) for decl in decls]
    expected = 'class Car\n{\npublic:\n    %s\n};' % '\n    '.join(expected_decls)
    actual = str(classcode)
    assert actual == expected


def test_class_array_access():
    vname = 'phi'
    shape = (3, 4)
    component = ('i0', 2)
    decl = ArrayDecl('double', vname, shape)
    dcode = str(decl)
    acode = str(ArrayAccess(decl, component))
    assert dcode == 'double phi[3][4];'
    assert acode == 'phi[i0][2]'


def test_while_loop():
    code = While(VerbatimExpr("--k < 3"), [])
    actual = str(code)
    expected = "while (--k < 3)\n{\n}"
    assert actual == expected

    code = While(VerbatimExpr("--k < 3"), body=["ting;", "tang;"])
    actual = str(code)
    expected = "while (--k < 3)\n{\n    ting;\n    tang;\n}"
    assert actual == expected


def test_for_loop():
    code = For("int i = 0", "i < 3", "++i", [])
    actual = str(code)
    expected = "for (int i = 0; i < 3; ++i)\n{\n}"
    assert actual == expected

    code = For("int i = 0", "i < 3", "++i", body=["ting;", "tang;"])
    actual = str(code)
    expected = "for (int i = 0; i < 3; ++i)\n{\n    ting;\n    tang;\n}"
    assert actual == expected


def test_for_range():
    code = ForRange("i", 0, 3, [])
    actual = str(code)
    expected = "for (int i = 0; i < 3; ++i)\n{\n}"
    assert actual == expected

    code = ForRange("i", 0, 3, body=["ting;", "tang;"])
    actual = str(code)
    expected = "for (int i = 0; i < 3; ++i)\n{\n    ting;\n    tang;\n}"
    assert actual == expected
