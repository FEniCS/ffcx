# -*- coding: utf-8 -*-
"""
Tests of CNode formatting.
"""
from ffc.uflacs.language.cnodes import *


def test_cnode_expression_precedence():
    assert str(Add(1, 2)) == "1 + 2"
    assert str(Add(Add(1, 2), 3)) == "1 + 2 + 3"
    assert str(Add(1, Add(2, 3))) == "1 + (2 + 3)"
    assert str(Add(Add(1, 2), Add(3, 4))) == "1 + 2 + (3 + 4)"
    assert str(Add(Add(1, Add(2, 3)), 4)) == "1 + (2 + 3) + 4"

    assert str(Div(Mul(1, 2), 3)) == "1 * 2 / 3"
    assert str(Mul(Div(1, 2), 3)) == "1 / 2 * 3"
    assert str(Mul(1, Div(2, 3))) == "1 * (2 / 3)"
    assert str(Mul(1, Mul(2, Mul(3, 4)))) == "1 * (2 * (3 * 4))"
    assert str(Mul(Mul(Mul(1, 2), 3), 4)) == "1 * 2 * 3 * 4"
    assert str(Mul(Mul(1, Mul(2, 3)), 4)) == "1 * (2 * 3) * 4"

    assert str(Mul(Add(1, 2), Add(3, 4))) == "(1 + 2) * (3 + 4)"


def test_cnode_expressions():
    A = Symbol("A")
    B = Symbol("B")
    i = Symbol("i")
    j = Symbol("j")
    x = Symbol("x")
    y = Symbol("y")

    # Literals
    assert str(LiteralInt(123)) == "123"
    assert str(LiteralFloat(0.0)) == "0.0"
    assert str(LiteralFloat(1.0)) == "1.0"
    assert str(LiteralFloat(12.3)) == "12.3"  # 1.23e+01"

    # Variables
    # TODO: VariableAccess

    # Arrays
    assert str(ArrayAccess("A", (1,))) == "A[1]"
    assert str(ArrayAccess(A, (1, 2))) == "A[1][2]"
    assert str(A[1, 2, 3]) == "A[1][2][3]"
    assert str(ArrayAccess(ArrayDecl("double", "A", (2,)), 1)) == "A[1]"
    assert str(ArrayDecl("double", A, (2, 3))[1, 2]) == "A[1][2]"

    # FlattenedArray
    n = Symbol("n")
    decl = ArrayDecl("double", A, (4,))
    assert str(FlattenedArray(decl, strides=(2,), offset=3)[0]) == "A[3]"  # "A[3 + 2 * 0]"
    assert str(FlattenedArray(decl, strides=(2,))[0]) == "A[0]"  # "A[2 * 0]"
    decl = ArrayDecl("double", A, (2, 3, 4))
    flattened = FlattenedArray(decl, strides=(7, 8 * n, n - 1))
    #assert str(flattened[0, n, n * 7]) == "A[7 * 0 + 8 * n * n + (n - 1) * (n * 7)]"
    #assert str(flattened[0, n][n * 7]) == "A[7 * 0 + 8 * n * n + (n - 1) * (n * 7)]"
    #assert str(flattened[0][n][n * 7]) == "A[7 * 0 + 8 * n * n + (n - 1) * (n * 7)]"
    assert str(flattened[0, n, n * 7]) == "A[8 * n * n + (n - 1) * (n * 7)]"
    assert str(flattened[0, n][n * 7]) == "A[8 * n * n + (n - 1) * (n * 7)]"
    assert str(flattened[0][n][n * 7]) == "A[8 * n * n + (n - 1) * (n * 7)]"

    # Unary operators
    assert str(Pos(1)) == "+1"
    assert str(Neg(1)) == "-1"
    assert str(Not(1)) == "!1"
    assert str(BitNot(1)) == "~1"
    assert str(PreIncrement(i)) == "++i"
    assert str(PostIncrement(i)) == "i++"
    assert str(PreDecrement(i)) == "--i"
    assert str(PostDecrement(i)) == "i--"
    assert str(Call(Symbol("f"), [])) == "f()"
    assert str(Call("f", [x, 3, Mul(y, 8.0)])) == "f(x, 3, y * 8.0)"
    assert str(Call("sin", Mul(B, 3.0))) == "sin(B * 3.0)"

    # Binary operators
    assert str(Add(1, 2)) == "1 + 2"
    assert str(Mul(1, 2)) == "1 * 2"
    assert str(Div(1, 2)) == "1 / 2"
    assert str(Sub(1, 2)) == "1 - 2"
    assert str(Mod(1, 2)) == "1 % 2"
    assert str(EQ(1, 2)) == "1 == 2"
    assert str(NE(1, 2)) == "1 != 2"
    assert str(LT(1, 2)) == "1 < 2"
    assert str(LE(1, 2)) == "1 <= 2"
    assert str(GT(1, 2)) == "1 > 2"
    assert str(GE(1, 2)) == "1 >= 2"
    assert str(And(1, 2)) == "1 && 2"
    assert str(Or(1, 2)) == "1 || 2"
    assert str(BitAnd(1, 2)) == "1 & 2"
    assert str(BitXor(1, 2)) == "1 ^ 2"
    assert str(BitOr(1, 2)) == "1 | 2"

    # Binary operators translated from python
    assert str(A + B) == "A + B"
    assert str(A * B) == "A * B"
    assert str(A / B) == "A / B"
    assert str(A - B) == "A - B"

    # Ternary operator
    assert str(Conditional(1, 2, 3)) == "1 ? 2 : 3"

    # N-ary "operators" simplify code generation
    assert str(Sum([1, 2, 3, 4])) == "1 + 2 + 3 + 4"
    assert str(Product([1, 2, 3, 4])) == "1 * 2 * 3 * 4"
    assert str(Product([Sum([1, 2, 3]), Sub(4, 5)])) == "(1 + 2 + 3) * (4 - 5)"

    # Custom expression
    assert str(Mul(VerbatimExpr("1 + std::foo(3)"), Add(5, 6))) == "(1 + std::foo(3)) * (5 + 6)"


def test_cnode_assignments():
    x = Symbol("x")
    y = Symbol("y")
    assert str(Assign(x, y)) == "x = y"
    assert str(AssignAdd(x, y)) == "x += y"
    assert str(AssignSub(x, y)) == "x -= y"
    assert str(AssignMul(x, y)) == "x *= y"
    assert str(AssignDiv(x, y)) == "x /= y"
    assert str(AssignMod(x, y)) == "x %= y"
    assert str(AssignLShift(x, y)) == "x <<= y"
    assert str(AssignRShift(x, y)) == "x >>= y"
    assert str(AssignAnd(x, y)) == "x &&= y"
    assert str(AssignOr(x, y)) == "x ||= y"
    assert str(AssignBitAnd(x, y)) == "x &= y"
    assert str(AssignBitXor(x, y)) == "x ^= y"
    assert str(AssignBitOr(x, y)) == "x |= y"


def test_cnode_variable_declarations():
    y = Symbol("y")
    assert str(VariableDecl("foo", "x")) == "foo x;"
    assert str(VariableDecl("int", "n", 1)) == "int n = 1;"
    assert str(VariableDecl("double", "x", Mul(y, 3.0))) == "double x = y * 3.0;"


def test_1d_initializer_list():
    fmt = lambda v, s: format_indented_lines(build_initializer_lists(v, s, 0, str))
    assert fmt([], (0,)) == "{  }"
    assert fmt([1], (1,)) == "{ 1 }"
    assert fmt([1, 2], (2,)) == "{ 1, 2 }"
    assert fmt([1, 2, 3], (3,)) == "{ 1, 2, 3 }"


def test_nd_initializer_list_oneitem():
    fmt = lambda v, s: format_indented_lines(build_initializer_lists(v, s, 0, str))
    assert fmt([1], (1,)) == "{ 1 }"
    assert fmt([[1]], (1, 1)) == "{ { 1 } }"
    assert fmt([[[1]]], (1, 1, 1)) == "{ { { 1 } } }"
    assert fmt([[[[1]]]], (1, 1, 1, 1)) == "{ { { { 1 } } } }"


def test_nd_initializer_list_twoitems():
    fmt = lambda v, s: format_indented_lines(build_initializer_lists(v, s, 0, str))
    assert fmt([1, 2], (2,)) == "{ 1, 2 }"
    assert fmt([[1, 2]], (1, 2)) == "{ { 1, 2 } }"
    assert fmt([[[1, 2]]], (1, 1, 2)) == "{ { { 1, 2 } } }"
    assert fmt([[[[1, 2]]]], (1, 1, 1, 2)) == "{ { { { 1, 2 } } } }"
    # transpose it:
    assert fmt([[1], [2]], (2, 1)) == "{ { 1 },\n  { 2 } }"
    assert fmt([[[1], [2]]], (1, 2, 1)) == "{ { { 1 },\n    { 2 } } }"


def test_nd_initializer_list_twobytwoitems():
    fmt = lambda v, s: format_indented_lines(build_initializer_lists(v, s, 0, str))
    assert fmt([[1, 2], [3, 4]], (2, 2)) == "{ { 1, 2 },\n  { 3, 4 } }"
    assert fmt([[[1, 2], [3, 4]]], (1, 2, 2)) == "{ { { 1, 2 },\n    { 3, 4 } } }"
    assert fmt([[[[1, 2], [3, 4]]]], (1, 1, 2, 2)) == "{ { { { 1, 2 },\n      { 3, 4 } } } }"


def test_2d_initializer_list():
    assert format_indented_lines(build_initializer_lists([[1, 2, 3], [4, 5, 6]], (2, 3), 0, str)) == "{ { 1, 2, 3 },\n  { 4, 5, 6 } }"

    values = [[[1], [2]], [[3], [4]], [[5], [6]]]
    reference = """\
{ { { 1 },
    { 2 } },
  { { 3 },
    { 4 } },
  { { 5 },
    { 6 } } }"""
    assert format_indented_lines(build_initializer_lists(values, (3, 2, 1), 0, str)) == reference


def test_2d_numpy_initializer_list():
    import numpy
    values = [[1, 2, 3], [4, 5, 6]]
    array = numpy.asarray(values)
    sh = (2, 3)
    assert array.shape == sh
    fmt = lambda v, s: format_indented_lines(build_initializer_lists(v, s, 0, str))
    assert fmt(values, sh) == "{ { 1, 2, 3 },\n  { 4, 5, 6 } }"
    assert fmt(array, sh) == "{ { 1, 2, 3 },\n  { 4, 5, 6 } }"

    values = [[[1], [2]], [[3], [4]], [[5], [6]]]
    array = numpy.asarray(values)
    sh = (3, 2, 1)
    assert sh == array.shape
    reference = """\
{ { { 1 },
    { 2 } },
  { { 3 },
    { 4 } },
  { { 5 },
    { 6 } } }"""
    assert fmt(values, sh) == reference
    assert fmt(array, sh) == reference


def test_cnode_array_declarations():
    assert str(ArrayDecl("double", "x", 3)) == "double x[3];"
    assert str(ArrayDecl("double", "x", (3,))) == "double x[3];"
    assert str(ArrayDecl("double", "x", (3, 4))) == "double x[3][4];"

    assert str(ArrayDecl("double", "x", 3, [1., 2., 3.])) == "double x[3] = { 1.0, 2.0, 3.0 };"
    assert str(ArrayDecl("double", "x", (3,), [1., 2., 3.])) == "double x[3] = { 1.0, 2.0, 3.0 };"
    reference = """\
{
    double x[2][3] =
        { { 1.0, 2.0, 3.0 },
          { 4.0, 5.0, 6.0 } };
}"""
    assert str(Scope(ArrayDecl("double", "x", (2, 3), [[1., 2., 3.], [4., 5., 6.]]))) == reference


def test_cnode_comments():
    assert str(Comment("hello world")) == "// hello world"
    assert str(Comment("  hello\n world  ")) == "// hello\n// world"
    assert format_indented_lines(Indented(Comment("  hello\n world  ").cs_format())) == "    // hello\n    // world"


def test_cnode_statements():
    assert str(Break()) == "break;"
    assert str(Continue()) == "continue;"
    assert str(Return(Add(1, 2))) == "return 1 + 2;"
    assert str(Case(3)) == "case 3:"
    assert str(Default()) == "default:"

    code = "for (std::vector<int>::iterator it = v.begin(); it != v.end(); ++it)\n{    /* foobar */\n}"
    assert str(VerbatimStatement(code)) == code


def test_cnode_loop_statements():
    i, j, x, y, A = [Symbol(ch) for ch in "ijxyA"]

    body = [Assign("x", 3), AssignAdd(x, 5)]
    body_fmt = "{\n    x = 3;\n    x += 5;\n}"
    assert str(Scope(body)) == body_fmt
    assert str(Namespace("foo", body)) == "namespace foo\n" + body_fmt
    assert str(If(LT(x, 4.0), body)) == "if (x < 4.0)\n" + body_fmt
    assert str(ElseIf(LT(x, 4.0), body)) == "else if (x < 4.0)\n" + body_fmt
    assert str(Else(body)) == "else\n" + body_fmt
    assert str(While(LT(x, 4.0), body)) == "while (x < 4.0)\n" + body_fmt
    assert str(Do(LT(x, 4.0), body)) == "do\n" + body_fmt + " while (x < 4.0);"
    assert str(For(VariableDecl("int", "i", 0), LT(i, 4), PreIncrement(i), body)) == "for (int i = 0; i < 4; ++i)\n" + body_fmt

    assert str(ForRange("i", 3, 7, Comment("body"))) == "for (int i = 3; i < 7; ++i)\n{\n    // body\n}"

    # Using assigns as both statements and expressions
    assert str(While(LT(AssignAdd("x", 4.0), 17.0), AssignAdd("A", y))) == "while ((x += 4.0) < 17.0)\n{\n    A += y;\n}"
    assert str(ForRange("i", 3, 7, AssignAdd("A", i))) == "for (int i = 3; i < 7; ++i)\n    A += i;"


def test_cnode_loop_helpers():
    i = Symbol("i")
    j = Symbol("j")
    A = Symbol("A")
    B = Symbol("B")
    C = Symbol("C")
    dst = A[i + 4 * j]
    src = 2.0 * B[j] * C[i]
    ranges = [(i, 0, 2), (j, 1, 3)]
    assert str(assign_loop(dst, src, ranges)) == """\
for (int i = 0; i < 2; ++i)
    for (int j = 1; j < 3; ++j)
        A[i + 4 * j] = 2.0 * B[j] * C[i];"""
    assert str(scale_loop(dst, src, ranges)) == """\
for (int i = 0; i < 2; ++i)
    for (int j = 1; j < 3; ++j)
        A[i + 4 * j] *= 2.0 * B[j] * C[i];"""
    assert str(accumulate_loop(dst, src, ranges)) == """\
for (int i = 0; i < 2; ++i)
    for (int j = 1; j < 3; ++j)
        A[i + 4 * j] += 2.0 * B[j] * C[i];"""


def test_cnode_switch_statements():
    i, j, x, y, A = [Symbol(ch) for ch in "ijxyA"]
    assert str(Switch("x", [])) == "switch (x)\n{\n}"
    assert str(Switch("x", [], default=Assign("i", 3))) == "switch (x)\n{\ndefault:\n    {\n        i = 3;\n    }\n}"

    reference_switch = """switch (x)
{
case 1:
    {
        y = 3;
    }
    break;
case 2:
    {
        y = 4;
    }
    break;
default:
    {
        y = 5;
    }
}"""
    cnode_switch = str(Switch("x",
                              [(1, Assign("y", 3)), (2, Assign("y", 4)), ],
                              default=Assign("y", 5)))
    assert cnode_switch == reference_switch

    reference_switch = """switch (x)
{
case 1:
    y = 3;
case 2:
    y = 4;
default:
    y = 5;
}"""
    cnode_switch = str(Switch("x",
                              [(1, Assign("y", 3)), (2, Assign("y", 4)), ],
                              default=Assign("y", 5),
                              autobreak=False, autoscope=False))
    assert cnode_switch == reference_switch


def test_conceptual_tabulate_tensor():
    A = ArrayDecl("double", "A", (4, 6), values=0.0)
    code = StatementList([
        A,
        ForRange("q", 0, 2, [
            VariableDecl("double", "x", ArrayAccess("quadpoints", "q")),
            ForRange("i", 0, 4, [
                ForRange("j", 0, 6, [
                    AssignAdd(ArrayAccess(A, ("i", "j")),
                              Mul(ArrayAccess("FE0", ("q", "i")),
                                  ArrayAccess("FE1", ("q", "j"))))
                ])
            ])
        ])
    ])
    print(str(code))
