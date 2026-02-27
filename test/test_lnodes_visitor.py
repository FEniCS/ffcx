"""Tests for LNode immutability, visitor, and transformer infrastructure."""

import pytest

from ffcx.codegeneration import lnodes as L
from ffcx.codegeneration.visitor import LNodeTransformer, LNodeVisitor

# --- Immutability tests ---


class TestImmutability:
    """Verify frozen dataclass nodes raise on attribute assignment."""

    def test_literal_float_frozen(self):
        f = L.LiteralFloat(3.14)
        with pytest.raises(AttributeError):
            f.value = 2.0

    def test_literal_int_frozen(self):
        i = L.LiteralInt(42)
        with pytest.raises(AttributeError):
            i.value = 0

    def test_symbol_frozen(self):
        s = L.Symbol("x", L.DataType.REAL)
        with pytest.raises(AttributeError):
            s.name = "y"

    def test_add_frozen(self):
        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        expr = L.Add(a, b)
        with pytest.raises(AttributeError):
            expr.lhs = b

    def test_product_frozen(self):
        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        p = L.Product([a, b])
        with pytest.raises(AttributeError):
            p.args = (b, a)

    def test_product_args_is_tuple(self):
        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        p = L.Product([a, b])
        assert isinstance(p.args, tuple)

    def test_neg_frozen(self):
        a = L.Symbol("a", L.DataType.REAL)
        n = L.Neg(a)
        with pytest.raises(AttributeError):
            n.arg = L.Symbol("b", L.DataType.REAL)

    def test_array_access_frozen(self):
        a = L.Symbol("a", L.DataType.REAL)
        aa = L.ArrayAccess(a, [L.LiteralInt(0)])
        with pytest.raises(AttributeError):
            aa.array = L.Symbol("b", L.DataType.REAL)

    def test_comment_frozen(self):
        c = L.Comment("hello")
        with pytest.raises(AttributeError):
            c.comment = "world"

    def test_statement_frozen(self):
        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        s = L.Statement(L.Assign(a, b))
        with pytest.raises(AttributeError):
            s.expr = L.Assign(b, a)

    def test_statement_list_frozen(self):
        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        sl = L.StatementList([L.Assign(a, b)])
        with pytest.raises(AttributeError):
            sl.statements = ()

    def test_statement_list_statements_is_tuple(self):
        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        sl = L.StatementList([L.Assign(a, b)])
        assert isinstance(sl.statements, tuple)

    def test_section_frozen(self):
        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        sec = L.Section("test", [L.Assign(a, b)], [])
        with pytest.raises(AttributeError):
            sec.name = "other"

    def test_section_statements_is_tuple(self):
        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        sec = L.Section("test", [L.Assign(a, b)], [])
        assert isinstance(sec.statements, tuple)

    def test_for_range_frozen(self):
        i = L.Symbol("i", L.DataType.INT)
        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        fr = L.ForRange(i, 0, 10, [L.Assign(a, b)])
        with pytest.raises(AttributeError):
            fr.index = L.Symbol("j", L.DataType.INT)

    def test_variable_decl_frozen(self):
        a = L.Symbol("a", L.DataType.REAL)
        vd = L.VariableDecl(a)
        with pytest.raises(AttributeError):
            vd.symbol = L.Symbol("b", L.DataType.REAL)

    def test_array_decl_frozen(self):
        a = L.Symbol("a", L.DataType.REAL)
        ad = L.ArrayDecl(a, sizes=(3,))
        with pytest.raises(AttributeError):
            ad.sizes = (4,)


# --- Equality tests ---


class TestEquality:
    """Verify custom __eq__ still works correctly."""

    def test_symbol_eq_by_name(self):
        s1 = L.Symbol("x", L.DataType.REAL)
        s2 = L.Symbol("x", L.DataType.INT)
        assert s1 == s2  # Symbol compares only name

    def test_symbol_neq(self):
        s1 = L.Symbol("x", L.DataType.REAL)
        s2 = L.Symbol("y", L.DataType.REAL)
        assert not (s1 == s2)

    def test_literal_float_eq(self):
        f1 = L.LiteralFloat(3.14)
        f2 = L.LiteralFloat(3.14)
        assert f1 == f2

    def test_literal_int_eq(self):
        i1 = L.LiteralInt(42)
        i2 = L.LiteralInt(42)
        assert i1 == i2

    def test_add_eq(self):
        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        assert L.Add(a, b) == L.Add(a, b)

    def test_add_neq(self):
        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        assert not (L.Add(a, b) == L.Sub(a, b))

    def test_product_eq(self):
        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        assert L.Product([a, b]) == L.Product([a, b])

    def test_array_decl_eq_fixed_bug(self):
        """Verify the ArrayDecl.__eq__ bug fix (was comparing self to self)."""
        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        ad1 = L.ArrayDecl(a, sizes=(3,))
        ad2 = L.ArrayDecl(b, sizes=(3,))
        assert not (ad1 == ad2)  # Different symbols should not be equal


# --- children() / named_children() tests ---


class TestChildren:
    """Verify correct child enumeration for each node type."""

    def test_literal_no_children(self):
        f = L.LiteralFloat(1.0)
        assert f.named_children() == []
        assert f.children() == []

    def test_symbol_no_children(self):
        s = L.Symbol("x", L.DataType.REAL)
        assert s.named_children() == []
        assert s.children() == []

    def test_neg_children(self):
        a = L.Symbol("a", L.DataType.REAL)
        n = L.Neg(a)
        assert len(n.named_children()) == 1
        assert n.named_children()[0] == ("arg", a)
        assert n.children() == [a]

    def test_add_children(self):
        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        expr = L.Add(a, b)
        named = expr.named_children()
        assert len(named) == 2
        assert named[0] == ("lhs", a)
        assert named[1] == ("rhs", b)

    def test_product_children(self):
        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        p = L.Product([a, b])
        named = p.named_children()
        assert len(named) == 1
        assert named[0][0] == "args"
        assert list(named[0][1]) == [a, b]

    def test_array_access_children(self):
        a = L.Symbol("a", L.DataType.REAL)
        idx = L.LiteralInt(0)
        aa = L.ArrayAccess(a, [idx])
        named = aa.named_children()
        assert len(named) == 2
        assert named[0] == ("array", a)
        assert named[1][0] == "indices"

    def test_conditional_children(self):
        cond = L.LiteralInt(1)
        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        c = L.Conditional(cond, a, b)
        named = c.named_children()
        assert len(named) == 3

    def test_section_children(self):
        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        sec = L.Section("test", [L.Assign(a, b)], [])
        named = sec.named_children()
        assert len(named) == 2
        assert named[0][0] == "statements"
        assert named[1][0] == "declarations"

    def test_for_range_children(self):
        i = L.Symbol("i", L.DataType.INT)
        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        fr = L.ForRange(i, 0, 10, [L.Assign(a, b)])
        named = fr.named_children()
        assert len(named) == 4
        names = [n[0] for n in named]
        assert names == ["index", "begin", "end", "body"]

    def test_comment_no_children(self):
        c = L.Comment("hello")
        assert c.named_children() == []
        assert c.children() == []


# --- replace() tests ---


class TestReplace:
    """Verify functional node replacement."""

    def test_symbol_replace(self):
        s = L.Symbol("x", L.DataType.REAL)
        s2 = s.replace(name="y")
        assert s2.name == "y"
        assert s.name == "x"  # Original unchanged

    def test_add_replace_lhs(self):
        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        c = L.Symbol("c", L.DataType.REAL)
        expr = L.Add(a, b)
        expr2 = expr.replace(lhs=c)
        assert expr2.lhs == c
        assert expr2.rhs == b
        assert expr.lhs == a  # Original unchanged

    def test_product_replace_args(self):
        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        c = L.Symbol("c", L.DataType.REAL)
        p = L.Product([a, b])
        p2 = p.replace(args=[a, b, c])
        assert len(p2.args) == 3
        assert len(p.args) == 2  # Original unchanged

    def test_for_range_replace_body(self):
        i = L.Symbol("i", L.DataType.INT)
        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        fr = L.ForRange(i, 0, 10, [L.Assign(a, b)])
        fr2 = fr.replace(end=L.LiteralInt(20))
        assert fr2.end == L.LiteralInt(20)
        assert fr.end == L.LiteralInt(10)  # Original unchanged


# --- Visitor tests ---


class TestVisitor:
    """Test LNodeVisitor."""

    def test_visitor_collects_symbols(self):
        """Visitor that collects all Symbol names in a tree."""

        class SymbolCollector(LNodeVisitor):
            def __init__(self):
                self.names = []

            def visit_Symbol(self, node):
                self.names.append(node.name)

        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        expr = L.Add(a, L.Mul(b, a))

        collector = SymbolCollector()
        collector.visit(expr)
        assert sorted(collector.names) == ["a", "a", "b"]


# --- Transformer tests ---


class TestTransformer:
    """Test LNodeTransformer."""

    def test_identity_transform(self):
        """Default transformer returns same object (no unnecessary copies)."""
        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        expr = L.Add(a, b)

        t = LNodeTransformer()
        result = t.visit(expr)
        assert result is expr

    def test_symbol_renamer(self):
        """Transformer that renames symbol 'a' to 'z'."""

        class SymbolRenamer(LNodeTransformer):
            def visit_Symbol(self, node):
                if node.name == "a":
                    return L.Symbol("z", node.dtype)
                return node

        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        expr = L.Add(a, b)

        renamer = SymbolRenamer()
        result = renamer.visit(expr)

        assert isinstance(result, L.Add)
        assert result.lhs == L.Symbol("z", L.DataType.REAL)
        assert result.rhs == b
        # Original unchanged
        assert expr.lhs == a

    def test_transform_nested(self):
        """Transformer works on nested expressions."""

        class DoubleConstants(LNodeTransformer):
            def visit_LiteralFloat(self, node):
                return L.LiteralFloat(node.value * 2)

        a = L.Symbol("a", L.DataType.REAL)
        expr = L.Add(a, L.LiteralFloat(3.0))

        result = DoubleConstants().visit(expr)
        assert isinstance(result, L.Add)
        assert result.lhs is a  # Symbol unchanged
        assert result.rhs == L.LiteralFloat(6.0)

    def test_identity_on_product(self):
        """Identity transform on Product returns same object."""
        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        p = L.Product([a, b])

        t = LNodeTransformer()
        result = t.visit(p)
        assert result is p

    def test_transform_product_args(self):
        """Transformer that renames symbols in a Product."""

        class SymbolRenamer(LNodeTransformer):
            def visit_Symbol(self, node):
                if node.name == "a":
                    return L.Symbol("z", node.dtype)
                return node

        a = L.Symbol("a", L.DataType.REAL)
        b = L.Symbol("b", L.DataType.REAL)
        p = L.Product([a, b])

        result = SymbolRenamer().visit(p)
        assert isinstance(result, L.Product)
        assert result.args[0] == L.Symbol("z", L.DataType.REAL)
        assert result.args[1] == b
