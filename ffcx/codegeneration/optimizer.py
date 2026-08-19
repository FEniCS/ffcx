"""Optimizer."""

import math
from collections import defaultdict

import ffcx.codegeneration.lnodes as L
from ffcx.ir.representationutils import QuadratureRule


def _pow2_exponent(node: L.LExpr) -> int | None:
    """Return e such that node == 2**e, if node is an exact literal power of two."""
    if isinstance(node, (L.LiteralFloat, L.LiteralInt)):
        if isinstance(node.value, complex):
            return None
        value = float(node.value)
        if value <= 0.0:
            return None
        mantissa, exponent = math.frexp(value)
        if mantissa == 0.5:
            return exponent - 1
    return None


def optimize(code: list[L.LNode], quadrature_rule: QuadratureRule) -> list[L.LNode]:
    """Optimize code.

    Args:
        code: List of LNodes to optimize.
        quadrature_rule: TODO.

    Returns:
        Optimized list of LNodes.
    """
    # Fuse sections with the same name and same annotations
    code = fuse_sections(code, "Coefficient")
    code = fuse_sections(code, "Jacobian")
    for i, section in enumerate(code):
        if isinstance(section, L.Section):
            if section.name == "Jacobian":
                section = reciprocal_cse(section)
            if L.Annotation.fuse in section.annotations:
                section = fuse_loops(section)
            if L.Annotation.licm in section.annotations:
                section = licm(section, quadrature_rule)
            code[i] = section

    return code


def reciprocal_cse_statements(
    statements: list[L.LNode], name_prefix: str = "recip"
) -> list[L.LNode]:
    """Replace repeated divisions by the same divisor with one reciprocal and multiplies.

    An affine cell's pseudo-inverse Jacobian divides every cofactor by the same
    determinant. FFCx's existing CSE already recognises the shared divisor but
    still emits one hardware division per entry; a reciprocal computed once and
    reused as a multiply is exact for the same divide (same IEEE rounding as any
    other division) and division is markedly more expensive than multiplication
    on typical hardware.

    Args:
        statements: A flat list of statements (no nested loops expected --
            this targets piecewise/varying scalar partitions, not quadrature
            loop bodies).
        name_prefix: Prefix for the generated reciprocal symbols, kept distinct
            per call site so two calls in the same kernel can't clash.

    Returns:
        A new statement list with repeated divisions rewritten.
    """
    # Keyed by the divisor's structural identity (see `_key`), not the LExpr
    # object itself: a literal divisor (e.g. dividing by a plain constant)
    # is a `LiteralFloat`, which -- unlike `LiteralInt` -- has no `__hash__`,
    # so using it directly as a dict key crashes.
    divisor_counts: dict[tuple, int] = defaultdict(int)
    for stmt in statements:
        if isinstance(stmt, L.VariableDecl) and isinstance(stmt.value, L.Div):
            divisor_counts[_key(stmt.value.rhs)] += 1

    recip_symbols: dict[tuple, L.Symbol] = {}
    counter = 0
    new_statements: list[L.LNode] = []
    for stmt in statements:
        if (
            isinstance(stmt, L.VariableDecl)
            and isinstance(stmt.value, L.Div)
            and divisor_counts[_key(stmt.value.rhs)] >= 2
        ):
            divisor = stmt.value.rhs
            divisor_key = _key(divisor)
            if divisor_key not in recip_symbols:
                recip_symbol = L.Symbol(f"{name_prefix}_{counter}", L.DataType.SCALAR)
                counter += 1
                recip_symbols[divisor_key] = recip_symbol
                new_statements.append(
                    L.VariableDecl(recip_symbol, L.Div(L.LiteralFloat(1.0), divisor))
                )
            recip_symbol = recip_symbols[divisor_key]
            new_statements.append(L.VariableDecl(stmt.symbol, L.Mul(stmt.value.lhs, recip_symbol)))
        else:
            new_statements.append(stmt)

    return new_statements


def reciprocal_cse(section: L.Section) -> L.Section:
    """Apply :func:`reciprocal_cse_statements` to a Section's statements in place."""
    section.statements = reciprocal_cse_statements(section.statements)
    return section


def _referenced_symbol_names(expr: L.LExpr, out: set[str]) -> None:
    """Collect the names of every Symbol referenced within an expression tree."""
    if isinstance(expr, L.Symbol):
        out.add(expr.name)
    elif isinstance(expr, L.ArrayAccess):
        out.add(expr.array.name)
        for i in expr.indices:
            _referenced_symbol_names(i, out)
    elif isinstance(expr, L.Conditional):
        _referenced_symbol_names(expr.condition, out)
        _referenced_symbol_names(expr.true, out)
        _referenced_symbol_names(expr.false, out)
    elif isinstance(expr, L.BinOp):
        _referenced_symbol_names(expr.lhs, out)
        _referenced_symbol_names(expr.rhs, out)
    elif isinstance(expr, L.PrefixUnaryOp):
        _referenced_symbol_names(expr.arg, out)
    elif isinstance(expr, (L.NaryOp, L.MathFunction)):
        for a in expr.args:
            _referenced_symbol_names(a, out)


def power_of_two_cse_statements(
    statements: list[L.LNode], speculative_out: set[str] | None = None
) -> list[L.LNode]:
    """Fold scalar chains that round-trip through exact powers of two.

    UFL's expansion of e.g. a symmetric gradient introduces factors of 2 and
    1/2 (Voigt-style off-diagonal terms, a doubled derivative combined with a
    later halving, etc.) that don't originate from the same place and so
    aren't recognised as equal by ordinary CSE, even though multiplying and
    dividing by an exact power of two never rounds in IEEE-754: ``(x + x) *
    0.5`` is bit-identical to ``x`` for any finite, non-subnormal x. This walks
    a flat statement list tracking, for each declared scalar, the (base
    symbol, power-of-two exponent) pair its value is exactly equal to, and
    replaces a statement whose computed value is already known under a
    different name with a plain reference to that name instead of
    recomputing it.

    Rerouting a statement straight back to `base` (or to an earlier scaled
    symbol) can leave the statement it used to depend on with no remaining
    reader -- if nothing else ever needed that same intermediate value, it
    becomes dead code, which GCC rejects under ``-Werror=unused-variable``.
    Whether that's actually true can depend on code well outside this one
    statement list (another partition's statements reading this one's
    output, or the tensor contraction reading it directly), so this
    function does not decide the question itself: every symbol it keeps
    only speculatively (in case a later statement in *this* list needs the
    same value) is recorded into `speculative_out` for the caller to check
    against the fully assembled kernel body, via :func:`prune_dead_scalars`.

    Args:
        statements: A flat list of statements (piecewise/varying partition,
            no nested loops).
        speculative_out: If given, updated with the names of symbols kept
            only speculatively by this call.

    Returns:
        A new statement list with round-trip duplicates folded to aliases.
    """
    if speculative_out is None:
        speculative_out = set()
    # symbol name -> (base symbol, exponent) with value == base * 2**exponent
    canonical: dict[str, tuple[L.Symbol, int]] = {}
    # (base symbol name, exponent) -> symbol already known to hold that value
    seen: dict[tuple[str, int], L.Symbol] = {}

    def resolve(node: L.LExpr) -> tuple[L.Symbol, int] | None:
        """Look up a symbol's tracked (base, exponent), if any."""
        if isinstance(node, L.Symbol) and node.name in canonical:
            return canonical[node.name]
        return None

    new_statements: list[L.LNode] = []
    for stmt in statements:
        if not isinstance(stmt, L.VariableDecl):
            new_statements.append(stmt)
            continue

        value = stmt.value
        result: tuple[L.Symbol, int] | None = None

        if isinstance(value, L.Add) and value.lhs == value.rhs:
            operand = resolve(value.lhs)
            if operand is not None:
                base, exponent = operand
                result = (base, exponent + 1)
        elif isinstance(value, L.Mul):
            for factor, other in ((value.lhs, value.rhs), (value.rhs, value.lhs)):
                exponent_shift = _pow2_exponent(other)
                operand = resolve(factor)
                if exponent_shift is not None and operand is not None:
                    base, exponent = operand
                    result = (base, exponent + exponent_shift)
                    break
        elif isinstance(value, L.Div):
            exponent_shift = _pow2_exponent(value.rhs)
            operand = resolve(value.lhs)
            if exponent_shift is not None and operand is not None:
                base, exponent = operand
                result = (base, exponent - exponent_shift)

        if result is not None:
            base, exponent = result
            canonical[stmt.symbol.name] = (base, exponent)
            key = (base.name, exponent)
            if exponent == 0:
                # Exactly equal to `base` itself: a pure alias.
                new_statements.append(L.VariableDecl(stmt.symbol, base))
                continue
            elif key in seen:
                # Exactly equal to some earlier scaled symbol: also an alias.
                new_statements.append(L.VariableDecl(stmt.symbol, seen[key]))
                continue
            else:
                seen[key] = stmt.symbol
                speculative_out.add(stmt.symbol.name)
        else:
            # Opaque: this symbol is its own base, unrelated to anything earlier.
            canonical[stmt.symbol.name] = (stmt.symbol, 0)
            seen.setdefault((stmt.symbol.name, 0), stmt.symbol)

        new_statements.append(stmt)

    return new_statements


def _collect_referenced_names(node, out: set[str]) -> None:
    """Recursively collect every symbol name read anywhere in a code tree."""
    if isinstance(node, list):
        for n in node:
            _collect_referenced_names(n, out)
    elif isinstance(node, L.VariableDecl):
        if node.value is not None:
            _referenced_symbol_names(node.value, out)
    elif isinstance(node, L.ArrayDecl):
        pass
    elif isinstance(node, (L.Section, L.StatementList)):
        _collect_referenced_names(node.statements, out)
    elif isinstance(node, L.ForRange):
        _referenced_symbol_names(node.begin, out)
        _referenced_symbol_names(node.end, out)
        _collect_referenced_names(node.body, out)
    elif isinstance(node, L.Comment):
        pass
    elif isinstance(node, L.Statement):
        _referenced_symbol_names(node.expr, out)


def _drop_unreferenced(node, dead: set[str]):
    """Recursively rebuild a code tree, dropping VariableDecls named in `dead`."""
    if isinstance(node, list):
        return [
            _drop_unreferenced(n, dead)
            for n in node
            if not (isinstance(n, L.VariableDecl) and n.symbol.name in dead)
        ]
    if isinstance(node, L.Section):
        node.statements = _drop_unreferenced(node.statements, dead)
    elif isinstance(node, L.StatementList):
        node.statements = _drop_unreferenced(node.statements, dead)
    elif isinstance(node, L.ForRange):
        # Mutate the existing StatementList's statements in place rather
        # than wrapping the result in a new one: `body` may itself be a
        # single-element list holding one multi-statement StatementList,
        # and re-wrapping that a second time is a shape `as_statement`
        # cannot unwrap again.
        node.body.statements = _drop_unreferenced(node.body.statements, dead)
    return node


def prune_dead_scalars(code: list[L.LNode], candidates: set[str]) -> list[L.LNode]:
    """Remove scalar declarations in `candidates` that end up with no reader.

    :func:`power_of_two_cse_statements` can reroute a statement past an
    earlier one it used to depend on (e.g. straight back to the value it
    was doubled from), leaving that earlier statement with no remaining
    use -- dead code that GCC rejects under ``-Werror=unused-variable``.
    Call this once the whole kernel body is assembled (definitions,
    intermediates, tensor contraction and all), so that a declaration used
    only by code outside its own partition's own statement list (a
    different quadrature rule's varying partition consuming a piecewise
    value, say, or the tensor contraction reading it directly) is
    correctly recognised as live.

    Args:
        code: The complete generated kernel body.
        candidates: Names that are safe to drop if nothing reads them --
            i.e. ones `power_of_two_cse_statements` itself introduced
            speculatively, never anything the caller didn't flag as such.

    Returns:
        `code`, with any now-dead candidate declarations removed. `code`'s
        Sections/ForRanges are mutated in place; the returned list is a new
        top-level list.
    """
    candidates = set(candidates)
    while True:
        referenced: set[str] = set()
        _collect_referenced_names(code, referenced)
        dead = candidates - referenced
        if not dead:
            return code
        code = _drop_unreferenced(code, dead)
        candidates -= dead


def fuse_sections(code: list[L.LNode], name: str) -> list[L.LNode]:
    """Fuse sections with the same name.

    Args:
        code: List of LNodes to fuse.
        name: Common name used by the sections that should be fused

    Returns:
        Fused list of LNodes.
    """
    statements: list[L.LNode] = []
    indices: list[int] = []
    input: list[L.Symbol] = []
    output: list[L.Symbol] = []
    declarations: list[L.Declaration] = []
    annotations: list[L.Annotation] = []

    for i, section in enumerate(code):
        if isinstance(section, L.Section):
            if section.name == name:
                declarations.extend(section.declarations)
                statements.extend(section.statements)
                indices.append(i)
                input.extend(section.input)
                output.extend(section.output)
                annotations = section.annotations

    # Remove duplicated inputs
    input = list(set(input))
    # Remove duplicated outputs
    output = list(set(output))

    section = L.Section(name, statements, declarations, input, output, annotations)

    # Replace the first section with the fused section
    code = code.copy()
    if indices:
        code[indices[0]] = section
        # Remove the other sections
        code = [c for i, c in enumerate(code) if i not in indices[1:]]

    return code


def fuse_loops(code: L.Section) -> L.Section:
    """Fuse loops with the same range and same annotations.

    Args:
        code: List of LNodes to fuse.

    Returns:
        Fused list of LNodes.
    """
    loops = defaultdict(list)
    output_code = []
    for statement in code.statements:
        if isinstance(statement, L.ForRange):
            id = (statement.index, statement.begin, statement.end)
            loops[id].append(statement.body)
        else:
            output_code.append(statement)

    for range, body in loops.items():
        output_code.append(L.ForRange(*range, body))

    return L.Section(code.name, output_code, code.declarations, code.input, code.output)


def get_statements(statement: L.Statement | L.StatementList) -> list[L.LNode]:
    """Get statements from a statement list.

    Args:
        statement: Statement list.

    Returns:
        List of statements.
    """
    if isinstance(statement, L.StatementList):
        return [statement.expr for statement in statement.statements]
    else:
        return [statement.expr]


def check_dependency(statement: L.Statement, index: L.Symbol) -> bool:
    """Check if a statement depends on a given index.

    Args:
        statement: Statement to check.
        index: Index to check.

    Returns:
        True if statement depends on index, False otherwise.
    """
    if isinstance(statement, L.ArrayAccess):
        if index in statement.indices:
            return True
        else:
            for i in statement.indices:
                if isinstance(i, L.Sum) or isinstance(i, L.Product):
                    if index in i.args:
                        return True
    elif isinstance(statement, L.Symbol):
        return False
    elif isinstance(statement, L.LiteralFloat) or isinstance(statement, L.LiteralInt):
        return False
    else:
        raise NotImplementedError(f"Statement {statement} not supported.")

    return False


def _key(expr) -> tuple:
    """Hashable identity of an expression, used to group equal factors.

    LNodes expressions are not generally hashable, so build a structural key
    that compares equal exactly when two factors are the same expression.

    Args:
        expr: Expression appearing as a factor in a product.

    Returns:
        Hashable structural key.
    """
    if isinstance(expr, L.ArrayAccess):
        return ("array", expr.array.name, *(_key(i) for i in expr.indices))
    elif isinstance(expr, L.Symbol):
        return ("symbol", expr.name)
    elif isinstance(expr, L.LiteralFloat | L.LiteralInt):
        return ("literal", expr.value)
    args = getattr(expr, "args", None)
    if args is not None:
        return (type(expr).__name__, *(_key(a) for a in args))
    return (type(expr).__name__, repr(expr))


def licm(section: L.Section, quadrature_rule: QuadratureRule) -> L.Section:
    """Perform loop invariant code motion.

    Args:
        section: List of LNodes to optimize.
        quadrature_rule: TODO.

    Returns:
        Optimized list of LNodes.
    """
    assert L.Annotation.licm in section.annotations

    counter = 0

    # Check depth of loops
    depth = L.depth(section.statements[0])
    if depth != 2:
        return section

    # Get statements in the inner loop
    outer_loop = section.statements[0]
    inner_loop = outer_loop.body.statements[0]

    # Collect all expressions in the inner loop by corresponding RHS
    expressions = defaultdict(list)
    for body in inner_loop.body.statements:
        statements = get_statements(body)
        assert isinstance(statements, list)
        for statement in statements:
            assert isinstance(statement, L.AssignAdd)  # Expecting AssignAdd
            rhs = statement.rhs
            assert isinstance(rhs, L.Product)  # Expecting Sum
            lhs = statement.lhs
            assert isinstance(lhs, L.ArrayAccess)  # Expecting ArrayAccess
            expressions[lhs].append(rhs)

    pre_loop: list[L.LNode] = []
    inner_body: list[L.LNode] = []
    size = outer_loop.end.value - outer_loop.begin.value
    for lhs, rhs in expressions.items():
        # Group the terms of this entry by the factors that stay in the inner
        # loop. Terms sharing them differ only in what is hoisted, so one
        # temporary can hold their sum and the entry needs a single update.
        groups: dict[tuple, tuple[list, list]] = {}
        terms: list[L.LExpr] = []
        for r in rhs:
            hoist_candidates = []
            keep = []
            for arg in r.args:
                if check_dependency(arg, inner_loop.index):
                    keep.append(arg)
                else:
                    hoist_candidates.append(arg)
            if len(hoist_candidates) > 1:
                key = tuple(_key(k) for k in keep)
                groups.setdefault(key, (keep, []))[1].append(L.Product(hoist_candidates))
            else:
                # Nothing worth hoisting, keep the term unchanged
                terms.append(r)

        for keep, hoisted in groups.values():
            # Create new temp holding the sum of the hoisted terms
            name = f"temp_{counter}"
            counter += 1
            temp = L.Symbol(name, L.DataType.SCALAR)
            pre_loop.append(L.ArrayDecl(temp, size, [0]))
            rhs_hoisted = hoisted[0] if len(hoisted) == 1 else L.Sum(hoisted)
            body = L.Assign(L.ArrayAccess(temp, [outer_loop.index]), rhs_hoisted)
            pre_loop.append(L.ForRange(outer_loop.index, outer_loop.begin, outer_loop.end, [body]))
            access = L.ArrayAccess(temp, [outer_loop.index])
            terms.append(L.Product([*keep, access]))

        # One read-modify-write of the entry instead of one per term
        if terms:
            inner_body.append(L.AssignAdd(lhs, terms[0] if len(terms) == 1 else L.Sum(terms)))

    # Rebuild the loop nest around the regrouped body
    new_inner = L.ForRange(inner_loop.index, inner_loop.begin, inner_loop.end, inner_body)
    new_outer = L.ForRange(outer_loop.index, outer_loop.begin, outer_loop.end, [new_inner])
    section.statements = pre_loop + [new_outer] + section.statements[1:]

    return section
