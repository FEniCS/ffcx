"""Optimizer."""

import math
from collections import Counter, defaultdict

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

    N divisions by the same value cost more than one reciprocal plus N
    multiplications, since division is markedly more expensive than
    multiplication on typical hardware. This is the same reassociation
    `-ffast-math`/`-freciprocal-math` would apply automatically, done
    unconditionally here: `a * (1.0 / b)` is not in general bit-identical to
    `a / b` (two roundings instead of one) unless `b` is an exact power of
    two, but the difference is at the ULP level, negligible next to
    quadrature/discretisation error.

    Args:
        statements: A flat list of statements (no nested loops expected --
            this targets piecewise/varying scalar partitions, not quadrature
            loop bodies).
        name_prefix: Prefix for the generated reciprocal symbols, kept distinct
            per call site so two calls in the same kernel can't clash.

    Returns:
        A new statement list with repeated divisions rewritten.
    """
    # Keyed by the divisor's structural identity (`expr_key`), not the LExpr
    # itself: a literal divisor is a `LiteralFloat`, which has no `__hash__`.
    divisor_counts: dict[tuple, int] = defaultdict(int)
    for stmt in statements:
        if isinstance(stmt, L.VariableDecl) and isinstance(stmt.value, L.Div):
            divisor_counts[expr_key(stmt.value.rhs)] += 1

    recip_symbols: dict[tuple, L.Symbol] = {}
    counter = 0
    new_statements: list[L.LNode] = []
    for stmt in statements:
        if (
            isinstance(stmt, L.VariableDecl)
            and isinstance(stmt.value, L.Div)
            and divisor_counts[expr_key(stmt.value.rhs)] >= 2
        ):
            divisor = stmt.value.rhs
            divisor_key = expr_key(divisor)
            if divisor_key not in recip_symbols:
                recip_expr = L.Div(L.LiteralFloat(1.0), divisor)
                recip_symbol = L.Symbol(f"{name_prefix}_{counter}", recip_expr.dtype)
                counter += 1
                recip_symbols[divisor_key] = recip_symbol
                new_statements.append(L.VariableDecl(recip_symbol, recip_expr))
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
    1/2 from different places (Voigt-style off-diagonal terms, a doubled
    derivative later halved, etc.), so ordinary CSE doesn't recognise them
    as equal -- even though multiplying/dividing by an exact power of two
    never rounds in IEEE-754: ``(x + x) * 0.5`` is bit-identical to ``x``.
    This walks a flat statement list tracking, for each declared scalar, the
    (base symbol, power-of-two exponent) pair its value equals, and replaces
    a statement whose value is already known under a different name with a
    reference to that name.

    Rerouting a statement back to `base` (or an earlier scaled symbol) can
    leave the statement it depended on with no remaining reader -- dead
    code, which GCC rejects under ``-Werror=unused-variable``. Whether
    that's true can depend on code outside this statement list, so this
    function doesn't decide it: symbols kept only speculatively are recorded
    into `speculative_out` for the caller to check against the fully
    assembled kernel body, via :func:`prune_dead_scalars`.

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
    elif isinstance(node, L.Section):
        for declaration in node.declarations:
            if isinstance(declaration, L.VariableDecl) and declaration.value is not None:
                _referenced_symbol_names(declaration.value, out)
        _collect_referenced_names(node.statements, out)
    elif isinstance(node, L.StatementList):
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
        node.declarations = [
            declaration
            for declaration in node.declarations
            if not (isinstance(declaration, L.VariableDecl) and declaration.symbol.name in dead)
        ]
        node.output = [symbol for symbol in node.output if symbol.name not in dead]
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


def scalar_declaration_names(code: list[L.LNode]) -> set[str]:
    """Return the names of every scalar declaration in a code tree.

    Generated scalar declarations are pure computations, so an unread one can
    be removed safely. Keeping this analysis on LNodes, rather than parsing
    formatted C, preserves lexical structure and works for every backend.
    """
    names: set[str] = set()

    def visit(node) -> None:
        if isinstance(node, list):
            for child in node:
                visit(child)
        elif isinstance(node, L.VariableDecl):
            names.add(node.symbol.name)
        elif isinstance(node, L.Section):
            visit(node.declarations)
            visit(node.statements)
        elif isinstance(node, L.StatementList):
            visit(node.statements)
        elif isinstance(node, L.ForRange):
            visit(node.body)

    visit(code)
    return names


def referenced_symbol_names(code: list[L.LNode]) -> set[str]:
    """Return the symbols read in a code tree."""
    referenced: set[str] = set()
    _collect_referenced_names(code, referenced)
    return referenced


def prune_unused_tables(code: list[L.LNode]) -> list[L.LNode]:
    """Remove top-level constant tables without any remaining accesses."""
    referenced = referenced_symbol_names(code)
    return [
        node
        for node in code
        if not (isinstance(node, L.ArrayDecl) and node.const and node.symbol.name not in referenced)
    ]


def prune_dead_scalars(code: list[L.LNode], candidates: set[str]) -> list[L.LNode]:
    """Remove scalar declarations in `candidates` that end up with no reader.

    :func:`power_of_two_cse_statements` can reroute a statement past an
    earlier one it depended on, leaving that one dead -- code GCC rejects
    under ``-Werror=unused-variable``. Call this once the whole kernel body
    is assembled, so a declaration read only by code outside its own
    partition (e.g. the tensor contraction) is still recognised as live.

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


def expr_key(expr) -> tuple:
    """Hashable identity of an expression, used to recognise equal factors.

    LNodes expressions are not generally hashable, so build a structural key
    that compares equal exactly when two expressions are the same expression.

    Args:
        expr: Expression to key.

    Returns:
        Hashable structural key.
    """
    if isinstance(expr, L.ArrayAccess):
        return ("array", expr.array.name, *(expr_key(i) for i in expr.indices))
    elif isinstance(expr, L.Symbol):
        return ("symbol", expr.name)
    elif isinstance(expr, L.LiteralFloat | L.LiteralInt):
        return ("literal", expr.value)
    elif isinstance(expr, L.PrefixUnaryOp):
        # PrefixUnaryOp (Neg, Not) has no __repr__ override, so it would
        # otherwise fall through to the identity-based default below.
        return (type(expr).__name__, expr_key(expr.arg))
    elif isinstance(expr, L.MathFunction):
        return ("math", expr.function, *(expr_key(arg) for arg in expr.args))
    args = getattr(expr, "args", None)
    if args is not None:
        return (type(expr).__name__, *(expr_key(a) for a in args))
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
    # Every hoisted temp's assignment is collected here and emitted as one
    # shared loop below, rather than one ForRange per temp -- many small
    # loops back to back give the vectoriser nothing to work with.
    hoist_loop_body: list[L.LNode] = []
    # Structural key of an already-hoisted expression -> its temp symbol, so
    # identical hoisted expressions (e.g. vector Poisson's diagonal blocks)
    # are computed once instead of once per block entry.
    seen_hoisted: dict[tuple, L.Symbol] = {}
    for lhs, rhs in expressions.items():
        # Group terms by the factors that stay in the inner loop: terms
        # sharing them differ only in what's hoisted, so one temporary can
        # hold their sum and the entry needs a single update.
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
                key = tuple(expr_key(k) for k in keep)
                groups.setdefault(key, (keep, []))[1].append(L.Product(hoist_candidates))
            else:
                # Nothing worth hoisting, keep the term unchanged
                terms.append(r)

        for keep, hoisted in groups.values():
            rhs_hoisted = hoisted[0] if len(hoisted) == 1 else L.Sum(hoisted)
            hoisted_key = expr_key(rhs_hoisted)
            temp = seen_hoisted.get(hoisted_key)
            if temp is None:
                # Create new temp holding the sum of the hoisted terms
                name = f"temp_{counter}"
                counter += 1
                temp = L.Symbol(name, L.DataType.SCALAR)
                # No initialiser: every entry is written unconditionally by
                # the loop just below (a plain Assign, not AssignAdd), so a
                # zero-fill here would be a pure dead store.
                pre_loop.append(L.ArrayDecl(temp, size))
                hoist_loop_body.append(
                    L.Assign(L.ArrayAccess(temp, [outer_loop.index]), rhs_hoisted)
                )
                seen_hoisted[hoisted_key] = temp
            access = L.ArrayAccess(temp, [outer_loop.index])
            terms.append(L.Product([*keep, access]))

        # One read-modify-write of the entry instead of one per term
        if terms:
            inner_body.append(L.AssignAdd(lhs, terms[0] if len(terms) == 1 else L.Sum(terms)))

    if hoist_loop_body:
        pre_loop.append(
            L.ForRange(outer_loop.index, outer_loop.begin, outer_loop.end, hoist_loop_body)
        )

    # Rebuild the loop nest around the regrouped body
    new_inner = L.ForRange(inner_loop.index, inner_loop.begin, inner_loop.end, inner_body)
    new_outer = L.ForRange(outer_loop.index, outer_loop.begin, outer_loop.end, [new_inner])
    section.statements = pre_loop + [new_outer] + section.statements[1:]

    return section


def factor_scalar_sums(statements: list[L.LNode]) -> list[L.LNode]:
    """Factor bounded scalar polynomials and reuse equivalent expressions.

    Expand single-use arithmetic temporaries and short products, then extract
    the most frequent factor from each sum. Shared sums stay as symbols, and
    expansion stops at 32 terms to bound work on large differentiated forms.
    Prefer factors used elsewhere when frequencies tie, exposing shared
    subexpressions to the native compiler. Sorted monomials make this
    independent of Python's hash seed.

    Like reciprocal CSE, this reassociates floating-point arithmetic and can
    change rounding. Calls, array accesses and nonliteral divisors are opaque.
    The caller must prune scalar declarations after assembling the whole body.
    """
    max_terms = 32
    max_factors = 8
    uses: Counter[str] = Counter()

    def count(expr):
        if isinstance(expr, L.Symbol):
            uses[expr.name] += 1
        elif isinstance(expr, L.BinOp):
            count(expr.lhs)
            count(expr.rhs)
        elif isinstance(expr, L.NaryOp | L.MathFunction):
            for arg in expr.args:
                count(arg)
        elif isinstance(expr, L.PrefixUnaryOp):
            count(expr.arg)

    for stmt in statements:
        if isinstance(stmt, L.VariableDecl):
            count(stmt.value)

    # A polynomial maps a sorted tuple of atomic factor keys to a coefficient.
    expanded: dict[str, dict[tuple, complex | float | int]] = {}
    atoms: dict[tuple, L.LExpr] = {}
    output: list[L.LNode] = []
    seen: dict[tuple, L.Symbol] = {}

    def atom(expr):
        key = expr_key(expr)
        atoms[key] = expr
        return {(key,): 1}

    def polynomial(expr):
        if isinstance(expr, L.Symbol):
            return expanded[expr.name] if expr.name in expanded else atom(expr)
        if isinstance(expr, L.LiteralInt | L.LiteralFloat):
            return {(): expr.value}
        if isinstance(expr, L.Neg):
            return {key: -value for key, value in polynomial(expr.arg).items()}
        if isinstance(expr, L.Add):
            result = polynomial(expr.lhs).copy()
            for key, value in polynomial(expr.rhs).items():
                result[key] = result.get(key, 0) + value
            if len(result) > max_terms:
                return atom(expr)
            return {key: value for key, value in result.items() if value != 0}
        if isinstance(expr, L.Mul):
            lhs, rhs = polynomial(expr.lhs), polynomial(expr.rhs)
            if len(lhs) * len(rhs) > max_terms:
                return atom(expr)
            result = {}
            for lk, lv in lhs.items():
                for rk, rv in rhs.items():
                    key = tuple(sorted(lk + rk, key=repr))
                    result[key] = result.get(key, 0) + lv * rv
            return {key: value for key, value in result.items() if value != 0}
        if isinstance(expr, L.Div) and isinstance(expr.rhs, L.LiteralInt | L.LiteralFloat):
            return {key: value / expr.rhs.value for key, value in polynomial(expr.lhs).items()}
        return atom(expr)

    def literal(value):
        return L.LiteralFloat(value) if isinstance(value, complex) else L.as_lexpr(value)

    def product(args):
        return L.as_lexpr(L.float_product(args))

    def build(terms):
        if not terms:
            return L.LiteralFloat(0.0)
        if len(terms) == 1:
            ((key, value),) = terms.items()
            return product([literal(value), *[atoms[arg] for arg in key]])
        counts = Counter(arg for key in terms for arg in sorted(set(key), key=repr))

        def priority(key):
            reuse = uses[atoms[key].name] if isinstance(atoms[key], L.Symbol) else 0
            return counts[key], reuse

        best = max(counts, key=priority, default=None)
        if best is not None and counts[best] > 1:
            group, rest = {}, {}
            for key, value in terms.items():
                if best in key:
                    factors = list(key)
                    factors.remove(best)
                    group[tuple(factors)] = value
                else:
                    rest[key] = value
            term = product([atoms[best], build(group)])
            return L.Sum([term, build(rest)]) if rest else term
        return L.Sum(
            [
                product([literal(value), *[atoms[arg] for arg in key]])
                for key, value in terms.items()
            ]
        )

    for stmt in statements:
        if not isinstance(stmt, L.VariableDecl) or not isinstance(
            stmt.value, L.Add | L.Mul | L.Neg | L.Div
        ):
            output.append(stmt)
            continue
        terms = polynomial(stmt.value)
        value = build(terms)
        key = expr_key(value)
        if key in seen:
            value = seen[key]
        else:
            seen[key] = stmt.symbol
        output.append(L.VariableDecl(stmt.symbol, value))
        if (uses[stmt.symbol.name] == 1 and len(terms) <= max_terms) or (
            len(terms) == 1 and len(next(iter(terms))) <= max_factors
        ):
            expanded[stmt.symbol.name] = terms
    return output
