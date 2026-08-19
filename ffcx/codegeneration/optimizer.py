"""Optimizer."""

from collections import defaultdict

import ffcx.codegeneration.lnodes as L
from ffcx.ir.representationutils import QuadratureRule


def _key(expr) -> tuple:
    """Hashable structural identity for an LExpr, for recognising equal factors.

    LNodes expressions are not generally hashable in a way that reflects
    their mathematical content, so build a key that compares equal exactly
    when two expressions are the same expression.

    Args:
        expr: Expression to key.

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


def _get_mul_operands(stmt: L.LNode) -> tuple[L.LExpr, L.LExpr] | None:
    """Return (lhs, rhs) if stmt is a simple scalar ``VariableDecl(sym, Mul(lhs, rhs))``."""
    if isinstance(stmt, L.VariableDecl) and isinstance(stmt.value, L.Mul):
        return stmt.value.lhs, stmt.value.rhs
    return None


def _detect_dense_product_run(statements: list[L.LNode], start: int):
    """Look for a maximal run of Mul declarations forming a dense outer product.

    Many block-coupled bilinear forms (e.g. one tensor/vector component
    multiplying another) produce, for every pair drawn from two small sets
    of "left" and "right" scalar factors, a separately declared scalar
    holding their product -- n0 * n1 declarations where n0 + n1 would do,
    since the two factor sets are each already computed once. This detects
    such a run when it's laid out in nested (left, right) order: n1
    consecutive statements sharing one left operand and enumerating every
    right operand, repeated for each of n0 distinct left operands, with the
    same right-operand sequence every time.

    Args:
        statements: Flat statement list to scan.
        start: Index to start looking for a run at.

    Returns:
        (n0, n1, left_factors, right_factors, entry_symbols) if a dense run
        of at least 2x2 is found starting at `start`, else None.
    """
    first = _get_mul_operands(statements[start])
    if first is None:
        return None
    first_left, first_right = first
    left_key0 = _key(first_left)

    j = start
    right_factors: list[L.LExpr] = []
    right_keys = []
    while j < len(statements):
        mul = _get_mul_operands(statements[j])
        if mul is None or _key(mul[0]) != left_key0:
            break
        right_factors.append(mul[1])
        right_keys.append(_key(mul[1]))
        j += 1
    n1 = j - start
    if n1 < 2:
        return None

    left_factors = [first_left]
    entry_symbols = [statements[k].symbol for k in range(start, j)]
    pos = j
    while pos < len(statements):
        row = _get_mul_operands(statements[pos])
        if row is None:
            break
        row_left_key = _key(row[0])
        if row_left_key == left_key0:
            break
        row_end = pos + n1
        if row_end > len(statements):
            break
        ok = True
        for k in range(n1):
            m = _get_mul_operands(statements[pos + k])
            if m is None or _key(m[0]) != row_left_key or _key(m[1]) != right_keys[k]:
                ok = False
                break
        if not ok:
            break
        left_factors.append(row[0])
        entry_symbols.extend(statements[pos + k].symbol for k in range(n1))
        pos = row_end

    n0 = len(left_factors)
    if n0 < 2:
        return None
    return n0, n1, left_factors, right_factors, entry_symbols


def _substitute_symbols(expr: L.LNode, remap: dict[str, L.LExpr]):
    """Rebuild expr, replacing any Symbol whose name is a key of remap.

    Statements later in the same list can already reference (by name) a
    scalar this pass is about to fold into an array, so removing its
    declaration alone would leave a dangling reference; this rewrites those
    uses in place.

    Args:
        expr: Expression or statement to rebuild.
        remap: Symbol name -> replacement LExpr.

    Returns:
        A structurally equivalent node with substitutions applied.
    """
    if isinstance(expr, L.Symbol):
        return remap.get(expr.name, expr)
    if isinstance(expr, (L.LiteralFloat, L.LiteralInt, L.MultiIndex)):
        return expr
    if isinstance(expr, L.ArrayAccess):
        array = _substitute_symbols(expr.array, remap)
        indices = [_substitute_symbols(idx, remap) for idx in expr.indices]
        return L.ArrayAccess(array, indices)
    if isinstance(expr, L.VariableDecl):
        return L.VariableDecl(expr.symbol, _substitute_symbols(expr.value, remap))
    if isinstance(expr, L.AssignOp):
        return type(expr)(
            _substitute_symbols(expr.lhs, remap), _substitute_symbols(expr.rhs, remap)
        )
    if isinstance(expr, L.BinOp):
        return type(expr)(
            _substitute_symbols(expr.lhs, remap), _substitute_symbols(expr.rhs, remap)
        )
    if isinstance(expr, L.Conditional):
        return L.Conditional(
            _substitute_symbols(expr.condition, remap),
            _substitute_symbols(expr.true, remap),
            _substitute_symbols(expr.false, remap),
        )
    if isinstance(expr, L.MathFunction):
        return L.MathFunction(expr.function, [_substitute_symbols(a, remap) for a in expr.args])
    if isinstance(expr, L.NaryOp):
        return type(expr)([_substitute_symbols(a, remap) for a in expr.args])
    if isinstance(expr, L.PrefixUnaryOp):
        return type(expr)(_substitute_symbols(expr.arg, remap))
    # Unknown/opaque node (e.g. a table Symbol used directly): nothing to do.
    return expr


def fuse_outer_products(
    statements: list[L.LNode], name_prefix: str = "outer"
) -> tuple[list[L.LNode], dict[str, L.LExpr]]:
    """Replace dense n0 x n1 product runs with two small arrays and one loop.

    A block-coupled bilinear form (e.g. linear elasticity's stress-strain
    coupling between vector components) typically needs, for every pair of
    "left" and "right" scalar factors drawn from two small sets, their
    product as a separate scalar -- n0 * n1 declarations, one per pair. Since
    each factor is itself already computed once, this is n0 * n1 scalar
    multiplies where n0 + n1 loads plus one small nested loop would do, and
    a flat run of n0 * n1 unrelated-looking declarations is exactly the kind
    of code GCC's vectoriser and scheduler have nothing to work with, versus
    one tight loop.

    Args:
        statements: A flat list of statements (piecewise/varying partition,
            no nested loops expected on input).
        name_prefix: Prefix for generated array/loop-index symbols, kept
            distinct per call site so two calls in the same kernel can't
            clash.

    Returns:
        (new_statements, remap) where remap maps the name of every scalar
        symbol whose declaration was replaced to the L.ArrayAccess expression
        that now holds its value -- callers that track a UFL-expression to
        LNode-access mapping for these symbols need to redirect it.
    """
    new_statements: list[L.LNode] = []
    remap: dict[str, L.LExpr] = {}
    i = 0
    counter = 0
    while i < len(statements):
        run = _detect_dense_product_run(statements, i)
        if run is None:
            new_statements.append(statements[i])
            i += 1
            continue

        n0, n1, left_factors, right_factors, entry_symbols = run

        left_arr = L.Symbol(f"{name_prefix}_l{counter}", L.DataType.SCALAR)
        right_arr = L.Symbol(f"{name_prefix}_r{counter}", L.DataType.SCALAR)
        out_arr = L.Symbol(f"{name_prefix}_{counter}", L.DataType.SCALAR)
        i0 = L.Symbol(f"{name_prefix}_i{counter}", L.DataType.INT)
        i1 = L.Symbol(f"{name_prefix}_j{counter}", L.DataType.INT)
        counter += 1

        new_statements.append(L.ArrayDecl(left_arr, n0))
        for k0, factor in enumerate(left_factors):
            new_statements.append(L.Assign(L.ArrayAccess(left_arr, [k0]), factor))
        new_statements.append(L.ArrayDecl(right_arr, n1))
        for k1, factor in enumerate(right_factors):
            new_statements.append(L.Assign(L.ArrayAccess(right_arr, [k1]), factor))
        new_statements.append(L.ArrayDecl(out_arr, (n0, n1)))
        body = L.Assign(
            L.ArrayAccess(out_arr, [i0, i1]),
            L.Mul(L.ArrayAccess(left_arr, [i0]), L.ArrayAccess(right_arr, [i1])),
        )
        new_statements.append(L.ForRange(i0, 0, n0, [L.ForRange(i1, 0, n1, [body])]))

        run_remap: dict[str, L.LExpr] = {}
        for k0 in range(n0):
            for k1 in range(n1):
                old_symbol = entry_symbols[k0 * n1 + k1]
                run_remap[old_symbol.name] = L.ArrayAccess(out_arr, [k0, k1])
        remap.update(run_remap)

        # Statements after this run may already reference one of the folded
        # symbols by name (e.g. a later accumulation combining this run's
        # entries with another quadrature/index-sum contribution); rewrite
        # those uses in place rather than leaving a dangling reference. A
        # symbol can only be used after its own declaration in this
        # topologically-ordered list, so the run's own statements (already
        # emitted above) never need this.
        tail_start = i + n0 * n1
        statements[tail_start:] = [
            _substitute_symbols(stmt, run_remap) for stmt in statements[tail_start:]
        ]

        i = tail_start

    return new_statements, remap


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
            if L.Annotation.fuse in section.annotations:
                section = fuse_loops(section)
            if L.Annotation.licm in section.annotations:
                section = licm(section, quadrature_rule)
            code[i] = section

    return code


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
    for lhs, rhs in expressions.items():
        for r in rhs:
            hoist_candidates = []
            for arg in r.args:
                dependency = check_dependency(arg, inner_loop.index)
                if not dependency:
                    hoist_candidates.append(arg)
            if len(hoist_candidates) > 1:
                # create new temp
                name = f"temp_{counter}"
                counter += 1
                temp = L.Symbol(name, L.DataType.SCALAR)
                for h in hoist_candidates:
                    r.args.remove(h)
                # update expression with new temp
                r.args.append(L.ArrayAccess(temp, [outer_loop.index]))
                # create code for hoisted term
                size = outer_loop.end.value - outer_loop.begin.value
                pre_loop.append(L.ArrayDecl(temp, size, [0]))
                body = L.Assign(
                    L.ArrayAccess(temp, [outer_loop.index]), L.Product(hoist_candidates)
                )
                pre_loop.append(
                    L.ForRange(outer_loop.index, outer_loop.begin, outer_loop.end, [body])
                )

    section.statements = pre_loop + section.statements

    return section
