"""Experimental reciprocal-scaling transformations for generated kernels."""

# Generated C fragments retain their original source layout.

import ffcx.codegeneration.lnodes as L
from ffcx.codegeneration.optimizer import expr_key


class Scaling:
    """Track homogeneous determinant factors in experimental scalar code."""

    def __init__(self):
        """Initialise the experimental scaling state."""
        self.groups = {}
        self.number = 0

    def __call__(self, statements):
        """Factor reciprocal powers in a flat scalar statement list."""

        def walk(e):
            """Apply one step of the experimental transformation."""
            yield e
            if isinstance(e, L.BinOp):
                yield from walk(e.lhs)
                yield from walk(e.rhs)
            elif isinstance(e, L.PrefixUnaryOp):
                yield from walk(e.arg)
            elif isinstance(e, (L.NaryOp, L.MathFunction)):
                for a in e.args:
                    yield from walk(a)

        for s in statements:
            if not isinstance(s, L.VariableDecl) or s.value is None:
                continue
            for e in walk(s.value):
                if (
                    isinstance(e, L.MathFunction)
                    and e.function == "abs"
                    and e.args[0].dtype == L.DataType.REAL
                ):
                    key = expr_key(e.args[0])
                    if key not in self.groups:
                        self.groups[key] = {
                            "divisor": e.args[0],
                            "absnames": set(),
                            "normal": {},
                            "values": {},
                            "recip": None,
                            "geometry": set(),
                        }
                    if e is s.value:
                        self.groups[key]["absnames"].add(s.symbol.name)
        output = statements
        for group in self.groups.values():
            output = self.run(output, group)
        return output

    def run(self, statements, g):
        """Rewrite scalar statements for one determinant."""
        norm = g["normal"]
        values = g["values"]
        absnames = g["absnames"]
        key = expr_key(g["divisor"])
        output = []

        def mul(a, b):
            """Apply one step of the experimental transformation."""
            return L.as_lexpr(L.float_product([a, b]))

        def normalize(e):
            """Apply one step of the experimental transformation."""
            if isinstance(e, L.Symbol):
                return norm.get(e.name, (e, 0))
            if (
                isinstance(e, L.Div)
                and isinstance(e.lhs, (L.LiteralFloat, L.LiteralInt))
                and e.lhs.value == 1
                and expr_key(e.rhs) == key
            ):
                return L.LiteralFloat(1.0), 1
            if isinstance(e, L.Product):
                terms = [normalize(a) for a in e.args]
                degree = sum(d for a, d in terms)
                return (
                    (L.as_lexpr(L.float_product([a for a, d in terms])), degree)
                    if degree <= 4
                    else (e, 0)
                )
            if isinstance(e, L.Sum):
                terms = [normalize(a) for a in e.args]
                if len({d for a, d in terms}) == 1:
                    return L.Sum([a for a, d in terms]), terms[0][1]
                return e, 0
            if isinstance(e, L.Mul):
                a, da = normalize(e.lhs)
                b, db = normalize(e.rhs)
                return (mul(a, b), da + db) if da + db <= 4 else (e, 0)
            if isinstance(e, L.Neg):
                a, d = normalize(e.arg)
                return -a, d
            if isinstance(e, L.Add):
                a, da = normalize(e.lhs)
                b, db = normalize(e.rhs)
                if da == db:
                    return a + b, da
            return e, 0

        def factors(e, depth=0):
            """Apply one step of the experimental transformation."""
            if isinstance(e, L.Product):
                return [f for a in e.args for f in factors(a, depth)]
            if isinstance(e, L.Mul):
                return factors(e.lhs, depth) + factors(e.rhs, depth)
            if isinstance(e, L.Symbol) and e.name in values and depth < 4:
                v = values[e.name]
                if isinstance(v, (L.Mul, L.Product, L.Symbol)):
                    return factors(v, depth + 1)
            return [e]

        def geom(e):
            """Apply one step of the experimental transformation."""
            if isinstance(e, L.Symbol):
                return e.name.startswith(("J", "K")) or e.name in g["geometry"]
            if isinstance(e, (L.LiteralFloat, L.LiteralInt)):
                return True
            if isinstance(e, L.BinOp):
                return geom(e.lhs) and geom(e.rhs)
            if isinstance(e, (L.NaryOp, L.MathFunction)):
                return all(geom(a) for a in e.args)
            if isinstance(e, L.Neg):
                return geom(e.arg)
            return False

        for s in statements:
            assigned = (
                isinstance(s, L.Statement)
                and not isinstance(s, L.VariableDecl)
                and isinstance(getattr(s, "expr", None), L.Assign)
                and isinstance(s.expr.lhs, L.Symbol)
            )
            if assigned:
                s = L.VariableDecl(s.expr.lhs, s.expr.rhs)
            if not isinstance(s, L.VariableDecl) or s.value is None:
                output.append(s)
                continue
            if geom(s.value):
                g["geometry"].add(s.symbol.name)
            e, degree = normalize(s.value)
            if degree == 1 and not geom(e):
                # Preserve precomputed inverse-Jacobian reuse in coefficient gradients.
                output.append(L.Statement(L.Assign(s.symbol, s.value)) if assigned else s)
                norm[s.symbol.name] = (s.symbol, 0)
                continue
            if not degree:
                output.append(L.Statement(L.Assign(s.symbol, s.value)) if assigned else s)
                norm[s.symbol.name] = (s.symbol, 0)
                if s.symbol.name in g["geometry"]:
                    values[s.symbol.name] = s.value
                continue
            if g["recip"] is None:
                g["recip"] = L.Symbol(f"scale_recip_{self.number}", L.DataType.REAL)
                self.number += 1
                output.append(L.VariableDecl(g["recip"], L.Div(L.LiteralFloat(1.0), g["divisor"])))
            name = L.Symbol("scale_" + s.symbol.name, e.dtype)
            if s.symbol.name in g["geometry"]:
                g["geometry"].add(name.name)
            output.append(L.VariableDecl(name, e))
            values[name.name] = e
            norm[s.symbol.name] = (name, degree)
            fs = factors(e)
            at = next(
                (j for j, a in enumerate(fs) if isinstance(a, L.Symbol) and a.name in absnames),
                None,
            )
            if degree >= 2 and at is not None:
                fs.pop(at)
                scale = L.MathFunction("abs", [g["recip"]])
                degree -= 2
                for _ in range(degree):
                    scale = mul(scale, g["recip"])
                value = L.as_lexpr(L.float_product([scale, *fs]))
            else:
                scale = g["recip"]
                for _ in range(degree - 1):
                    scale = mul(scale, g["recip"])
                value = mul(scale, name)
            output.append(
                L.Statement(L.Assign(s.symbol, value))
                if assigned
                else L.VariableDecl(s.symbol, value)
            )
        return output


def install():
    """Install the experimental generator monkeypatch."""
    from ffcx.codegeneration.integral_generator import IntegralGenerator

    original = IntegralGenerator.generate_partition

    def generate(self, *args, **kwargs):
        """Apply one step of the experimental transformation."""
        definitions, intermediates = original(self, *args, **kwargs)
        if not hasattr(self, "_scaling"):
            self._scaling = Scaling()
        return definitions, self._scaling(intermediates)

    IntegralGenerator.generate_partition = generate


def install_full():
    """Install the experimental whole-kernel scaling monkeypatch."""
    from ffcx.codegeneration.integral_generator import IntegralGenerator
    from ffcx.codegeneration.optimizer import prune_dead_scalars, scalar_declaration_names

    original = IntegralGenerator.generate

    def generate(self, *args, **kwargs):
        """Apply one step of the experimental transformation."""
        code = original(self, *args, **kwargs)
        engine = Scaling()

        def flatten(nodes):
            """Apply one step of the experimental transformation."""
            out = []
            for n in nodes:
                if isinstance(n, (L.Section, L.StatementList)):
                    out += flatten(n.statements)
                elif isinstance(n, L.ForRange):
                    out += flatten(n.body.statements)
                else:
                    out.append(n)
            return out

        engine(flatten(code.statements))
        for g in engine.groups.values():
            g.update(normal={}, values={}, geometry=set(), recip=None)
        engine.number = 0

        def block(nodes, scoped=False):
            """Apply one step of the experimental transformation."""
            saved = {
                k: {
                    **g,
                    "normal": g["normal"].copy(),
                    "values": g["values"].copy(),
                    "geometry": g["geometry"].copy(),
                }
                for k, g in engine.groups.items()
            }
            output = []
            flat = []

            def flush():
                """Apply one step of the experimental transformation."""
                if flat:
                    output.extend(engine(flat))
                    flat.clear()

            for n in nodes:
                if isinstance(n, L.Section):
                    flush()
                    n.statements = block(n.statements, True)
                    output.append(n)
                elif isinstance(n, L.ForRange):
                    flush()
                    n.body.statements = block(n.body.statements, True)
                    output.append(n)
                elif isinstance(n, L.StatementList):
                    flush()
                    n.statements = block(n.statements)
                    output.append(n)
                else:
                    flat.append(n)
            flush()
            if scoped:
                for k, g in engine.groups.items():
                    if k in saved:
                        g.update(saved[k])
            return output

        parts = block(code.statements)
        code.statements = prune_dead_scalars(parts, scalar_declaration_names(parts))
        return code

    IntegralGenerator.generate = generate
