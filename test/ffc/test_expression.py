import ufl

import ffc

ffc.logger.setLevel(ffc.logging.DEBUG)


def test_scalar_expression():
    e = ufl.FiniteElement("P", ufl.triangle, 2)
    v = ufl.Coefficient(e)

    expr = v ** 2 + ufl.exp(v)

    code_c, code_h = ffc.compiler.compile_ufl_objects([expr], prefix="JIT")
    print(code_c)
