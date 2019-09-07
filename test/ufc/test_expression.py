import ufl
import ffc


def test_expression():
    e = ufl.VectorElement("P", "triangle", 1)

    f = ufl.Coefficient(e)
    g = ufl.Coefficient(e)

    Te = ufl.TensorElement("P", "triangle", 1)
    u = ufl.TrialFunction(Te)

    a = ufl.as_vector([1.0, 2.0])
    expr = ufl.dot(u, a) * ufl.inner(f, g) + a * u[0, 0]

    code_h, code_c = ffc.compiler.compile_ufl_objects(expr)
