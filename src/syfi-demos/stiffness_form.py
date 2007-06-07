import SyFi, SFC

from SFC.symbolic_utils import grad, inner 
from SFC import MatrixForm, compile_form

SFC.options.include_from = "pycc"


def stiffness_form(u, v, G, Ginv):
    grad_v = grad(v, Ginv.transpose())                               
    grad_u = grad(u, Ginv.transpose())                               
    return inner(grad_v, grad_u) 

A  = MatrixForm(stiffness_form)

SyFi.initSyFi(2)
polygon = SyFi.ReferenceTriangle()
fe = SyFi.Lagrange(polygon,1)

compiled_stiffness_form = compile_form(A, fe)

