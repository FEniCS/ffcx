import SyFi, SFC

from SFC.symbolic_utils import div, grad, inner 
from SFC import MatrixForm, compile_form

SFC.options.include_from = "pycc"


def conv_form(v, u, w, rho, G, Ginv):
    return rho*inner(inner(w, grad(u, Ginv)),v)

C  = MatrixForm(conv_form)

SyFi.initSyFi(2)
polygon = SyFi.ReferenceTriangle()
v_fe = SyFi.VectorLagrange(polygon,1)
u_fe = SyFi.VectorLagrange(polygon,1)
w_fe = SyFi.VectorP0(polygon,0)
rho_fe = SyFi.P0(polygon,1)

compiled_conv_form = compile_form(C, [v_fe, u_fe, w_fe, rho_fe])

