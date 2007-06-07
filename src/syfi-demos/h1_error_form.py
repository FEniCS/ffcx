import SyFi, SFC

from SFC.symbolic_utils import div, inner 
from SFC import ScalarForm, compile_form

SFC.options.include_from = "pycc"


def H1_error(u, u_h, G, Ginv):
    return inner(u-u_h, u-u_h) 

E  = ScalarForm(H1_error)

SyFi.initSyFi(2)
polygon = SyFi.ReferenceTriangle()
u_fe = SyFi.Lagrange(polygon,1)
uh_fe = SyFi.Lagrange(polygon,1)

compiled_error_form = compile_form(E, [u_fe, uh_fe] )

