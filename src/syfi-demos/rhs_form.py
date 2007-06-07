import SyFi, SFC

from SFC.symbolic_utils import div 
from SFC import VectorForm, compile_form

SFC.options.include_from = "pycc"


def rhs_form(v, f, G, Ginv):
    return v*f 

F  = VectorForm(rhs_form)

SyFi.initSyFi(2)
polygon = SyFi.ReferenceTriangle()
v_fe = SyFi.Lagrange(polygon,1)
f_fe = SyFi.Lagrange(polygon,1)

compiled_rhs_form = compile_form(F, [v_fe, f_fe])

