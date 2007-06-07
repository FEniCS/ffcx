import SyFi, SFC

from SFC.symbolic_utils import div 
from SFC import MatrixForm, compile_form

SFC.options.include_from = "pycc"


def div_form(v, q, G, Ginv):
    div_v = div(v, Ginv.transpose())                               
    return div_v*q  

B  = MatrixForm(div_form)

SyFi.initSyFi(2)
polygon = SyFi.ReferenceTriangle()
v_fe = SyFi.VectorLagrange(polygon,2)
q_fe = SyFi.Lagrange(polygon,1)

compiled_div_form = compile_form(B, [v_fe, q_fe])

