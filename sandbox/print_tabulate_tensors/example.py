import sys
sys.path.insert(0, "../../site-packages")
uflacs = __import__("uflacs")

import uflacs.backends.ffc as uffc
from ufl import *

cell = triangle

x = cell.x[0]
xi = cell.xi[0]
J = cell.J[0,0]
detJ = cell.detJ
K = cell.Jinv[0,0]

c = Constant(cell)
v = VectorConstant(cell)
t = TensorConstant(cell)#, symmetry=True) # symmetry is not fully implemented in uflacs

V0 = FiniteElement("CG", cell, 1)
u0 = Coefficient(V0)

V = VectorElement("CG", cell, 1)
u = Coefficient(V)

a1 = (x*xi)*dx
a2 = (detJ*J*K)*dx
a3 = c*v[0]*t[0,1]*t[1,0]*dx

forms = [a1, a2, a3]
forms = [a3]

#forms = [Constant(cell)*dx]
forms = [u0*dx]
forms = [(u[0] + u[1])*dx]

for form in forms:
    code = uffc.compile_tabulate_tensor_code(form, optimize=True)
    print
    print '/'*60
    print code
    print
