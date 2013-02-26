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

vt = TestFunction(V0)
ut = TrialFunction(V0)

T = TensorElement("CG", cell, 1)
TS = TensorElement("CG", cell, 1)#, symmetry=True) FIXME
V = VectorElement("CG", cell, 1, dim=4)
RT = FiniteElement("RT", cell, 1)
u = Coefficient(V)
ut = Coefficient(T)
uts = Coefficient(TS)

a1 = (x*xi)*dx
a2 = (detJ*J*K)*dx
a3 = c * (v[0]*v[1]) * (t[0,1]*t[1,0]*t[1,1]) * dx

forms = [a1, a2, a3]
forms = [a3]

#S = FiniteElement("CG", cell, 1)
#V = VectorElement("CG", cell, 1, dim=1)
#us, uvs = Coefficients(S*(V*S))
#forms = [us*uvs[0]*uvs[1]*dx]

#forms = [Constant(cell)*dx]
#forms = [u0*dx]
#forms = [u0.dx(0)*dx]
#forms = [(u[0] + u[1])*dx]
#forms = [(u[0] + u[3])*dx]

#ua, ub = Coefficients(V*V)
#forms = [(ua[0] + ub[3])*dx]

#ua, ub = Coefficients(V*RT)
#forms = [(ua[1] + ub[0] + ub[1])*dx]

#forms = [u0*vt*dx]
#forms = [u0*ut*vt*dx]

#forms = [(u[0] + (ut[0,0] + uts[0,0]) + (ut[1,1] + uts[1,1]))*dx] # FIXME

for form in forms:
    code = uffc.compile_tabulate_tensor_code(form, optimize=True)
    print
    print '/'*60
    print code
    print
