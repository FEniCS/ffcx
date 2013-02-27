import sys
sys.path.insert(0, "../../site-packages")
uflacs = __import__("uflacs")

import uflacs.backends.ffc as uffc
from ufl import *

cell = triangle

x = cell.x
xi = cell.xi
J = cell.J
detJ = cell.detJ
K = cell.Jinv

S0 = FiniteElement("DG", cell, 0)
S1 = FiniteElement("CG", cell, 1)
S2 = FiniteElement("CG", cell, 2)
V0 = VectorElement("DG", cell, 0)
V1 = VectorElement("CG", cell, 1)
V2 = VectorElement("CG", cell, 2)
T0 = TensorElement("DG", cell, 0)
T1 = TensorElement("CG", cell, 1)
T2 = TensorElement("CG", cell, 2)
TS0 = TensorElement("DG", cell, 0, symmetry=True)
TS1 = TensorElement("CG", cell, 1, symmetry=True)
TS2 = TensorElement("CG", cell, 2, symmetry=True)

RT1 = FiniteElement("RT", cell, 1)
RT2 = FiniteElement("RT", cell, 2)

cs = Constant(cell)
cv = VectorConstant(cell)
ct = TensorConstant(cell)
cts = TensorConstant(cell, symmetry=True)

s0 = Coefficient(S0)
s1 = Coefficient(S1)
s2 = Coefficient(S2)
v0 = Coefficient(V0)
v1 = Coefficient(V1)
v2 = Coefficient(V2)
t0 = Coefficient(T0)
t1 = Coefficient(T1)
t2 = Coefficient(T2)
ts0 = Coefficient(TS0)
ts1 = Coefficient(TS1)
ts2 = Coefficient(TS2)

us0 = TrialFunction(S0)
us1 = TrialFunction(S1)
us2 = TrialFunction(S2)
uv0 = TrialFunction(V0)
uv1 = TrialFunction(V1)
uv2 = TrialFunction(V2)
ut0 = TrialFunction(T0)
ut1 = TrialFunction(T1)
ut2 = TrialFunction(T2)
uts0 = TrialFunction(TS0)
uts1 = TrialFunction(TS1)
uts2 = TrialFunction(TS2)

vs0 = TestFunction(S0)
vs1 = TestFunction(S1)
vs2 = TestFunction(S2)
vv0 = TestFunction(V0)
vv1 = TestFunction(V1)
vv2 = TestFunction(V2)
vt0 = TestFunction(T0)
vt1 = TestFunction(T1)
vt2 = TestFunction(T2)
vts0 = TestFunction(TS0)
vts1 = TestFunction(TS1)
vts2 = TestFunction(TS2)

# Build list of forms
forms = []

# Geometry forms:
a_xi = (xi[0]*xi[1])*dx
a_x = (x[0]*x[1])*dx
a_J = (tr(J*K) / (3*detJ))*dx
forms += [a_x, a_xi, a_J]

# Constant access forms
a_c = (cs + (cv[0]*cv[1]) + ((cts[0,1]*cts[1,0])*(cts[0,0]*cts[1,1])) + ((ct[0,1]*ct[1,0])*(ct[0,0]*ct[1,1]))) * dx
forms += [a_c]

# Coefficient access forms
a_w0 = (s0 + (v0[0]*v0[1]) + ((ts0[0,1]*ts0[1,0])*(ts0[0,0]*ts0[1,1])) + ((t0[0,1]*t0[1,0])*(t0[0,0]*t0[1,1]))) * dx
a_w1 = (s1 + (v1[0]*v1[1]) + ((ts1[0,1]*ts1[1,0])*(ts1[0,0]*ts1[1,1])) + ((t1[0,1]*t1[1,0])*(t1[0,0]*t1[1,1]))) * dx
forms += [a_w0, a_w1]

# Mixed element coefficient access forms
ums, umvs = Coefficients(S1*(V0*S1))
a_m1 = ums*umvs[0]*umvs[1]*dx
umva, umvb = Coefficients(V1*V2)
a_m2 = (umva[0] + umvb[1])*dx
umvc, umrt = Coefficients(V1*RT1)
a_m3 = (umvc[1] + umrt[0] + umrt[1])*dx
forms += [a_m1, a_m2, a_m3]

# Coefficient derivative access forms
a_d0 = s1.dx(0)*dx
a_d1 = s1.dx(1)*dx
a_d2 = (s1 + s1.dx(0) + s1.dx(1))*dx
forms += [a_d0, a_d1, a_d2]

# Argument access forms
a_a0 = us0*dx
a_a1 = (vv1[0] + vv1[1])*dx
a_a2 = us2*(vts0[0,0] + vts0[0,1] + vts0[1,0] + vts0[1,1])*dx
forms += [a_a0, a_a1, a_a2]

# Argument derivative access forms
a_da0 = us0.dx(0)*dx
a_da1 = (vv1[0].dx(0) + vv1[1].dx(1))*dx
a_da2 = us2.dx(1)*(vts1[0,0].dx(0) + vts1[0,1].dx(0) + vts1[1,0].dx(1) + vts1[1,1].dx(1))*dx
forms += [a_da0, a_da1, a_da2]

for form in forms:
    code = uffc.compile_tabulate_tensor_code(form, optimize=True)
    print
    print ('/'*40), str(form)
    print code
    print

print "Done compiling %d example forms." % len(forms)

