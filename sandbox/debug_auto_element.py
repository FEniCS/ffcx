
"""
for n in 1 2; do instant-clean; python debug_auto_element.py $n; done
"""
from dolfin import *
import sys

parameters["form_compiler"]["representation"] = "uflacs"

mesh = UnitSquareMesh(16, 16)
V = FunctionSpace(mesh, 'CG', 2)
#w = TestFunction(V)
#Pv = TrialFunction(V)
#a = inner(w, Pv)*dx
#A = assemble(a)
#u = Function(V)

if "0" in sys.argv:
    f = Expression("x[0]*x[0] + 2*x[1]*x[1]", degree=2)
    L = f*dx
    print("Value:", assemble(L, mesh=mesh))
    print("Form data")
    print(L.form_data())

if "1" in sys.argv:
    f = Expression("x[0]*x[0] + 2*x[1]*x[1]")
    L = f*dx
    print("Value:", assemble(L, mesh=mesh))
    print("Form data")
    print(L.form_data())

