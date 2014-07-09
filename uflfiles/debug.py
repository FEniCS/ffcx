from dolfin import *
import sys

rep = sys.argv[2]
if rep == 'u':
    rep = 'uflacs'
elif rep == 'q':
    rep = 'quadrature'
n = int(sys.argv[1])
parameters["form_compiler"]["representation"] = rep

mesh = UnitSquareMesh(n, n)
domain = mesh.ufl_domain()
V = FunctionSpace(domain, "Lagrange", 1)
v = TestFunction(V)
u = TrialFunction(V)

a = dot(grad(u),grad(v))*dx(domain)
A = assemble(a)
print("A:", sum(sum(A.array()[:,:]**2)))
