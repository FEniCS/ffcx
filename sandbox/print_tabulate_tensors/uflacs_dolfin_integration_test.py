from dolfin import *
import numpy
import sys
if sys.argv[1] == 'u':
    parameters["form_compiler"]["representation"] = "uflacs"
print "Using form compiler representation:", parameters["form_compiler"]["representation"]

n = int(sys.argv[2])
mesh = UnitSquareMesh(n,3)
V = FunctionSpace(mesh, "CG", 1)
c = Constant(2.3)
f = Function(V)
e = Expression("x[0]+x[1]")
f.interpolate(e)
u = TrialFunction(V)
v = TestFunction(V)

print assemble(c*dx, mesh=mesh)
print assemble(f*dx)
print assemble(f**2*dx)
print assemble(c*f*dx)

print assemble(v*dx).norm('linf')
print assemble(u*v*dx).norm('linf')
print assemble(c*v*dx).norm('linf')

