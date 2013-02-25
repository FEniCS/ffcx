from dolfin import *
parameters["form_compiler"]["representation"] = "uflacs"

mesh = UnitSquare(3,3)
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
print assemble(v*dx)
print assemble(u*v*dx)
print assemble(c*v*dx)

