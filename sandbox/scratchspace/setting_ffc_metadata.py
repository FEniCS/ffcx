
rep = 'uflacs.backends.ffc'
r = __import__(rep, fromlist=['compute_integral_ir', 'optimize_integral_ir', 'generate_integral_code'])
print(r.optimize_integral_ir, r.compute_integral_ir, r.generate_integral_code)

from dolfin import *
mesh = UnitSquareMesh(10,10)
V = FunctionSpace(mesh, "CG", 1)
f = Function(V)
M = f*dx(metadata={'representation': 'uflacs.backends.ffc'})
print(assemble(M))
