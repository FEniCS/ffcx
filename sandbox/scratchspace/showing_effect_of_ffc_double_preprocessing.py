from dolfin import *

mesh = UnitSquareMesh(10, 10)

e = Expression(("x[0]", "x[1]"))
print("original degree", e.ufl_element().degree())

f = dot(e, e)*dx
val = assemble(f, mesh=mesh)
print("value", val)

fd = f.compute_form_data()
V2 = fd.renumbered_coefficients[0].element()
print("new degree", V2.degree())
print(fd.function_replace_map)

