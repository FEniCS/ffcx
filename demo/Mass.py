import ufl

cell = ufl.hexahedron
element = ufl.FiniteElement("Lagrange", cell, 3)
coords = ufl.VectorElement("Lagrange", cell, 1)
mesh = ufl.Mesh(coords)
V = ufl.FunctionSpace(mesh, element)
x = ufl.SpatialCoordinate(mesh)

v = ufl.TestFunction(V)
u = ufl.Coefficient(V)

a = ufl.inner(u, v) * ufl.dx
