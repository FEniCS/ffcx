import ufl

cell = ufl.hexahedron
element = ufl.VectorElement("Lagrange", cell, 2)
coords = ufl.VectorElement("Lagrange", cell, 1)
mesh = ufl.Mesh(coords)
V = ufl.FunctionSpace(mesh, element)
x = ufl.SpatialCoordinate(mesh)

v = ufl.TestFunction(V)
u = ufl.Coefficient(V)

a = ufl.inner(u, v) * ufl.dx
