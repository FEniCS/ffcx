import ufl

cell = ufl.hexahedron
element = ufl.VectorElement("Lagrange", cell, 5)
coords = ufl.VectorElement("Lagrange", cell, 1)
mesh = ufl.Mesh(coords)
function_space = ufl.FunctionSpace(mesh, element)
x = ufl.SpatialCoordinate(mesh)

v = ufl.TestFunction(function_space)
u = ufl.Coefficient(function_space)

a = ufl.inner(u, v) * ufl.dx
