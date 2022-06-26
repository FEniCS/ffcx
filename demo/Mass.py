import ufl

cell = ufl.hexahedron
element = ufl.FiniteElement("Lagrange", cell, 3)
coords = ufl.VectorElement("Lagrange", cell, 1)
mesh = ufl.Mesh(coords)
function_space = ufl.FunctionSpace(mesh, element)
x = ufl.SpatialCoordinate(mesh)

v = ufl.TestFunction(function_space)
u = ufl.TrialFunction(function_space)

a = ufl.inner(u, v) * ufl.dx
