import ufl

cell = ufl.tetrahedron
element = ufl.FiniteElement("Lagrange", cell, 2)
coords = ufl.VectorElement("Lagrange", cell, 1)
mesh = ufl.Mesh(coords)
function_space = ufl.FunctionSpace(mesh, element)
x = ufl.SpatialCoordinate(mesh)

v = ufl.TestFunction(function_space)
u = ufl.TrialFunction(function_space)

a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
