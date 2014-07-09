from dolfin import *

set_log_level(DEBUG)

parameters["form_compiler"]["cpp_optimize"] = True
#parameters["form_compiler"]["representation"] = "uflacs"
compiler_options = {"representation": "uflacs", "quadrature_degree": 2}

mesh = Mesh("lv_mesh.xml")
bi = MeshFunction("uint", mesh, "lv_mesh_bi.xml")

# Hack to convert boundary indicators from diffpack mesh
ff = FacetFunction("uint", mesh)
ff.set_all(0)
append7 = 0
append7b = 0
append8 = 0
for facet in facets(mesh):
    verts = vertices(facet)
    ve = [v.index() for v in verts]
    if 578 in ve and append7 < 3:
        if append7b > 2:
            ff[facet.index()] = 7
            append7 += 1
        else:
            append7b += 1
        continue
    if 560 in ve and append8 < 2:
        ff[facet.index()] = 8
        append8 += 1
        continue
    values = [bi.array()[l1] for l1 in ve]
    if all(values[i] == values[2] for i in xrange(2)):
        ff[facet.index()] = values[0]
    else:
        ff[facet.index()] = max(values)

# Function space and boundary conditions
V = VectorFunctionSpace(mesh, "CG", 1)
bc1 = DirichletBC(V.sub(0), 0.0, ff, 6)
bc2 = DirichletBC(V.sub(0), 0.0, ff, 7)
bc3 = DirichletBC(V.sub(1), 0.0, ff, 7)
bc4 = DirichletBC(V.sub(2), 0.0, ff, 7)
bc5 = DirichletBC(V.sub(0), 0.0, ff, 8)
bc6 = DirichletBC(V.sub(2), 0.0, ff, 8)
bcs = [bc1, bc2, bc3, bc4, bc5, bc6]


conductivity_code = """

class Conductivity : public Expression
{
public:

  Conductivity() : Expression(9) {}

  // Function for evaluating expression on each cell
  //void eval(Array<double>& values, const Data& data) const
  void eval(Array<double>& values, const Array<double>& x,  const ufc::cell& cell) const
  {
    const uint D = cell.topological_dimension;
    const uint cell_index = cell.entity_indices[D][0];

    values[0] = (*f1)[cell_index];
    values[1] = (*f2)[cell_index];
    values[2] = (*f3)[cell_index];
    values[3] = (*n1)[cell_index];
    values[4] = (*n2)[cell_index];
    values[5] = (*n3)[cell_index];

    values[6] =  values[1]*values[5] - values[2]*values[4];
    values[7] = -values[0]*values[5] + values[2]*values[3];
    values[8] =  values[0]*values[4] - values[1]*values[3];

  }

  // The data stored in mesh functions
  boost::shared_ptr<MeshFunction<double> > f1;
  boost::shared_ptr<MeshFunction<double> > f2;
  boost::shared_ptr<MeshFunction<double> > f3;
  boost::shared_ptr<MeshFunction<double> > n1;
  boost::shared_ptr<MeshFunction<double> > n2;
  boost::shared_ptr<MeshFunction<double> > n3;


};
"""

c00 = MeshFunction("double", mesh, "f1.xml.gz")
c01 = MeshFunction("double", mesh, "f2.xml.gz")
c11 = MeshFunction("double", mesh, "f3.xml.gz")
b00 = MeshFunction("double", mesh, "n1.xml.gz")
b01 = MeshFunction("double", mesh, "n2.xml.gz")
b11 = MeshFunction("double", mesh, "n3.xml.gz")

c = Expression(cppcode=conductivity_code)
c.f1 = c00
c.f2 = c01
c.f3 = c11
c.n1 = b00
c.n2 = b01
c.n3 = b11

R = as_matrix([[c[0], c[6], c[3]],
               [c[1], c[7], c[4]],
               [c[2], c[8], c[5]]])


# Functions
v  = TestFunction(V)
u  = Function(V)
du = TrialFunction(V)

# Parameters
K   = Constant(0.876)
bff = Constant(18.48)
bxx = Constant(3.58)
bfx = Constant(2.8)
Ccompr = Constant(10.0)

# Identity matrix and global deformation gradient
I = Identity(3)
F_glob = I + grad(u)
F = variable(R.T*F_glob*R)

E = 0.5*(F.T*F - I)
J = det(F)

# Material law
f=0; s=1; n=2
W = (bff*E[f,f]**2
     + bxx*(E[n,n]**2 + E[s,s]**2 + E[n,s]**2)
     + bfx*(E[f,n]**2 + E[n,f]**2 + E[f,s]**2 + E[s,f]**2))
psi = 0.5*K*(exp(W) - 1) + Ccompr*(J*ln(J) - J + 1)

# PK1
P = diff(psi, F)
P = R*P*R.T

# Neumann condition
N = FacetNormal(mesh)
sigma = Constant(-0.02)
Jg = det(F_glob)
T = dot(Jg*sigma*inv(F_glob.T), N)
ds = Measure("ds")[ff]

# Variational form
L = inner(P, grad(v))*dx - inner(T, v)*ds(1)

# Jacobian
a = derivative(L, u, du)

import time
for i in xrange(3):
	tic()
	A = assemble(a, form_compiler_parameters=compiler_options)
	print 'TIME UF', toc()
import time
ap = inner(grad(du),grad(v))*dx
for i in xrange(3):
	tic()
	A = assemble(ap, form_compiler_parameters=compiler_options)
	print 'TIME UP', toc()
import time
ap = inner(grad(du),grad(v))*dx
compiler_options["representation"] = "tensor"
for i in xrange(3):
	tic()
	A = assemble(ap, form_compiler_parameters=compiler_options)
	print 'TIME TP', toc()
import time
ap = inner(grad(du),grad(v))*dx
compiler_options["representation"] = "quadrature"
compiler_options["optimize"] = True
for i in xrange(3):
	tic()
	A = assemble(ap, form_compiler_parameters=compiler_options)
	print 'TIME QP', toc()
crash




problem = NonlinearVariationalProblem(L, u, bcs, J=a, form_compiler_parameters=compiler_options)
solver = NonlinearVariationalSolver(problem)
solver.parameters["newton_solver"]["maximum_iterations"] = 30

ufile = File("u.pvd")
ufile << u

dt = 0.1
def sigma_t(t):
    return -3.0*(0.1+(sin(t)**2)**0.5)

# Solve problem and plot solution
for i in xrange(1,100):
    print '='*40, "Timestep: ", i
    sigma.assign(sigma_t(i*dt))
    solver.solve()
    ufile << u
    plot(u, mode="displacement")
#interactive()

