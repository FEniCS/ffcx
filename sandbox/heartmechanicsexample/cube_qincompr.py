# :: PROBLEM ::
# uniaxial and biaxial stretch of a cube
# Holzapfel-Ogden material
# quasi-incompressibility
# zero mean displacement and rotation

from dolfin import *
from holzapfelodgen import *
import numpy as np

# ::::::::::::::::
# :: parameters ::
# ::::::::::::::::
# output directory
batchdir = "cube_tr_fix"
# mesh data
length = 1.0 # [cm]
subdiv = 2
# microstructure
f0 = Constant((1.0, 0.0, 0.0))
s0 = Constant((0.0, 1.0, 0.0))
# traction at the boundaries (symmetric)
fx = Constant(0.1) # [N/cm^2]
# material parameters
mattype = "tr_isotropic"
elptfix = False
# form parameters
_ffc_options = \
{
  "quadrature_degree" : 2,
  "optimize" : True,
  "eliminates_zeros" : True,
  "precompute_basis_const" : True,
  "cpp_optimize" : True
}
ffc_options = \
{
  "representation" : "uflacs",
  "quadrature_degree" : 2,
}
# solver parameters
solver_options = \
{
  "linear_solver" : "mumps",
  "newton_solver" : \
  {
    "maximum_iterations" : 50
  }
}

# :::::::::::
# :: utils ::
# :::::::::::
# creates mesh and boundaries
def create_mesh(L, N) :
  mesh = BoxMesh(-L/2, -L/2, -L/2, L/2, L/2, L/2, N, N, N)
  left, right, bottom, top, back, front = \
  compile_subdomains([
    "std::abs(x[0]+0.5) < DOLFIN_EPS && on_boundary",
    "std::abs(x[0]-0.5) < DOLFIN_EPS && on_boundary",
    "std::abs(x[1]+0.5) < DOLFIN_EPS && on_boundary",
    "std::abs(x[1]-0.5) < DOLFIN_EPS && on_boundary",
    "std::abs(x[2]+0.5) < DOLFIN_EPS && on_boundary",
    "std::abs(x[2]-0.5) < DOLFIN_EPS && on_boundary"])
  bfun = FacetFunctionSizet(mesh, 0)
  left.mark(bfun, 10)
  right.mark(bfun, 20)
  back.mark(bfun, 30)
  front.mark(bfun, 40)
  bottom.mark(bfun, 50)
  top.mark(bfun, 60)
  
  return mesh, bfun

# ::::::::::
# :: main ::
# ::::::::::
if __name__ == "__main__" :
  
  # geometry
  mesh, bfun = create_mesh(1.0, 32)
  dGamma = Measure("ds")[bfun]
  
  # function space
  #V = VectorFunctionSpace(mesh, "CG", 2)
  P1 = VectorFunctionSpace(mesh, "CG", 1)
  #B  = VectorFunctionSpace(mesh, "Bubble", 4)
  V = P1
  R = VectorFunctionSpace(mesh, "R", 0)
  W = MixedFunctionSpace([V, R, R])
  
  # variational problem
  w = Function(W)
  u, c1, c2 = split(w)
  v, d1, d2 = TestFunctions(W)
  
  I = Identity(3)
  F = I + grad(u)
  F = variable(F)
  J = det(F)
  C = F.T * F
  Fiso = pow(J, -float(1)/3) * F
   
  # 2nd Piola
  material = HolzapfelOgden(f0, s0, mattype, elptfix)
  #S = material.S(F) - p * inv(C)
  
  kappa = 350.e-1
  Pi = material.W(Fiso) + Constant(kappa/4.)*(J**2 - J - 2*ln(J))
  
  P = diff(Pi, F)
  
  N = FacetNormal(mesh)
  X = SpatialCoordinate(W.cell())
  
  # internal stress
  G  = inner(grad(v), P) * dx
  # external traction
  G -= fx * inner(cofac(F)*N, v) * dGamma(10)
  G -= fx * inner(cofac(F)*N, v) * dGamma(20)
  # mean displacement zero
  G += inner(v, c1) * dx + inner(d1, u) * dx
  # mean rotation zero
  G += inner(cross(X, v), c2) * dx + inner(d2, cross(X, u)) * dx
  
  # tangent problem
  dw = TrialFunction(W)
  dG = derivative(G, w, dw)

  # evaluation points of the strain
  X1 = np.array([length/2, 0.0, 0.0])
  X2 = np.array([0.0, length/2, 0.0])
  X3 = np.array([0.0, 0.0, length/2])
  
  # solve
  solve(G == 0, w, J = dG, \
    form_compiler_parameters = ffc_options, \
    solver_parameters = solver_options)
 
  # export the solution
  u_file = File(batchdir+"/u.pvd", "compressed")
  I_file = File(batchdir+"/I.pvd", "compressed")
  T_file = File(batchdir+"/T.pvd", "compressed")

  u, _, _ = w.split()
  u_file << u
 
  # new position
  x1 = X1 + u(X1)
  x2 = X2 + u(X2)
  x3 = X3 + u(X3)
  # strain F
  l1 = 2.0/length * x1[0]
  l2 = 2.0/length * x2[1]
  l3 = 2.0/length * x3[2]
  # strain E
  e1 = (l1**2 - 1)/2.0
  e2 = (l2**2 - 1)/2.0
  e3 = (l3**2 - 1)/2.0
  
  print("F = diag(%e, %e, %e)" % (l1, l2, l3))
  print("E = diag(%e, %e, %e)" % (e1, e2, e3))

  # invariants
  I_f  = material.I_4f(F)
  I_s  = material.I_4s(F)
  I_fs = material.I_8fs(F)

  ## cauchy stress
  T = material.T(F)
  ff = F*f0
  ss = F*s0
  
  T_f  = inner(T*ff, ff) / dot(ff, ff)
  T_s  = inner(T*ss, ss) / dot(ss, ss)
  T_fs = inner(T*ff, ss) / sqrt(dot(ff, ff)*dot(ss, ss))

  ## export
  Xh = VectorFunctionSpace(mesh, "DG", 0)
  Tset = project(as_vector([T_f, T_s, T_fs]), Xh, \
    form_compiler_parameters = ffc_options)
  Iset = project(as_vector([I_f, I_s, I_fs]), Xh, \
    form_compiler_parameters = ffc_options)

  u_file << u
  I_file << Iset
  T_file << Tset
