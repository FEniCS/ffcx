# Copyright (C) 2012 Harish Narayanan

# Library imports and settings
from dolfin import *
from numpy import array, arange
parameters["form_compiler"]["name"] = "sfc"

# Material parameters for Figure 7 in HolzapfelOgden2009
a    =  Constant(0.500)   #kPa
b    =  Constant(8.023)
a_f  =  Constant(16.472)  #kPa
b_f  =  Constant(16.026)
a_s  =  Constant(2.481)   #kPa
b_s  =  Constant(11.120)
a_fs =  Constant(0.356)   #kPa
b_fs =  Constant(11.436)

# Material parameters for compressibility
kappa = Constant(2.0e3)   #kPa
beta  = Constant(9.0)

# Parameters related to time-stepping
T = 10.0
dt = T/100
gamma_max = 0.5

# Parameters related to viscoelasticity
tau = 0.5
beta_inf = 0.25
xi = -dt/(2*tau)

# Strain energy functions for the passive myocardium
# Isochoric part
def psi_iso_inf(I1_bar, I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar):
    return(a/(2*b)*exp(b*(I1_bar - 3)) \
         + a_f/(2*b_f)*(exp(b_f*(I4_f_bar - 1)**2) - 1) \
         + a_s/(2*b_s)*(exp(b_s*(I4_s_bar - 1)**2) - 1) \
         + a_fs/(2*b_fs)*(exp(b_fs*I8_fs_bar**2) - 1))

# Volumetric part
def psi_vol_inf(J):
    return(kappa*(1/(beta**2)*(beta*ln(J) + 1/(J**beta) - 1)))

# Reference fibre, sheet and sheet-normal directions
f0 = Constant((1, 0, 0))
s0 = Constant((0, 1, 0))
n0 = Constant((0, 0, 1))

# Define kinematic measures in terms of the displacement
def kinematics(u):
    I = Identity(u.cell().d)    # Identity tensor
    F = I + grad(u)             # Deformation gradient
    C = F.T*F                   # Right Cauchy-Green tensor
    J = variable(det(F))        # Jacobian
    C_bar = J**(-2.0/3.0)*C     # Modified right Cauchy-Green tensor

    # Principle isotropic invariants
    I1_bar = variable(tr(C_bar))
    I2_bar = variable(0.5*(tr(C_bar)**2 - tr(C_bar*C_bar)))

    # Anisotropic (quasi) invariants
    I4_f_bar = variable(inner(f0, C_bar*f0))
    I4_s_bar = variable(inner(s0, C_bar*s0))
    I8_fs_bar = variable(inner(f0, C_bar*s0))
    I8_fn_bar = variable(inner(f0, C_bar*n0))

    return [I, F, C, J, C_bar, I1_bar, I2_bar, \
            I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar]

# Define the elastic response of the material
# Isochoric part of the second Piola-Kirchhoff stress
def S_iso_inf(u):
    # Define useful kinematic measures
    [I, F, C, J, C_bar, I1_bar, I2_bar, \
     I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar] = kinematics(u)

    # Strain energy functions
    psi_iso = psi_iso_inf(I1_bar, I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar)

    # Define the second Piola-Kirchhoff stress in terms of the invariants
    S_bar =   2*(diff(psi_iso, I1_bar) + diff(psi_iso, I2_bar))*I \
            - 2*diff(psi_iso, I2_bar)*C_bar \
            + 2*diff(psi_iso, I4_f_bar)*outer(f0, f0) \
            + 2*diff(psi_iso, I4_s_bar)*outer(s0, s0) \
            + diff(psi_iso, I8_fs_bar)*(outer(f0, s0) + outer(s0, f0)) \
            + diff(psi_iso, I8_fn_bar)*(outer(f0, n0) + outer(n0, f0))
    Dev_S_bar = S_bar - (1.0/3.0)*inner(S_bar, C)*inv(C)
    S_iso_inf = J**(-2.0/3.0)*Dev_S_bar
    return(S_iso_inf)

# Volumetric part of the second Piola-Kirchhoff stress
def S_vol_inf(u):
    # Define useful kinematic measures
    [I, F, C, J, C_bar, I1_bar, I2_bar, \
     I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar] = kinematics(u)
    psi_vol = psi_vol_inf(J)
    S_vol_inf = J*diff(psi_vol, J)*inv(C)
    return(S_vol_inf)

# Cauchy stress
def sigma(u):
    [I, F, C, J, C_bar, I1_bar, I2_bar, \
     I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar] = kinematics(u)
    return(1/J*P(u)*F.T)

# Dimensions and mesh density of the domain
width = 1
depth = 1
height = 1
n = 10
mesh = Box(0, width, 0, depth, 0, height, n*width, n*depth, n*height)

# Function spaces
scalar = FunctionSpace(mesh, "Lagrange", 1)
vector = VectorFunctionSpace(mesh, "Lagrange", 1)
tensor = TensorFunctionSpace(mesh, "Lagrange", 1)

# Functions
du = TrialFunction(vector)            # Incremental displacement
v  = TestFunction(vector)             # Test function
u  = Function(vector)                 # Displacement from previous iteration
S_iso_inf_p = Function(tensor)
S_vol_inf_p = Function(tensor)
Q_p = Function(tensor)

# Boundary conditions
back_condition   = "x[0] == 0.0 && on_boundary"
front_condition  = "x[0] == %g && on_boundary" % depth
left_condition   = "x[1] == 0.0 && on_boundary"
right_condition  = "x[1] == %g && on_boundary" % width
bottom_condition = "x[2] == 0.0 && on_boundary"
top_condition    = "x[2] == %g && on_boundary" % height

back, front = compile_subdomains([back_condition, front_condition])
left, right = compile_subdomains([left_condition, right_condition])
bottom, top = compile_subdomains([bottom_condition, top_condition])

hold = Expression(("0.0", "0.0", "0.0"))

# Simple shear along the fs plane
shear = Expression(("0.0", "gamma*depth", "0.0"), gamma=0.0, depth=depth)
hold_back = DirichletBC(vector, hold, back)
shear_front = DirichletBC(vector, shear, front)
bcs = [hold_back, shear_front]

# Create files to store output
u_store = TimeSeries("../output/visco/u")
S_vol_inf_store = TimeSeries("../output/visco/S_vol_inf")
S_iso_inf_store = TimeSeries("../output/visco/S_iso_inf")
Q_store = TimeSeries("../output/visco/Q")

# And store initial values
u_store.store(u.vector(), 0.0)
S_iso_inf_store.store(S_iso_inf_p.vector(), 0.0)
S_vol_inf_store.store(S_vol_inf_p.vector(), 0.0)
Q_store.store(Q_p.vector(), 0.0)

# Define the time range
times = arange(dt, T + dt, dt)

# Subject the body to a known strain protocol and record the stresses
for t_n in times:
    print "t_n = ", t_n

    # Load previous elastic stress states
    t_p = t_n - dt
    S_vol_inf_store.retrieve(S_vol_inf_p.vector(), t_p)
    S_iso_inf_store.retrieve(S_iso_inf_p.vector(), t_p)
    Q_store.retrieve(Q_p.vector(), t_p)

    # Compute current shear strain and update the boundary condition
    gamma_n = gamma_max*sin(2*t_n/T*float(pi))
    shear.gamma = gamma_n

    # Update stress state
    S_vol_inf_n = S_vol_inf(u)
    S_iso_inf_n = S_iso_inf(u)
    H_p = exp(xi)*(exp(xi)*Q_p - beta_inf*S_iso_inf_p)
    Q_n = beta_inf*exp(xi)*S_iso_inf_n + H_p
    S_n = S_vol_inf_n + S_iso_inf_n + Q_n

    # Define the variational form for the problem
    F = inner((Identity(u.cell().d) + grad(u))*S_n, grad(v))*dx
    J = derivative(F, u, du)

    # Solve the boundary value problem
    solve(F == 0, u, bcs, J=J)

    # Project the stress states to the tensor function space
    S_iso_inf_n = project(S_iso_inf(u), tensor)
    S_vol_inf_n = project(S_vol_inf(u), tensor)
    Q_n = project(Q_n, tensor)

    # Store the displacement and stress state at current time
    u_store.store(u.vector(), t_n)
    S_iso_inf_store.store(S_iso_inf_n.vector(), t_n)
    S_vol_inf_store.store(S_vol_inf_n.vector(), t_n)
    Q_store.store(Q_n.vector(), t_n)

    # Convert to Cauchy stress for comparison with Dokos et al.
#    S_n = F.subs({gamma:gamma_n})*S_n*F.T.subs({gamma:gamma_n})
