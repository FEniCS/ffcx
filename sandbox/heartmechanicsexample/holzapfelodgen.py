from dolfin import *

class HolzapfelOgden(object) :

  # initialization
  def __init__(self, f0, s0, mattype, elptfix = True) :

    self.mattype = mattype
    self.elptfix = elptfix
    
    # fibers and sheets
    self.f0 = f0
    self.s0 = s0
    
    # parameters [N/cm^2]
    if (mattype == "isotropic") :
      self.a = 0.44
      self.b = 9.726
      self.a_f  = 1.0
      self.b_f  = 0.0
      self.a_s  = 0.0
      self.b_s  = 1.0
      self.a_fs = 0.0
      self.b_fs = 1.0
    elif (mattype == "tr_isotropic" ) :
      self.a = 2.280e-1
      self.b = 9.726
      self.a_f  = 1.685e-1
      self.b_f  = 15.779
      self.a_s  = 0.0
      self.b_s  = 1.0
      self.a_fs = 0.0
      self.b_fs = 1.0
    elif (mattype == "goktepe") :
      self.a = 0.496e-1
      self.b = 7.209
      self.a_f  = 15.193e-1
      self.b_f  = 20.417
      self.a_s  = 3.283e-1
      self.b_s  = 11.176
      self.a_fs = 0.662e-1
      self.b_fs = 9.466
    elif (mattype == "orthotropic") :
      self.a = 0.059e-1
      self.b = 8.029
      self.a_f  = 18.472e-1
      self.b_f  = 16.026
      self.a_s  = 2.481e-1
      self.b_s  = 11.120
      self.a_fs = 0.216e-1
      self.b_fs = 11.436
    else :
      raise Exception("Material not defined!")

  def I_1(self, F) :

    return inner(F, F)

  def I_2(self, F) :

    C = F.T * F
    return 0.5*(I_1(F)**2 - tr(C*C))

  def I_4f(self, F) :
    
    f = F * self.f0
    return inner(f, f)
  
  def I_5f(self, F) :

    f0 = self.f0
    C = F.T * F
    return inner(C * f0, C * f0)

  def I_4s(self, F) :

    s = F *self.s0
    return inner(s, s)

  def I_5s(self, F) :

    s0 = self.s0
    C = F.T * F
    return inner(C * s0, C * s0)

  def I_8fs(self, F) :

    f = F * self.f0
    s = F * self.s0
    return 0.5*(inner(s, f) + inner(f, s))

  # total strain-energy function
  def W(self, F) :
    
    # directions
    f0 = self.f0
    s0 = self.s0
    # invariants
    I_1  = self.I_1(F)
    I_f  = self.I_4f(F)
    I_s  = self.I_4s(F)
    I_fs = self.I_8fs(F)
    # constants
    k_1  = Constant(0.5 * self.a / self.b)
    k_f  = Constant(0.5 * self.a_f / self.b_f)
    k_s  = Constant(0.5 * self.a_s / self.b_s)
    k_fs = Constant(0.5 * self.a_fs / self.b_fs)
    b_1  = Constant(self.b)
    b_f  = Constant(self.b_f)
    b_s  = Constant(self.b_s)
    b_fs = Constant(self.b_fs)
    # energies
    W_1  = k_1  * (exp(b_1  * (I_1 - 3)) - 1.0)
    W_f  = k_f  * (exp(b_f  * (I_f - 1)**2) - 1.0)
    W_s  = k_s  * (exp(b_s  * (I_s - 1)**2) - 1.0)
    W_fs = k_fs * (exp(b_fs * I_fs**2) - 1.0)
    # ellipticity fix
    if (self.elptfix == True) :
      W_f = conditional(ge(I_f, 1.0), W_f, Constant(0.0))
      W_s = conditional(ge(I_s, 1.0), W_s, Constant(0.0))
    # total energy
    if (self.mattype == "isotropic") :
      W = W_1
    elif (self.mattype == "tr_isotropic") :
      W = W_1 + W_f
    elif (self.mattype == "orthotropic" or self.mattype == "goktepe") :
      W = W_1 + W_f + W_s + W_fs
    
    return W

  # total second piola-kirchhoff
  def S(self, F) :
    # directions
    f0 = self.f0
    s0 = self.s0
    # invariants
    I_1   = self.I_1(F)
    I_4f  = self.I_4f(F)
    I_4s  = self.I_4s(F)
    I_8fs = self.I_8fs(F)
    # constants
    k_1  = Constant(0.5 * self.a)
    k_f  = Constant(self.a_f)
    k_s  = Constant(self.a_s)
    k_fs = Constant(self.a_fs)
    b_1  = Constant(self.b)
    b_f  = Constant(self.b_f)
    b_s  = Constant(self.b_s)
    b_fs = Constant(self.b_fs)
    # derivatives
    dW_1  = k_1  * exp(b_1*(I_1-3))
    dW_f  = k_f  * (I_4f-1)  * exp(b_f*(I_4f-1)**2)
    dW_s  = k_s  * (I_4s-1)  * exp(b_s*(I_4s-1)**2)
    dW_fs = k_fs * I_8fs * exp(b_fs*I_8fs**2)
    # ellipticity fix
    if (self.elptfix == True) :
      dW_f = conditional(ge(I_4f, 1.0), dW_f, 0.0)
      dW_s = conditional(ge(I_4s, 1.0), dW_s, 0.0)
    # total stress
    if (self.mattype == "isotropic") :
      dW_dC = dW_1 * Identity(3)
    elif (self.mattype == "tr_isotropic") :
      dW_dC = dW_1 * Identity(3) \
            + dW_f * outer(f0, f0)
    elif (self.mattype == "orthotropic" or self.mattype == "goktepe") :
      dW_dC = dW_1 * Identity(3) \
            + dW_f * outer(f0, f0) \
            + dW_s * outer(s0, s0) \
            + dW_fs * sym(outer(f0, s0))
    
    return 2.0 * dW_dC
  
  # total cauchy stress
  def T(self, F) :

    # deformed directions
    f = F * self.f0
    s = F * self.s0
    # invariants
    I_1  = self.I_1(F)
    I_f  = self.I_4f(F)
    I_s  = self.I_4s(F)
    I_fs = self.I_8fs(F)
    # tensors
    I  = Identity(3)
    B  = F * F.T
    ff = outer(f, f)
    ss = outer(s, s)
    fs = sym(outer(f, s))
    # constants
    k_1  = Constant(self.a)
    k_f  = Constant(2.0 * self.a_f)
    k_s  = Constant(2.0 * self.a_s)
    k_fs = Constant(2.0 * self.a_fs)
    b_1  = Constant(self.b)
    b_f  = Constant(self.b_f)
    b_s  = Constant(self.b_s)
    b_fs = Constant(self.b_fs)
    # stresses
    T_1  = k_1  * exp(b_1*(I_1-3))
    T_f  = k_f  * (I_f-1) * exp(b_f*(I_f-1)**2)
    T_s  = k_s  * (I_s-1) * exp(b_s*(I_s-1)**2)
    T_fs = k_fs * I_fs * exp(b_fs*I_fs**2)
    # ellipticity fix
    if (self.elptfix == True) :
      T_f = conditional(ge(I_f, 1.0), T_f, 0.0)
      T_s = conditional(ge(I_s, 1.0), T_s, 0.0)
    # total stress
    if (self.mattype == "isotropic") :
      T = T_1 * B
    if (self.mattype == "tr_isotropic") :
      T = T_1 * B + T_f * ff
    if (self.mattype == "orthotropic" or self.mattype == "goktepe") :
      T = T_1 * B + T_s * ss + T_fs * fs

    return T
    