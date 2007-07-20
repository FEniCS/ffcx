# This script calls the just-in-time compiler of FFC
# which generates code, compiles the code, wraps it
# back to Python (using Instant/SWIG) and returns it
# as a Python object

import sys
sys.path.append("../../")

from ffc import *

element = FiniteElement("Lagrange", "triangle", 1)

v = TestFunction(element)
u = TrialFunction(element)

a = dot(grad(v), grad(u))*dx

compiled_form = jit(a)

print compiled_form.rank(), compiled_form.num_coefficients()
