
import uflacs.backends.ffc as uffc
from ufl import *

cell = triangle

x = cell.x[0]
xi = cell.xi[0]
J = cell.J[0,0]
detJ = cell.detJ
K = cell.Jinv[0,0]

a = (x*xi)*(detJ*J*K)*dx

forms = [a]

for form in forms:
    fd = a.compute_form_data()
    ffc_data = None
    code = uffc.compile_tabulate_tensor_code(ffc_data, fd.integral_data[0], fd, {})
    print
    print '/'*60
    print code
    print

