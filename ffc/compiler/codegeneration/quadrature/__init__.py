from quadraturegenerator import *
try:
    from uflquadraturegenerator import QuadratureGenerator as UFLQuadratureGenerator
except:
    warning("UFL quadrature representation module did not load correctly")
    pass
