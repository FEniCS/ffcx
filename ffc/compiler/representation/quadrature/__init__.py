from quadraturerepresentation import *
try:
    from uflquadraturerepresentation import QuadratureRepresentation as UFLQuadratureRepresentation
except:
    warning("UFL quadrature representation module did not load correctly")
    pass
