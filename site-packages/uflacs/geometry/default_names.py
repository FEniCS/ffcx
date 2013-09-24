
class Names:
    def __init__(self):
        # Topology argument names
        self.vertex = "vertex"
        self.facet = "facet"

        # Geometry names
        self.vertex_coordinates = "vertex_coordinates"
        self.xi = "xi"
        self.x = "x"
        self.J = "J"
        self.K = "K"
        self.detJ = "detJ"
        self.det = "det"

        # Quadrature rule
        self.points = "points"
        self.weights = "weights"

        # Quadrature temps
        self.qw = "qw"
        self.D = "D"

        # (Base)name for intermediate registers
        self.s = "s"

        # Element tensor
        self.A = "A"

        # Coefficient dofs array
        self.w = "w"

        # Basenames for function components
        self.wbase = "w"
        self.vbase = "v"
        self.dwbase = "dw"
        self.dvbase = "dv"

        # Loop indices
        self.iq = "iq"   # Quadrature loop
        self.ic = "ic"   # Coefficient accumulation loop
        self.ia = "ia"   # Argument dof loop
        self.ild = "ild" # Local derivative accumulation loop

        # Rules, make functions?
        self.restriction_postfix = { "+": "_0", "-": "_1", None: "" } # TODO: Use this wherever we need it?

names = Names()
