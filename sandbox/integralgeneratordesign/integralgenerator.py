
from uflacs.codeutils.format_code_structure import format_code_structure, Indented, Block


class IntegralGenerator(object):
    def __init__(self):
        # TODO: This is just fake test data
        self._num_points = [1, 2]
        self._dof_ranges = {
            1: [
                ((0,3), (3,6)),
                ((1,3), (2,6))
                ],
            2: [
                ((0,2), (4,6)),
                ((0,5),)
                ],
            }
        self._num_arguments = 2

    def generate(self):
        """Generate entire tabulate_tensor body."""
        parts = []
        parts += [self.generate_body_setup()]
        parts += [self.generate_body_parts()]
        return format_code_structure(Indented(parts))

    def generate_body_setup(self):
        """Generate setup stage statements.

        This includes initial piecewise constant geometry computations and such.
        """
        parts = ["// Piecewise constant stage"]
        parts += ["// TODO: Implement this"]
        return parts

    def generate_body_parts(self):
        """Generate all quadrature loops."""
        parts = ["// Begin quadrature loops"]
        for num_points in self._num_points:
            parts += [self.generate_quadrature_loop(num_points)]
        parts += ["// End quadrature loops"]
        return parts

    def generate_quadrature_loop(self, num_points):
        """Generate single quadrature loop."""
        parts = ["// Quadrature loop {0}".format(num_points)]
        parts += ['for (iq{0}=0; iq{0}<nq{0}; ++iq{0})'.format(num_points)] # FIXME
        parts += [Block(self.generate_quadrature_body(num_points))]
        return parts

    def generate_quadrature_body_setup(self, num_points):
        """
        """
        parts = ["// Quadrature loop body setup {0}".format(num_points)]
        # FIXME: All leftover argument independent computations here
        return parts

    def generate_quadrature_body(self, num_points):
        """
        """
        parts = []
        parts += [self.generate_quadrature_body_setup(num_points)]

        drs = self._dof_ranges[num_points]

        if self._num_arguments == 0:
            dofblock = []
            parts += [self.generate_integrand_accumulation(num_points, dofblock)]

        elif self._num_arguments == 1:
            for dofrange0 in drs[0]:
                (b0,e0) = dofrange0
                dofblock = [dofrange0]
                parts += ["for ({ia} = {begin}; {ia} < {end}; ++{ia})".format(ia="ia0", begin=b0, end=e0)]
                body0 = [
                    self.generate_argument_partition(num_points, 0, dofrange0),
                    self.generate_integrand_accumulation(num_points, dofblock),
                    ]
                parts += Block(body0)

        elif self._num_arguments == 2:
            for dofrange0 in drs[0]:
                (b0,e0) = dofrange0
                loop0 = ["for ({ia} = {begin}; {ia} < {end}; ++{ia})".format(ia="ia0", begin=b0, end=e0)]
                body0 = [self.generate_argument_partition(num_points, 0, dofrange0)]
                for dofrange1 in drs[1]:
                    (b1,e1) = dofrange1
                    dofblock = [dofrange0, dofrange1]
                    loop1 = ["for ({ia} = {begin}; {ia} < {end}; ++{ia})".format(ia="ia1", begin=b1, end=e1)]
                    body1 = [self.generate_argument_partition(num_points, 1, dofrange1),
                             self.generate_integrand_accumulation(num_points, dofblock)]
                    body0 += [loop1, Block(body1)]
                parts += [loop0, Block(body0)]

        # TODO: Generalize this recursively to support num_arguments > 2

        return parts

    def generate_argument_partition(self, num_points, iarg, dofrange):
        parts = []
        parts += ["s[...] = ...; // {0} x {1}".format(iarg, dofrange)]
        return parts

    def generate_integrand_accumulation(self, num_points, dofblock):
        parts = []
        parts += ["A[{0}] = ...;".format(dofblock)]
        return parts


class CellIntegralGenerator(IntegralGenerator):
    pass

class ExteriorFacetIntegralGenerator(IntegralGenerator):
    pass

class InteriorFacetIntegralGenerator(IntegralGenerator):
    pass

class DiracIntegralGenerator(IntegralGenerator):
    pass

class SubcellIntegralGenerator(IntegralGenerator):
    pass

class CutcellIntegralGenerator(IntegralGenerator):
    pass
