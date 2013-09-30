
from uflacs.codeutils.format_code_structure import format_code_structure, Indented, Block, ForRange
from uflacs.geometry.default_names import names


class IntegralGenerator(object):
    def __init__(self):
        # TODO: This is just fake test data

        self._argument_space_dimensions = [7, 8] # ir["prim_idims"] # FIXME: *2 for dS?
        rank = len(self._argument_space_dimensions)
        self._idofs = ["%s%d" % (names.ia, i) for i in range(rank)]


        self._num_points = [1, 2]
        self._dof_ranges = {
            1: [
                ((0,3), (3,6)),
                ((1,2), (4,5))
                ],
            2: [
                ((1,2), (4,6)),
                ((0,5),)
                ],
            }
        self._num_arguments = 2

    def generate(self):
        """Generate entire tabulate_tensor body."""
        parts = []
        parts += [self.generate_pre_quadrature()]
        parts += [self.generate_quadrature_loops()]
        parts += [self.generate_post_quadrature()]
        return format_code_structure(Indented(parts))

    def generate_pre_quadrature(self): # FIXME
        """Generate initialization statements.

        This includes initial piecewise constant geometry computations and such.
        """
        parts = ["// Piecewise constant stage"]
        parts += ["// FIXME: Implement this"]
        return parts

    def generate_post_quadrature(self): # TODO
        """Generate finishing statements.

        This includes assigning to output array if there is no integration.
        """
        parts = ["// Finishing statements"]
        parts += ["// TODO: Implement this"]
        return parts

    def generate_quadrature_loops(self):
        """Generate all quadrature loops."""
        parts = [self.generate_quadrature_loop(num_points) for num_points in self._num_points]
        return parts

    def generate_quadrature_loop(self, num_points):
        """Generate single quadrature loop."""
        return ForRange(names.iq, 0, num_points,
                        body=self.generate_quadrature_body(num_points))

    def generate_quadrature_body_setup(self, num_points): # FIXME
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

        if self._num_arguments == 0: # Functional
            # Accumulate into element scalar
            dofblock = ()
            parts += [self.generate_integrand_accumulation(num_points, dofblock)]

        elif self._num_arguments == 1: # Linear form
            # Loop over dofranges of argument 0 and accumulate into element vector
            for dofrange0 in drs[0]:
                (b0,e0) = dofrange0
                dofblock = (dofrange0,)

                body0 = [self.generate_argument_partition(num_points, 0, dofblock)]
                body0 += [self.generate_integrand_accumulation(num_points, dofblock)]
                loop0 = ForRange(self._idofs[0], b0, e0, body=body0)

                parts += [loop0]

        elif self._num_arguments == 2: # Bilinear form
            # Loop over dofranges of argument 0 (rows of element matrix)
            for dofrange0 in drs[0]:
                (b0,e0) = dofrange0
                dofblock = (dofrange0,)

                body0 = [self.generate_argument_partition(num_points, 0, dofblock)]

                # Loop over dofranges of argument 1 and accumulate into element matrix
                for dofrange1 in drs[1]:
                    (b1,e1) = dofrange1
                    dofblock = (dofrange0,dofrange1)

                    body1 = [self.generate_argument_partition(num_points, 1, dofblock)]
                    body1 += [self.generate_integrand_accumulation(num_points, dofblock)]
                    loop1 = ForRange(self._idofs[1], b1, e1, body=body1)

                    body0 += [loop1]

                loop0 = ForRange(self._idofs[0], b0, e0, body=body0)

                parts += [loop0]
        else:
            error("Bonus points: Generalize this recursively to support num_arguments > 2.")

        return parts

    def generate_argument_partition(self, num_points, iarg, dofblock): # FIXME
        # FIXME: Get partition associated with (num_points, iarg, dofblock)
        # FIXME: Generate accumulation properly
        parts = []
        parts += ["s[...] = ...; // {0} x {1}".format(iarg, dofblock)]
        return parts

    def generate_integrand_accumulation(self, num_points, dofblock): # FIXME
        # FIXME: Get partition associated with (num_points, dofblock)
        # FIXME: Generate accumulation properly
        parts = []
        parts += ["A[{0}] += f * v * D;".format(dofblock)]
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
