
from six.moves import xrange as range
from ufl.common import product
from ufl.common import component_to_index
from ufl.permutation import build_component_numbering
from ufl.classes import (Coefficient, SpatialCoordinate,
                         FacetNormal, FacetArea,
                         CellVolume, Circumradius)

from uflacs.utils.log import uflacs_assert

class DolfinExpressionStatementFormatter(object):
    def __init__(self, dependency_handler, ir):
        self._dependency_handler = dependency_handler

    def define_registers(self, num_registers):
        return ['double s[%d];' % num_registers]

    def define_piecewise_geometry(self):
        # FIXME: Compute piecewise constant geometry and coeffs in DG0 and R.
        code = []
        for t, c, d, r in self._dependency_handler.terminal_data:
            if isinstance(t, (FacetNormal, FacetArea)):
                code += ["// Facet unknown in Expression eval: %s" % repr(t)]
            elif isinstance(t, (CellVolume, FacetArea, Circumradius)): # TODO: Add other geometry here
                code += ["// Cell unknown in Expression eval: %s" % repr(t)]
        return code

    def define_piecewise_coefficients(self):
        # FIXME: Compute coeffs in DG0 and R.
        evaluated = set()
        code = []
        for t, c, d, r in self._dependency_handler.terminal_data:
            if isinstance(t, Coefficient) and t.is_cellwise_constant():
                uflacs_assert(len(d) == 0, "grad(constant) should be removed at earlier stage by compiler.")
                uflacs_assert(r is None, "Restrictions not supported in dolfin expressions.")

                # Skip already evaluated functions
                if t in evaluated:
                    continue
                evaluated.add(t)

                # Render variable name
                basename = self._dependency_handler.coefficient_names[t]
                #basename = "w%d" % t.count()
                der = "v_"
                res = {"+":0, "-":1, None:""}[r]
                name = "%s%s%s" % (der, basename, res)

                n = product(t.shape())
                tdecl = "Array<double> %s(%d);" % (name, n)

                teval = "%s->eval(%s, x);" % (basename, name)

                code += [tdecl, teval]
        return code

    def define_coord_dependent_coefficients(self):
        # FIXME: Compute x dependent coefficients.
        evaluated = set()
        code = []
        for t, c, d, r in self._dependency_handler.terminal_data:
            if isinstance(t, Coefficient) and not t.is_cellwise_constant():
                uflacs_assert(r is None, "Restrictions not supported in dolfin expressions.")

                # Skip already evaluated functions
                if t in evaluated:
                    continue
                evaluated.add(t)

                # Render variable name
                basename = self._dependency_handler.coefficient_names[t]
                #basename = "w%d" % t.count()
                # FIXME: Handle derivative naming properly
                if d:
                    der = "d%s_" % "".join(derivatives_listing_to_counting(d))
                else:
                    der = "v_"
                res = {"+":0, "-":1, None:""}[r]
                name = "%s%s%s" % (der, basename, res)

                n = product(t.shape())
                tdecl = "Array<double> %s(%d);" % (name, n)

                if not d:
                    teval = "%s->eval(%s, x);" % (basename, name)
                else:
                    teval = "%s->eval_derivatives(%s, x); // FIXME" % (basename, name)

                code += [tdecl, teval]
        return code

    def output_variable_names(self, num_variables):
        return ['%s[%d]' % ("values", i,) for i in range(num_variables)]
