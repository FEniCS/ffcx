
#from uflacs.generation.integralgenerator import *
from integralgenerator import *

def debug(title, code):
    print('\n'*3)
    print(':::', title)
    print(code)
    print('\n'*2)

import unittest

class TestIntegralGenerator(unittest.TestCase):
    def test_cell_integral_generator(self):
        g = CellIntegralGenerator()
        code = g.generate()

        debug(g.__class__.__name__, code)

    def xtest_exterior_facet_integral_generator(self):
        g = ExteriorFacetIntegralGenerator()
        code = g.generate()

        debug(g.__class__.__name__, code)

    def xtest_interior_facet_integral_generator(self):
        g = InteriorFacetIntegralGenerator()
        code = g.generate()

        debug(g.__class__.__name__, code)

    def xtest_quadrature_integral_generator(self):
        g = QuadratureIntegralGenerator()
        code = g.generate()

        debug(g.__class__.__name__, code)

    def xtest_diract_integral_generator(self):
        g = DiracIntegralGenerator()
        code = g.generate()

        debug(g.__class__.__name__, code)
