
#from uflacs.generation.integralgenerator import *
from integralgenerator import *

def debug(title, code):
    print '\n'*3
    print ':::', title
    print code
    print '\n'*2

import unittest

class TestIntegralGenerator(unittest.TestCase):
    def test_cell_integral_generator(self):
        g = CellIntegralGenerator()
        code = g.generate()

        debug(g.__class__.__name__, code)

    def test_exterior_facet_integral_generator(self):
        g = ExteriorFacetIntegralGenerator()
        code = g.generate()

        debug(g.__class__.__name__, code)

    def test_interior_facet_integral_generator(self):
        g = InteriorFacetIntegralGenerator()
        code = g.generate()

        debug(g.__class__.__name__, code)

    def test_cutcell_integral_generator(self):
        g = CutcellIntegralGenerator()
        code = g.generate()

        debug(g.__class__.__name__, code)

    def test_subcell_integral_generator(self):
        g = SubcellIntegralGenerator()
        code = g.generate()

        debug(g.__class__.__name__, code)

    def test_diract_integral_generator(self):
        g = DiracIntegralGenerator()
        code = g.generate()

        debug(g.__class__.__name__, code)
