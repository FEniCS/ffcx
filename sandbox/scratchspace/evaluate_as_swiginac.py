
from collections import defaultdict
from itertools import izip, chain

import swiginac

from ufl import *
from ufl.classes import *
from ufl.common import some_key, product, Stack, StackDict
from ufl.algorithms.transformations import Transformer, MultiFunction
from ufl.permutation import compute_indices

from sfc.common import sfc_assert, sfc_error, sfc_warning
from sfc.symbolic_utils import symbol, symbols

class EvaluateAsSwiginac(MultiFunction):
    def __init__(self, formrep, itgrep, data, on_facet): # TODO: Remove unused arguments after implementing:
        MultiFunction.__init__(self)
        self.formrep = formrep
        self.itgrep = itgrep
        self.data = data
        self.on_facet = on_facet
        self.current_basis_function = (None,)*formrep.rank

    ### Fallback handlers:

    def expr(self, o, *ops):
        sfc_error("Evaluation not implemented for expr %s." % type(o).__name__)

    def terminal(self, o, *ops):
        sfc_error("Evaluation not implemented for terminal %s." % type(o).__name__)

    ### Terminals:

    def zero(self, o):
        return swiginac.numeric(0)

    def scalar_value(self, o):
        return swiginac.numeric(o.value())

    def spatial_coordinate(self, o, component=(), derivatives=()):
        # Called by indexed
        
        if component:
            # 2D, 3D
            c, = component
        else:
            # 1D
            c = 0
        
        if derivatives:
            if len(derivatives) > 1:
                return swiginac.numeric(0)
            d, = derivatives
            if d == c:
                return swiginac.numeric(1)
            return swiginac.numeric(0)
        
        return self.formrep.x_sym[c]

    def facet_normal(self, o, component=(), derivatives=()):
        # Called by indexed
        assert self.on_facet, "Expecting to be on a facet in facet_normal."

        if derivatives:
            return swiginac.numeric(0)
        
        if component:
            # 2D, 3D
            c, = component
        else:
            # 1D
            c = 0
        
        return self.formrep.n_sym[c]
    
    def argument(self, o, component=(), derivatives=()):
        
        # Assuming renumbered arguments!
        iarg = o.count()
        sfc_assert(iarg >= 0, "Argument count shouldn't be negative.")
        sfc_assert(isinstance(component, tuple), "Expecting tuple for component.")
        
        j = self.current_basis_function[iarg]
        
        if derivatives:
            s = self.formrep.Dv_sym(iarg, j, component, derivatives, self.on_facet)
            e = self.formrep.Dv_expr(iarg, j, component, derivatives, False, self.on_facet) # FIXME: use_symbols = False -> can do better 
        else:
            s = self.formrep.v_sym(iarg, j, component, self.on_facet)
            e = self.formrep.v_expr(iarg, j, component)
        
        if e.nops() == 0:
            return e # FIXME: Avoid generating code for s when not using it
        return s

    def coefficient(self, o, component=(), derivatives=()):
        # Assuming renumbered arguments!
        iarg = o.count()
        sfc_assert(iarg >= 0, "Coefficient count shouldn't be negative.")
        sfc_assert(isinstance(component, tuple), "Expecting tuple for component.")

        if derivatives:
            s = self.formrep.Dw_sym(iarg, component, derivatives)
            e = self.formrep.Dw_expr(iarg, component, derivatives, False, self.on_facet) # FIXME: use_symbols = False -> can do better
        else:
            # w^i_h(x) = \sum_j w[i][j] * phi^i_j(x)
            s = self.formrep.w_sym(iarg, component)
            e = self.formrep.w_expr(iarg, component, False, self.on_facet) # FIXME: use_symbols = False -> can do better
        
        if e.nops() == 0:
            return e # FIXME: Avoid generating code for s when not using it
        return s

    ### Indexing and derivatives:

    def multi_index(self, o):
        return tuple(map(int, o))

    def spatial_derivative(self, o, f, i, component=(), derivatives=()):
        derivatives = sorted(derivatives + self.multi_index(i))
        if isinstance(f, Indexed):
            # Since expand_indices moves Indexed in to the terminals,
            # SpatialDerivative can be outside an Indexed
            sfc_assert(component == (), "Expecting no component outside of Indexed!")
            A, ii = f.operands()
            component = self.multi_index(ii)
            return self(A, component, derivatives)
        return self(f, component, derivatives)
    
    def indexed(self, o, A, ii):
        # Passes on control to one of:
        #def argument(self, o, component):
        #def coefficient(self, o, component):
        #def facet_normal(self, o, component):
        #def spatial_coordinate(self, o, component):
        #def spatial_derivative(self, o, component):
        component = self.multi_index(ii)
        if isinstance(A, SpatialDerivative):
            f, i = A.operands()
            return self.spatial_derivative(A, f, i, component)
        return self(A, component)

    ### Algebraic operators:

    def power(self, o, a, b):
        return a**b

    def sum(self, o, *ops):
        return sum(ops)

    def product(self, o, *ops):
        return product(ops)

    def division(self, o, a, b):
        return a / b

    def abs(self, o, a):
        return swiginac.abs(a)

    ### Basic math functions:

    def sqrt(self, o, a):
        return swiginac.sqrt(a)

    def exp(self, o, a):
        return swiginac.exp(a)

    def ln(self, o, a):
        return swiginac.log(a)

    def cos(self, o, a):
        return swiginac.cos(a)

    def sin(self, o, a):
        return swiginac.sin(a)

    def variable(self, o):
        sfc_error("Should strip away variables before building graph.")

    # FIXME: Implement all missing operators here
