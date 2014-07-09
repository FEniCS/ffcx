
from six.moves import xrange, zip
from collections import defaultdict
from itertools import chain

import swiginac

from ufl import *
from ufl.classes import *
from ufl.common import some_key, product, Stack, StackDict
from ufl.algorithms.transformations import Transformer, MultiFunction
from ufl.permutation import compute_indices

from sfc.common import sfc_assert, sfc_error, sfc_warning
from sfc.symbolic_utils import symbol, symbols

class SwiginacEvaluator(Transformer):
    "Algorithm for evaluation of an UFL expression as a swiginac expression."
    def __init__(self, formrep, use_symbols, on_facet):
        Transformer.__init__(self)#, variable_cache)

        # input
        self._formrep = formrep
        self._use_symbols = use_symbols
        self._on_facet = on_facet

        # current basis function configuration
        self._current_basis_function = tuple(0 for i in xrange(formrep.rank))

        # current indexing status
        self._components = Stack()
        self._index2value = StackDict()

        # code and cache structures
        # FIXME: Need pre-initialized self._variable2symbol and self._tokens
        self._variable2symbol = {}
        self._tokens = []

        # convenience variables
        self.nsd = self._formrep.cell.nsd

    def pop_tokens(self):
        # TODO: Make a generator approach to this? Allow handlers to trigger a "token yield"?
        t = self._tokens
        self._tokens = []
        return t

    def update(self, iota):
        self._current_basis_function = tuple(iota)

    def component(self):
        "Return current component tuple."
        if len(self._components):
            return self._components.peek()
        return ()

    ### Fallback handlers:

    def expr(self, x):
        sfc_error("Missing ufl to swiginac handler for type %s" % str(type(x)))

    def terminal(self, x):
        sfc_error("Missing ufl to swiginac handler for terminal type %s" % str(type(x)))

    ### Handlers for basic terminal objects:

    def zero(self, x):
        #sfc_assert(len(self.component()) == len(x.shape()), "Index component length mismatch in zero tensor!")
        return swiginac.numeric(0)

    def scalar_value(self, x):
        #sfc_assert(self.component() == (), "Shouldn't have any component at this point.")
        return swiginac.numeric(x._value)

    def identity(self, x):
        c = self.component()
        v = 1 if c[0] == c[1] else 0
        return swiginac.numeric(v)

    def argument(self, x):
        iarg = x.count()
        sfc_assert(iarg >= 0, "Argument count shouldn't be negative.")
        j = self._current_basis_function[iarg]
        c = self.component()
        if self._use_symbols:
            return self._formrep.v_sym(iarg, j, c, self._on_facet)
        else:
            return self._formrep.v_expr(iarg, j, c)

    def coefficient(self, x):
        iarg = x.count()
        c = self.component()
        if self._use_symbols:
            return self._formrep.w_sym(iarg, c)
        else:
            # w^i_h(x) = \sum_j w[i][j] * phi^i_j(x)
            return self._formrep.w_expr(iarg, c, False, self._on_facet)

    def facet_normal(self, x):
        sfc_assert(self._on_facet, "Expecting to be on a facet in facet_normal.")
        c, = self.component()
        return self._formrep.n_sym[c]

    def spatial_coordinate(self, x):
        c, = self.component()
        return self._formrep.x_sym[c]

    ### Handler for variables:

    def variable(self, x):
        return self.visit(x._expression)

    def garbage(self, x): # TODO: Maybe some of the ideas here can be used in code generation
        c = self.component()
        index_values = tuple(self._index2value[k] for k in x._expression.free_indices())

        # TODO: Doesn't always depend on _current_basis_function, this is crap:
        key = (x.count(), c, index_values, self._current_basis_function)
        vsym = self._variable2symbol.get(key)

        if vsym is None:
            expr = self.visit(x._expression)
            # TODO: Doesn't always depend on _current_basis_function, this is crap:
            compstr = "_".join("%d" % k for k in chain(c, index_values, self._current_basis_function))
            vname = "_".join(("t_%d" % x.count(), compstr))
            vsym = symbol(vname)
            self._variable2symbol[key] = vsym
            self._tokens.append((vsym, expr))

    ### Handlers for basic algebra:

    def sum(self, x, *ops):
        return sum(ops)

    def index_sum(self, x):
        ops = []
        summand, multiindex = x.operands()
        index, = multiindex
        for i in xrange(x.dimension()):
            self._index2value.push(index, i)
            ops.append(self.visit(summand))
            self._index2value.pop()
        return sum(ops)

    def product(self, x):
        sfc_assert(not self.component(), "Non-empty indexing component in product!")
        ops = [self.visit(o) for o in x.operands()]
        return product(ops)
        # ...

    def division(self, x, a, b):
        return a / b

    def power(self, x, a, b):
        return a ** b

    def abs(self, x, a):
        return swiginac.abs(a)

    ### Basic math functions:
    def sqrt(self, x, y):
        return swiginac.sqrt(y)

    def exp(self, x, y):
        return swiginac.exp(y)

    def ln(self, x, y):
        return swiginac.log(y)

    def cos(self, x, y):
        return swiginac.cos(y)

    def sin(self, x, y):
        return swiginac.sin(y)

    ### Index handling:
    def multi_index(self, x):
        subcomp = []
        for i in x:
            if isinstance(i, FixedIndex):
                subcomp.append(i._value)
            elif isinstance(i, Index):
                subcomp.append(self._index2value[i])
        return tuple(subcomp)

    def indexed(self, x):
        A, ii = x.operands()
        self._components.push(self.visit(ii))
        result = self.visit(A)
        self._components.pop()
        return result

    ### Container handling:

    def old_list_tensor(self, x): # doesn't support e.g. building a matrix from vector rows
        component = self.component()
        sfc_assert(len(component) > 0 and \
                   all(isinstance(i, int) for i in component),
                   "Can't index tensor with %s." % repr(component))

        # Hide indexing when evaluating subexpression
        self._components.push(())

        # Get scalar UFL subexpression from tensor
        e = x
        for i in component:
            e = e._expressions[i]
        sfc_assert(e.shape() == (), "Expecting scalar expression "\
                   "after extracting component from tensor.")

        # Apply conversion to scalar subexpression
        r = self.visit(e)

        # Return to previous component state
        self._components.pop()
        return r

    def list_tensor(self, x):
        # Pick the right subtensor and subcomponent
        c = self.component()
        c0, c1 = c[0], c[1:]
        op = x.operands()[c0]
        # Evaluate subtensor with this subcomponent
        self._components.push(c1)
        r = self.visit(op)
        self._components.pop()
        return r

    def component_tensor(self, x):
        # this function evaluates the tensor expression
        # with indices equal to the current component tuple
        expression, indices = x.operands()
        sfc_assert(expression.shape() == (), "Expecting scalar base expression.")

        # update index map with component tuple values
        comp = self.component()
        sfc_assert(len(indices) == len(comp), "Index/component mismatch.")
        for i, v in zip(indices._indices, comp):
            self._index2value.push(i, v)
        self._components.push(())

        # evaluate with these indices
        result = self.visit(expression)

        # revert index map
        for i in xrange(len(comp)):
            self._index2value.pop()
        self._components.pop()
        return result

    ### Differentiation:

    def _ddx(self, f, i):
        """Differentiate swiginac expression f w.r.t. x_i, using
        df/dx_i = df/dxi_j dxi_j/dx_i."""
        Ginv = self._formrep.Ginv_sym
        xi = self._formrep.xi_sym
        return sum(Ginv[j, i] * swiginac.diff(f, xi[j]) for j in xrange(self.nsd))

    def spatial_derivative(self, x):
        # Assuming that AD has been applied, so
        # the expression to differentiate is always a Terminal.

        f, ii = x.operands()

        sfc_assert(isinstance(f, Terminal), \
            "Expecting to differentiate a Terminal object, you must apply AD first!") # The exception is higher order derivatives, ignoring for now

        # Get component and derivative directions
        c = self.component()
        der = self.visit(ii)

        # --- Handle derivatives of basis functions
        if isinstance(f, Argument):
            iarg = f.count()
            i = self._current_basis_function[iarg]
            if self._use_symbols:
                return self._formrep.Dv_sym(iarg, i, c, der, self._on_facet)
            else:
                return self._formrep.Dv_expr(iarg, i, c, der, False, self._on_facet)

        # --- Handle derivatives of coefficient functions
        if isinstance(f, Coefficient):
            iarg = f.count()
            if self._use_symbols:
                return self._formrep.Dw_sym(iarg, c, der)
            else:
                return self._formrep.Dw_expr(iarg, c, der, False, self._on_facet)

        # --- Handle derivatives of geometry objects
        if isinstance(f, FacetNormal):
            return swiginac.numeric(0.0)

        if isinstance(f, SpatialCoordinate):
            c, = c
            if der[0] == c:
                return swiginac.numeric(1.0)
            else:
                return swiginac.numeric(0.0)

        sfc_error("Eh?")

    def derivative(self, x):
        sfc_error("Derivative shouldn't occur here, you must apply AD first!")

    ### Interior facet stuff:

    def positive_restricted(self, x, y):
        sfc_error("TODO: Restrictions not implemented!")
        return y

    def negative_restricted(self, x, y):
        sfc_error("TODO: Restrictions not implemented!")
        return y


### These require code structure and thus shouldn't occur in SwiginacEvaluator
# (i.e. any conditionals should be handled externally)
#d[EQ] =
#d[NE] =
#d[LE] =
#d[GE] =
#d[LT] =
#d[GT] =

### These are replaced by expand_compounds, so we skip them here:
#d[Identity]   =
#d[Transposed] =
#d[Outer] =
#d[Inner] =
#d[Dot]   =
#d[Cross] =
#d[Trace] =
#d[Determinant] =
#d[Inverse]     =
#d[Deviatoric]  =
#d[Cofactor]    =
#d[Grad] =
#d[Div]  =
#d[Curl] =
#d[Rot]  =

