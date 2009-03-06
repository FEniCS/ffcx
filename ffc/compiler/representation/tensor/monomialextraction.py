"Extraction of monomial representations of UFL forms."

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2008-08-01 -- 2009-03-06"
__copyright__ = "Copyright (C) 2008-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Martin Alnes, 2008

# UFL modules
from ufl.basisfunction import BasisFunction
from ufl.function import Function
from ufl.constantvalue import ScalarValue, IntValue
from ufl.form import Form
from ufl.algorithms.transformations import purge_list_tensors
from ufl.algorithms.transformations import ReuseTransformer, apply_transformer
from ufl.algorithms.printing import tree_format

# FFC common modules
from ffc.common.log import ffc_assert

# Exception raised when monomial extraction fails
class MonomialException(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class MonomialIndex:
    pass

class MonomialTransform:

    def __init__(self):
        pass

class MonomialFactor:

    def __init__(self, arg=None):
        if isinstance(arg, MonomialFactor):
            self.function = arg.function
            self.component = arg.component
            self.derivative = arg.derivative
        elif isinstance(arg, (BasisFunction, Function)):
            self.function = arg
            self.component = []
            self.derivative = []
        elif arg is None:
            self.function = None
            self.component = []
            self.derivative = []
        else:
            raise MonomialException, ("Unable to create monomial from expression: " + str(arg))

    def element(self):
        return self.function.element()

    def apply_derivative(self, indices):
        self.derivative += indices

    def replace_indices(self, old_indices, new_indices):
        if old_indices is None:
            self.component = new_indices
        else:
            _replace_indices(self.component, old_indices, new_indices)
            _replace_indices(self.derivative, old_indices, new_indices)

    def __str__(self):
        c = ""
        if len(self.component) == 0:
            c = ""
        else:
            c = "[%s]" % ", ".join(str(c) for c in self.component)
        if len(self.derivative) == 0:
            d0 = ""
            d1 = ""
        else:
            d0 = "(" + " ".join("d/dx_%s" % str(d) for d in self.derivative) + " "
            d1 = ")"
        return d0 + str(self.function) + c + d1

class Monomial:
    
    def __init__(self, arg=None):
        if isinstance(arg, Monomial):
            self.float_value = arg.float_value
            self.factors = [MonomialFactor(v) for v in arg.factors]
            self.index_slots = arg.index_slots
        elif isinstance(arg, (MonomialFactor, BasisFunction, Function)):
            self.float_value = 1.0
            self.factors = [MonomialFactor(arg)]
            self.index_slots = None
        elif isinstance(arg, ScalarValue):
            self.float_value = float(arg)
            self.factors = []
            self.index_slots = None
        elif arg is None:
            self.float_value = 1.0
            self.factors = []
            self.index_slots = None
        else:
            raise MonomialException, ("Unable to create monomial from expression: " + str(arg))

    def apply_derivative(self, indices):
        if not len(self.factors) == 1:
            raise MonomialException, "Expecting a single factor."
        self.factors[0].apply_derivative(indices)

    def apply_tensor(self, indices):
        if not self.index_slots is None:
            raise MonomialException, "Expecting scalar-valued expression."
        self.index_slots = indices
    
    def apply_indices(self, indices):
        print "Applying indices:", self.index_slots, "-->", indices
        for v in self.factors:
            v.replace_indices(self.index_slots, indices)
        self.index_slots = None

    def __mul__(self, other):
        m = Monomial()
        m.float_value = self.float_value * other.float_value
        m.factors = self.factors + other.factors
        return m

    def __str__(self):
        if self.float_value == 1.0:
            float_value = ""
        else:
            float_value = "%g * " % self.float_value
        return float_value + " * ".join(str(v) for v in self.factors)

class MonomialSum:

    def __init__(self, arg=None):
        if isinstance(arg, MonomialSum):
            self.monomials = [Monomial(m) for m in arg.monomials]
        elif arg is None:
            self.monomials = []
        else:
            self.monomials = [Monomial(arg)]

    def apply_derivative(self, indices):
        for m in self.monomials:
            m.apply_derivative(indices)

    def apply_tensor(self, indices):
        for m in self.monomials:
            m.apply_tensor(indices)

    def apply_indices(self, indices):
        for m in self.monomials:
            m.apply_indices(indices)

    def __add__(self, other):
        sum = MonomialSum()
        sum.monomials = [Monomial(m) for m in self.monomials] + [Monomial(m) for m in other.monomials]
        return sum

    def __mul__(self, other):
        sum = MonomialSum()
        for m0 in self.monomials:
            for m1 in other.monomials:
                sum.monomials.append(m0 * m1)
        return sum

    def __str__(self):
        return " + ".join(str(m) for m in self.monomials)

class MonomialForm:

    def __init__(self):
        self.integrals = []

    def append(self, integral, measure):
        self.integrals.append((integral, measure))

    def __len__(self):
        return len(self.integrals)

    def __getitem__(self, i):
        return self.integrals[i]

    def __iter__(self):
        return iter(self.integrals)

    def __str__(self):
        s  = "Monomial form of %d integral(s)\n" % len(self.integrals)
        s += len(s) * "-" + "\n"
        for (integrand, measure) in self.integrals:
            s += "Integrand: " + str(integrand) + "\n"
            s += "Measure:   " + str(measure) + "\n"
        return s

class MonomialTransformer(ReuseTransformer):

    def __init__(self):
        ReuseTransformer.__init__(self)
    
    def expr(self, o, *ops):
        raise MonomialException, ("No handler defined for expression %s." % o._uflclass.__name__)

    def terminal(self, o):
        raise MonomialException, ("No handler defined for terminal %s." % o._uflclass.__name__)

    def variable(self, o):
        return self.visit(o.expression())

    #--- Operator handles ---

    def sum(self, o, s0, s1):
        print "\nSum"
        s = s0 + s1
        print "Result:", s
        return s

    def product(self, o, s0, s1):
        print "\nProduct: [%s] * [%s]" % (str(s0), str(s1))
        s = s0 * s1
        print "Result:", s
        return s

    def index_sum(self, o, s, index):
        print "\nIgnoring IndexSum expression for now"
        print "Result:", s
        return s

    def indexed(self, o, s, indices): 
        print "\nIndexed", s, indices
        s = MonomialSum(s)
        s.apply_indices(indices)
        print "Result:", s
        return s

    def component_tensor(self, o, s, indices):
        print "\nComponentTensor", s, indices
        s = MonomialSum(s)
        s.apply_tensor(indices)
        print "Result:", s
        return s

    def spatial_derivative(self, o, s, indices):
        print "\nSpatialDerivative", s, indices
        s = MonomialSum(s)
        s.apply_derivative(indices)
        print "Result:", s
        return s

    def power(self, o, s, ignored_exponent_expressed_as_sum):
        print "\nPower", s, ignored_exponent_expressed_as_sum
        (expr, exponent) = o.operands()
        if not isinstance(exponent, IntValue):
            raise MonomialException, "Cannot handle non-integer exponents."
        p = MonomialSum(Monomial())
        for i in range(int(exponent)):
            p = p * s
        return p

    #--- Terminal handlers ---

    def multi_index(self, multi_index):
        print "\nMultiIndex"
        indices = [index for index in multi_index]
        print indices
        return indices

    def index(self, o):
        raise MonomialException, "Not expecting to see an Index terminal."

    def basis_function(self, v):
        print "\nBasisFunction", v
        s = MonomialSum(v)
        print "Result:", s
        return s

    def function(self, v):
        print "\nFunction", v
        s = MonomialSum(v)
        print "Result:", s
        return s

    def scalar_value(self, x):
        print "\nScalarValue", x
        s = MonomialSum(x)
        print "Result:", s
        return s

def extract_monomial_form(form):
    """Extract monomial representation of form (if possible). When
    successful, the form is represented as a sum of products of scalar
    components of basis functions or derivatives of basis functions.
    The sum of products is represented as a tuple of tuples of basis
    functions. If unsuccessful, MonomialException is raised."""

    # Check that we get a Form
    ffc_assert(isinstance(form, Form), "Expecting a UFL form.")

    print ""
    print "Extracting monomials"
    print "--------------------"
    print ""
    
    # Extract processed form
    form_data = form.form_data()
    form = form_data.form

    # Purge list tensors from expression tree
    form = purge_list_tensors(form)

    # Iterate over all integrals
    monomial_form = MonomialForm()
    for integral in form.cell_integrals():

        # Get measure and integrand
        measure = integral.measure()
        integrand = integral.integrand()
        print tree_format(integrand)

        # Extract monomial representation if possible
        integrand = apply_transformer(integrand, MonomialTransformer())
        monomial_form.append(integrand, measure)

    return monomial_form


def _replace_indices(indices, old_indices, new_indices):
    "Handle replacement of subsets of multi indices."

    # Old and new indices must match
    if not len(old_indices) == len(new_indices):
        raise MonomialException, "Unable to replace indices, mismatching index dimensions."

    # Build index map
    index_map = {}
    for (i, index) in enumerate(old_indices):
        index_map[index] = new_indices[i]

    # Check all indices and replace
    for (i, index) in enumerate(indices):
        if index in old_indices:
            indices[i] = index_map[index]
