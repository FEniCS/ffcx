"Extraction of monomial representations of UFL forms."

# Copyright (C) 2008-2013 Anders Logg
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Martin Alnaes, 2008, 2013
# Modified by Kristian B. Oelgaard
#
# First added:  2008-08-01
# Last changed: 2013-01-08

# UFL modules
from ufl.classes import Form, Argument, Coefficient, ScalarValue, IntValue
from ufl.algorithms import purge_list_tensors, apply_transformer, ReuseTransformer

# FFC modules
from ffc.log import info, debug, ffc_assert

# Cache for computed integrand representations
#_cache = {}

def extract_monomial_form(integrals, function_replace_map):
    """
    Extract monomial representation of form (if possible). When
    successful, the form is represented as a sum of products of scalar
    components of basis functions or derivatives of basis functions.
    If unsuccessful, MonomialException is raised.
    """

    info("Extracting monomial form representation from UFL form")

    # Iterate over all integrals
    monomial_form = MonomialForm()
    for integral in integrals:

        # Get measure and integrand
        measure = integral.measure()
        if measure.domain_type() == "point":
            raise MonomialException("Point integrals are not supported by tensor representation.")

        integrand = integral.integrand()

        # Extract monomial representation if possible
        integrand = extract_monomial_integrand(integrand, function_replace_map)
        monomial_form.append(integrand, measure)

    return monomial_form

def extract_monomial_integrand(integrand, function_replace_map):
    "Extract monomial integrand (if possible)."

    # Check cache
    #if integrand in _cache:
    #    debug("Reusing monomial integrand from cache")
    #    return _cache[integrand]

    # Purge list tensors
    integrand = purge_list_tensors(integrand)

    # Apply monomial transformer
    monomial_integrand = apply_transformer(integrand, MonomialTransformer(function_replace_map))

    # Store in cache
    #_cache[integrand] = monomial_integrand

    return monomial_integrand

class MonomialException(Exception):
    "Exception raised when monomial extraction fails."
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class MonomialFactor:
    """
    This class represents a monomial factor, that is, a derivative of
    a scalar component of a basis function.
    """

    def __init__(self, arg=None):
        if isinstance(arg, MonomialFactor):
            self.function = arg.function
            self.components = arg.components
            self.derivatives = arg.derivatives
            self.restriction = arg.restriction
        elif isinstance(arg, (Argument, Coefficient)):
            self.function = arg
            self.components = []
            self.derivatives = []
            self.restriction = None
        elif arg is None:
            self.function = None
            self.components = []
            self.derivatives = []
            self.restriction = None
        else:
            raise MonomialException, ("Unable to create monomial from expression: " + str(arg))

    def element(self):
        return self.function.element()

    def count(self):
        return self.function.count()

    def apply_derivative(self, indices):
        self.derivatives += indices

    def apply_restriction(self, restriction):
        self.restriction = restriction

    def replace_indices(self, old_indices, new_indices):
        if old_indices is None:
            self.components = new_indices
        else:
            _replace_indices(self.components, old_indices, new_indices)
            _replace_indices(self.derivatives, old_indices, new_indices)

    def __str__(self):
        if len(self.components) == 0:
            c = ""
        else:
            c = "[%s]" % ", ".join(str(c) for c in self.components)
        if len(self.derivatives) == 0:
            d0 = ""
            d1 = ""
        else:
            d0 = "(" + " ".join("d/dx_%s" % str(d) for d in self.derivatives) + " "
            d1 = ")"
        if self.restriction is None:
            r = ""
        else:
            r = "(%s)" % str(self.restriction)
        return d0 + str(self.function) + r + c + d1

class Monomial:
    "This class represents a product of monomial factors."

    def __init__(self, arg=None):
        if isinstance(arg, Monomial):
            self.float_value = arg.float_value
            self.factors = [MonomialFactor(v) for v in arg.factors]
            self.index_slots = arg.index_slots
        elif isinstance(arg, (MonomialFactor, Argument, Coefficient)):
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
        for v in self.factors:
            v.replace_indices(self.index_slots, indices)
        self.index_slots = None

    def apply_restriction(self, restriction):
        for v in self.factors:
            v.apply_restriction(restriction)

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
    "This class represents a sum of monomials."

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

    def apply_restriction(self, restriction):
        for m in self.monomials:
            m.apply_restriction(restriction)

    def __add__(self, other):
        m0 = [Monomial(m) for m in self.monomials]
        m1 = [Monomial(m) for m in other.monomials]
        sum = MonomialSum()
        sum.monomials = m0 + m1
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
    """
    This class represents a monomial form, that is, a sum of
    integrals, each represented as a MonomialSum.
    """

    def __init__(self):
        self.integrals = []

    def append(self, integrand, measure):
        self.integrals.append((integrand, measure))

    def __len__(self):
        return len(self.integrals)

    def __getitem__(self, i):
        return self.integrals[i]

    def __iter__(self):
        return iter(self.integrals)

    def __str__(self):
        if len(self.integrals) == 0:
            return "<Empty form>"
        s  = "Monomial form of %d integral(s)\n" % len(self.integrals)
        s += len(s) * "-" + "\n"
        for (integrand, measure) in self.integrals:
            s += "Integrand: " + str(integrand) + "\n"
            s += "Measure:   " + str(measure) + "\n"
        return s

class MonomialTransformer(ReuseTransformer):
    """
    This class defines the transformation rules for extraction of a
    monomial form represented as a MonomialSum from a UFL integral.
    """

    def __init__(self, function_replace_map=None):
        ReuseTransformer.__init__(self)
        self._function_replace_map = function_replace_map or {}

    def expr(self, o, *ops):
        raise MonomialException("No handler defined for expression %s." % o._uflclass.__name__)

    def terminal(self, o):
        raise MonomialException("No handler defined for terminal %s." % o._uflclass.__name__)

    def variable(self, o):
        return self.visit(o.expression())

    #--- Operator handles ---

    def division(self, o):

        # Handle division by scalars as multiplication by inverse
        denominator = o.operands()[1]
        if not isinstance(denominator, ScalarValue):
            raise MonomialException("No handler defined for expression %s."
                                    % o._uflclass.__name__)
        inverse = self.scalar_value(ScalarValue(1.0/denominator.value()))
        numerator = self.visit(o.operands()[0])

        return self.product(o, inverse, numerator)

    def sum(self, o, s0, s1):
        s = s0 + s1
        return s

    def product(self, o, s0, s1):
        s = s0 * s1
        return s

    def index_sum(self, o, s, index):
        return s

    def indexed(self, o, s, indices):
        s = MonomialSum(s)
        s.apply_indices(indices)
        return s

    def component_tensor(self, o, s, indices):
        s = MonomialSum(s)
        s.apply_tensor(indices)
        return s

    def grad(self, o, s):
        # The representation
        #   o = Grad(s)
        # is equivalent to
        #   o = as_tensor(s[ii].dx(i), ii+(i,))
        # with ii a tuple of free indices and i a single free index.

        # Make some unique utility indices
        from ufl import indices
        on = o.rank()
        ind = list(indices(on))

        # Get underlying expression to take gradient of
        f, = o.operands()
        fn = f.rank()

        ffc_assert(on == fn + 1, "Assuming grad always adds one axis.")

        s = MonomialSum(s)
        s.apply_indices(list(ind[:-1]))

        s = MonomialSum(s)
        s.apply_derivative([ind[-1]])

        s = MonomialSum(s)
        s.apply_tensor(ind)
        return s

    def positive_restricted(self, o, s):
        s.apply_restriction("+")
        return s

    def negative_restricted(self, o, s):
        s.apply_restriction("-")
        return s

    def power(self, o, s, ignored_exponent_expressed_as_sum):
        (expr, exponent) = o.operands()
        if not isinstance(exponent, IntValue):
            raise MonomialException("Cannot handle non-integer exponents.")
        e = int(exponent)
        if e < 0:
            raise MonomialException("Cannot handle negative exponents.")
        p = MonomialSum(Monomial())
        for i in range(e):
            p = p * s
        return p

    #--- Terminal handlers ---

    def multi_index(self, multi_index):
        indices = [index for index in multi_index]
        return indices

    def argument(self, v):
        s = MonomialSum(self._function_replace_map.get(v,v))
        return s

    def coefficient(self, v):
        s = MonomialSum(self._function_replace_map.get(v,v))
        return s

    def scalar_value(self, x):
        s = MonomialSum(x)
        return s

def _replace_indices(indices, old_indices, new_indices):
    "Handle replacement of subsets of multi indices."

    # Old and new indices must match
    if not len(old_indices) == len(new_indices):
        raise MonomialException("Unable to replace indices, mismatching index dimensions.")

    # Build index map
    index_map = {}
    for (i, index) in enumerate(old_indices):
        index_map[index] = new_indices[i]

    # Check all indices and replace
    for (i, index) in enumerate(indices):
        if index in old_indices:
            indices[i] = index_map[index]
