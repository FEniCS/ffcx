"Transformation of monomial representations of UFL forms."

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2009-03-06 -- 2009-03-06"
__copyright__ = "Copyright (C) 2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# UFL modules
from ufl import BasisFunction, Function
from ufl.indexing import FixedIndex

# FFC common modules
from ffc.common.log import ffc_assert

# FFC fem modules
from ffc.fem import create_element

# FFC tensor representation modules
from monomialextraction import MonomialForm

# Index counters
_current_secondary_index = 0
_current_auxiliary_index = 0

def next_secondary_index():
    global _current_secondary_index
    _current_secondary_index += 1
    return _current_secondary_index - 1

def next_auxiliary_index():
    global _current_auxiliary_index
    _current_auxiliary_index += 1
    return _current_auxiliary_index - 1

def reset_indices():
    global _current_secondary_index
    global _current_auxiliary_index
    _current_secondary_index = 0
    _current_auxiliary_index = 0

class MonomialIndex:

    PRIMARY = "primary"
    SECONDARY = "secondary"
    AUXILIARY = "auxiliary"

    def __init__(self, type, index=None):
        if type == MonomialIndex.SECONDARY and index is None:
            index = next_secondary_index()
        elif type == MonomialIndex.AUXILIARY and index is None:
            index = next_auxiliary_index()
        self.type = type
        self.index = index

    def __str__(self):
        if self.type == MonomialIndex.PRIMARY:
            return "i_" + str(self.index)
        elif self.type == MonomialIndex.SECONDARY:
            return "a_" + str(self.index)
        elif self.type == MonomialIndex.AUXILIARY:
            return "b_" + str(self.index)
        elif self.type == MonomialIndex.FIXED:
            return str(self.index)

class MonomialDeterminant:

    def __init__(self):
        self.power = 1

    def __str__(self):
        if self.power == 1:
            return "(det F')"
        else:
            return "(det F')^%s" + str(self.power)

class MonomialCoefficient:

    def __init__(self, index):
        self.index = index

    def __str__(self):
        return "c_" + str(self.index)

class MonomialTransform:

    def __init__(self, index0, index1):
        self.index0 = index0
        self.index1 = index1

    def __str__(self):
        return "dX_%s/dx_%s" % (str(self.index0), str(self.index1))

class MonomialBasisFunction:

    def __init__(self, element, index, components, derivatives):
        self.element = element
        self.index = index
        self.components = components
        self.derivatives = derivatives

    def __str__(self):
        c = ""
        if len(self.components) == 0:
            c = ""
        else:
            c = "[%s]" % ", ".join(str(c) for c in self.components)
        if len(self.derivatives) == 0:
            d0 = ""
            d1 = ""
        else:
            d0 = "(" + " ".join("d/dX_%s" % str(d) for d in self.derivatives) + " "
            d1 = ")"
        v = "V_" + str(self.index)
        return d0 + v + c + d1

class TransformedMonomial:

    def __init__(self, monomial):

        # Reset monomial data
        self.float_value = 1.0
        self.determinant = MonomialDeterminant()
        self.coefficients = []
        self.transforms = []
        self.basis_functions = []

        # Reset index counters
        reset_indices()

        # Initialize index map
        index_map = {}

        # Iterate over factors
        for f in monomial.factors:

            # Extract basis function index and coefficients
            if isinstance(f.function, BasisFunction):
                vindex = f.function.count()
            elif isinstance(f.function, Function):
                index = MonomialIndex(MonomialIndex.SECONDARY)
                coefficient = MonomialCoefficient(index)
                vindex = index
                self.coefficients.append(coefficient)

            # Extract components
            components = []
            for c in f.components:
                if c in index_map:
                    index = index_map[c]
                elif isinstance(c, FixedIndex):
                    index = MonomialIndex(type=MonomialIndex.FIXED, value=c.count())
                else:
                    index = MonomialIndex(MonomialIndex.AUXILIARY)
                index_map[c] = index
                components.append(index)

            # Extract derivatives / transforms
            derivatives = []
            for d in f.derivatives:
                index0 = MonomialIndex(MonomialIndex.SECONDARY)
                if d in index_map:
                    index1 = index_map[d]
                elif isinstance(d, FixedIndex):
                    index1 = MonomialIndex(type=MonomialIndex.FIXED, value=d.count())
                else:
                    index1 = MonomialIndex(MonomialIndex.AUXILIARY)
                index_map[d] = index1
                transform = MonomialTransform(index0, index1)
                self.transforms.append(transform)                
                derivatives.append(index0)

            # Extract element
            element = create_element(f.element())

            # Create basis function
            v = MonomialBasisFunction(element, vindex, components, derivatives)
            self.basis_functions.append(v)

    def __str__(self):
        factors = []
        if not self.float_value == 1.0:
            factors.append(float_value)
        factors.append(self.determinant)
        factors += self.coefficients
        factors += self.transforms
        return " * ".join([str(f) for f in factors]) + " | " + " * ".join([str(v) for v in self.basis_functions])
    
def transform_monomial_form(monomial_form):
    "Transform monomial form to reference element."

    # Check that we get a Form
    ffc_assert(isinstance(monomial_form, MonomialForm), "Expecting a MonomialForm.")

    # Transform each monomial
    for (integrand, measure) in monomial_form:
        for (i, monomial) in enumerate(integrand.monomials):
            integrand.monomials[i] = TransformedMonomial(monomial)
