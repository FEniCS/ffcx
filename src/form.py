__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-09-27"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
from Numeric import *

# FFC modules
import dolfin
import latex
from rank import Rank
from index import *
from algebra import *
from reassign import *
from integrator import Integrator
from finiteelement import FiniteElement

class Form:

    """A Form represents a multi-linear form typically appearing in
    the variational formulation of partial differential equation.

    A Form holds the following data:

        sum   - the representation of the form as a Sum
        ranks - a list of auxiliary data for each Product"""

    def __init__(self, form):
        "Create Form."

        if isinstance(form, Form):
            self.sum = reassign_indices(Sum(form.sum))
            self.ranks = [Rank(p) for p in self.sum.products]
        else:
            self.sum = reassign_indices(Sum(form))
            self.ranks = [Rank(p) for p in self.sum.products]

        # Check that all Products have the same primary rank,
        # otherwise it's not a multi-linear form.
        for i in range(len(self.ranks) - 1):
            if not self.ranks[i].r0 == self.ranks[i + 1].r0:
                raise RuntimeError, "Form must be linear in each of its arguments."

        print "Created form: " + str(self)
        
        return

    def compile(self, language = "C++"):
        "Generate code for evaluation of the variational form."

        # Compute the reference tensor for each product
        A0s = []
        for i in range(len(self.sum.products)):
            A0s += [self.compute_reference_tensor(i)]
            #print "A0 = " + str(A0s[i])

        # Choose language
        if language == "C++":
            dolfin.compile(self.sum.products, A0s, self.ranks)
        elif language == "LaTeX":
            latex.compile(self.sum.products, A0s, self.ranks)
        else:
            print "Unknown language " + str(language)
        return

    def compute_reference_tensor(self, i):
        "Compute the integrals of the reference tensor."

        product = self.sum.products[i]
        rank = self.ranks[i]
        indices0 = [] + rank.indices0
        indices1 = [] + rank.indices1

        # Make sure that the iteration is not empty
        if not indices1:
            indices1 = [[]]

        # Create reference tensor
        A0 = zeros(rank.dims0 + rank.dims1, Float)

        # Create quadrature rule
        integrate = Integrator(product)

        # Iterate over all combinations of indices
        for i in indices0:
            for a in indices1:
                A0[i + a] = integrate(product, i, a)
            
        return A0

    def __repr__(self):
        "Print nicely formatted representation of Form."
        output = "a("
        r0 = self.ranks[0].r0 # All primary ranks are equal
        for i in range(r0):
            if i < (r0 - 1):
                output += "v" + str(i) + ", "
            else:
                output += "v" + str(i) + ") = "
        output += self.sum.__repr__()
        return output

if __name__ == "__main__":

    print "Testing form compiler"
    print "---------------------"

    element = FiniteElement("Lagrange", 1, "triangle")
    
    u = BasisFunction(element)
    v = BasisFunction(element)
    i = Index()
    
    a = Form(u.dx(i)*v.dx(i) + u.dx(0)*v + u*v)
    a.compile()
