__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-09-27"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
from Numeric import *

# FFC modules
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

        return

        # Check that all Products have the same primary rank,
        # otherwise it's not a multi-linear form.
        for i in range(len(self.ranks) - 1):
            if not self.ranks[i].r0 == self.ranks[i + 1].r0:
                raise RuntimeError, "Form must be linear in each of its arguments."
            
        return

    def compile(self, language = "C++"):
        "Generate code for evaluation of the variational form."

        # Compute the reference tensor for each product
        for i in range(len(self.sum.products)):
            A0 = self.compute_reference_tensor(i)
            print "A0 = " + str(A0)

        # Choos language
        if language == "C":
            self.compile_c()
        elif language == "C++":
            self.compile_cpp()
        elif language == "LaTeX":
            self.compile_latex()
        else:
            print "Unknown language " + str(language)
        return

    def compile_c(self):
        "Generate C code."
        print "Compiling multi-linear form for C."
        print "Not yet supported."
        return

    def compile_cpp(self):
        "Generate C++ code."
        print "Compiling multi-linear form for C++."
        print "Not yet supported."
        return

    def compile_latex(self):
        "Generate LaTeX code (paste into your report :-)."
        print "Compiling multi-linear form for LaTeX."
        print "Not yet supported."
        return

    def compute_reference_tensor(self, i):
        "Compute the integrals of the reference tensor."

        product = self.sum.products[i]
        r0 = self.ranks[i].r0
        r1 = self.ranks[i].r1
        dims = self.ranks[i].dims

        # Create reference tensor and a list of all indices
        A0 = zeros(dims, Float)
        indices = build_indices(dims)

        # Create quadrature rule
        integrate = Integrator(product)

        # Iterate over all combinations of indices
        for index in indices:
            A0[index] = integrate(product, index, r0, r1)
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
    
    a = Form(u.dx(i)*v.dx(i) + u*v)

    print a
    a.compile()
