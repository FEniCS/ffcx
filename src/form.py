__author__ = "Anders Logg (logg@tti-c.org)"
__version__ = "0.0.1"
__date__ = "2004-09-27"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

from Numeric import *
from algebra import *

class Form:
    """A Form represents a multi-linear form typically appearing in
    the variational formulation of partial differential equation."""

    def __init__(self, form):

        # Make sure that we get a Sum
        if isinstance(form, Form):
            self.sum = form.sum
        elif isinstance(form, BasisFunction):
            self.sum = Sum(form)
        elif isinstance(form, Factor):
            self.sum = Sum(form)
        elif isinstance(form, Product):
            self.sum = Sum(form)
        elif isinstance(form, Sum):
            self.sum = form
        else:
            raise RuntimeError, "Cannot create Form from " + str(form)

        # Compute the rank of the tensor for each Product. The
        # secondary rank may differ for each Product.
        self.r0 = []
        self.r1 = []
        for p in self.sum.products:
            [tmp0, tmp1] = p.rank()
            self.r0 += [tmp0]
            self.r1 += [tmp1]
        print "Primary rank is " + str(self.r0) + ", secondary rank is " + str(self.r1)

        # Check that all Products have the same primary rank,
        # otherwise it's not a multi-linear form.
        for i in range(len(self.r0) - 1):
            if not self.r0[i] == self.r0[i + 1]:
                raise RuntimeError, "Form must be linear in each of its arguments."

        # Compute dimensions of the tensor for each Product. Note that
        # we have one list of dimensions for each Product in the Sum.
        self.dims = []
        for i in range(len(self.sum.products)):
            self.dims += [self.sum.products[i].dims(self.r0[i], self.r1[i])]
        
        return

    def compile(self, language = "C++"):
        "Generate code for evaluation of the variational form."
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

        # Compute the reference tensors for each product
        for i in range(len(self.sum.products)):
            A0 = self.compute_reference_tensor(self.sum.products[i],
                                               self.r0[i], self.r1[i],
                                               self.dims[i])
            print A0

        # Pass the reference tensor to FErari for simplification
        # Make call to FErari here
        
        return

    def compile_latex(self):
        "Generate LaTeX code (paste into your report :-)."
        print "Compiling multi-linear form for LaTeX."
        print "Not yet supported."
        return

    def compute_reference_tensor(self, product, r0, r1, dims):
        "Compute the integrals of the reference tensor using FIAT."
        A0 = zeros(dims, Float)
        indices = self.build_indices(dims)
        for index in indices:
            A0[index] = self.integrate(product, index, r0, r1)
        return A0

    def integrate(self, product, index, r0, r1):
        return sum(index)

    def build_indices(self, dims):
        """Create a list of all indices of the reference tensor.
        Someone please tell me if there is a better way to iterate
        over a multi-dimensional array. The list of indices is
        constucted by iterating over all integer numbers that can be
        represented with the base of each position determined by the
        given dimension for that dimension."""
        current = zeros(len(dims))
        indices = []
        posvalue = [1] + list(cumproduct(dims)[:-1])
        for i in range(product(dims)):
            sum = i
            for pos in range(len(dims)):
                current[pos] = (i / posvalue[pos]) % dims[pos]
                sum -= current[pos] * posvalue[pos]
            indices += [list(current)]
        return array(indices)

    def __repr__(self):
        output = "a("
        for i in range(self.r0[0]):
            if i < (self.r0[0] - 1):
                output += "v" + str(i) + ", "
            else:
                output += "v" + str(i) + ") = "
        output += self.sum.__repr__()
        return output

if __name__ == "__main__":

    print "Testing form compiler"
    print "---------------------"

    u = BasisFunction()
    v = BasisFunction()
    i = Index()

    a = Form(u.dx(i)*v.dx(i) + u*v)
    print a
    a.compile()
