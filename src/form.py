__author__ = "Anders Logg (logg@tti-c.org)"
__version__ = "0.0.1"
__date__ = "2004-09-27"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

from algebra import *

class Form:
    """A Form represents a multi-linear form typically appearing in
    the variational formulation of partial differential equation."""

    def __init__(self, form):
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
        print "Compiling multi-linear for for LaTeX."
        print "Not yet supported."
        return

    def compile_cpp(self):
        "Generate C++ code."
        print "Compiling multi-linear for for C++."

        # Compute the reference tensor
        A0 = self.compute_reference_tensor()

        # Pass the reference tensor to FErari for simplification
        # Make call to FErari here
        
        return

    def compile_latex(self):
        "Generate LaTeX code (paste into your report :-)."
        print "Compiling multi-linear for for LaTeX."
        print "Not yet supported."
        return

    def compute_reference_tensor(self):
        "Compute the integrals of the reference tensors using FIAT."

        # Compute the rank of the tensor
#        for p in self.sum.products:
#            [r0, r1, rtot] = p.rank()

if __name__ == "__main__":

    print "Testing form compiler"
    print "---------------------"

    u = BasisFunction()
    v = BasisFunction()
    i = Index()

    a = Form(u.dx(i)*v.dx(i))
    a.compile()
