__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-09-27"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
import sys
from Numeric import *

# FFC format modules
sys.path.append("../../")
from ffc.format import dolfin
from ffc.format import latex

# FFC compiler modules
from index import *
from algebra import *
from integral import *
from reassign import *
from finiteelement import *
from elementtensor import *

class Form:

    """A Form represents a multi-linear form typically appearing in
    the variational formulation of partial differential equation.
    
    A Form holds the following data:

        sum     - a Sum representing the multi-linear form
        name    - a string, the name of the multi-linear form
        AKi     - interior ElementTensor
        AKb     - boundary ElementTensor
        rank    - primary rank of the multi-linear form
        dims    - list of primary dimensions
        indices - list of primary indices

    A multi-linear form is first expressed as an element of the
    algebra (a Sum) and is then post-processed to be expressed as a
    sum of Terms, where each Term is the product of a ReferenceTensor
    and a GeometryTensor."""

    def __init__(self, sum, name = "MyPDE"):
        "Create Form."

        # Initialize Form
        self.sum     = Sum(sum)
        self.name    = name
        self.AKi     = None
        self.AKb     = None
        self.rank    = None
        self.dims    = None
        self.indices = None

        return

    def compile(self, language = None):
        "Generate code for evaluation of the variational form."

        # Choose language
        if not language:
            compiler = dolfin
        elif language == "C++" or language == "c++":
            compiler = dolfin
        elif language == "LaTeX" or language == "latex":
            compiler = latex
        else:
            raise "RuntimeError", "Unknown language " + str(language)

        # Reassign indices
        reassign_indices(self.sum)
        print "Compiling form: " + str(self.sum)
        
        # Create element tensors
        self.AKi = ElementTensor(self.sum, "interior", compiler.format)
        self.AKb = ElementTensor(self.sum, "boundary", compiler.format)

        # Check primary ranks
        self.__check_primary_ranks()

        # Generate output
        compiler.compile(self)

        return

    def __check_primary_ranks(self):
        "Check that all primary ranks are equal."
        terms = self.AKi.terms + self.AKb.terms
        ranks = [term.A0.i.rank for term in terms]
        if not ranks[1:] == ranks[:-1]:
            "Form must be linear in each of its arguments."
        self.rank = ranks[0]
        self.dims = terms[0].A0.i.dims
        self.indices = terms[0].A0.i.indices
        return

    def __repr__(self):
        "Print nicely formatted representation of Form."
        v = ", ".join(["vi" + str(i) for i in range(self.rank)])
        return "a(%s) = %s" % (v, self.sum.__repr__())

if __name__ == "__main__":

    print "Testing form compiler"
    print "---------------------"

    element = FiniteElement("Lagrange", "triangle", 1)
    
    u = BasisFunction(element)
    v = BasisFunction(element)
    i = Index()
    dx = Integral("interior")
    ds = Integral("boundary")
    
    a = Form(u.dx(i)*v.dx(i)*dx, "FFCPoisson")
    a.compile("C++")
    a.compile("LaTeX")
