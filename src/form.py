__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-09-27"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
from Numeric import *

# FFC modules
import dolfin
import latex
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

    def compile(self, language = "C++"):
        "Generate code for evaluation of the variational form."

        # Choose format
        if language == "C++":
            format = dolfin.format
        elif language == "LaTeX":
            format = latex.format
        else:
            print "Unknown language " + str(language)

        # Reassign indices
        reassign_indices(self.sum)

        # Create element tensors
        self.AKi = ElementTensor(self.sum, "interior", format)
        self.AKb = ElementTensor(self.sum, "boundary", format)

        # Check ranks
        self.__check_primary_ranks()
            
        # Generate output
        if language == "C++":
            dolfin.compile(self)
        elif language == "LaTeX":
            latex.compile(self)
        else:
            print "Unknown language " + str(language)
            
        print "Compiled form: " + str(self)
            
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

    element = FiniteElement("Lagrange", 1, "tetrahedron")
    
    u = BasisFunction(element)
    v = BasisFunction(element)
    i = Index()
    dx = Integral("interior")
    ds = Integral("boundary")
    
    a = Form(u.dx(i)*v.dx(i)*dx, "FFCPoisson")
    a.compile("C++")
    #a.compile("LaTeX")
    
    #element = FiniteElement("Vector Lagrange", 1, "triangle")
    #u = BasisFunction(element)
    #v = BasisFunction(element)
    #i = Index()
    #j = Index()
    #a = Form(u[i].dx(j)*v[i].dx(j), "FFCPoissonSystem")
    #a.compile()
