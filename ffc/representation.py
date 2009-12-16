__author__ = "Anders Logg (logg@simula.no) and friends"
__date__ = "2009-12-16"
__copyright__ = "Copyright (C) 2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2009-12-16

def form_representation(form, form_data, method):
    "Compute and return intermediate representation of form."

    # FIXME: Call correct method in quadrature or tensor module

    if method == "quadrature":
        return {}
    else:
        return {}

def element_representation(ufl_element):
    "Compute and return intermediate representation of element."

    return {}

def dofmap_representation(ufl_element):
    "Compute and return intermediate representation of dofmap."

    return {}
