__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2008-09-04 -- 2008-09-11"
__copyright__ = "Copyright (C) 2008 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

from hashlib import sha1
from instant import get_swig_version

# FFC compiler modules
from ffc.compiler.analysis import simplify, analyze
from ffc.common.constants  import FFC_VERSION

try:
    from ufl.algorithms import FormData
except:
    pass

class JITObject:
    """This class is a wrapper for a compiled object in the context of
    specific compiler options. A JITObject is identified either by its
    hash value or by its signature. The hash value is valid only in a
    single instance of an application (at runtime). The signature is
    persistent and may be used for caching modules on disk."""

    def __init__(self, form, options):
        "Create JITObject for given form and options"
        self.form = form
        self.options = options
        self._form_data = None
        self._hash = None
        self._signature = None

    def __hash__(self):
        "Return unique integer for form + options"

        # Check if we have computed the hash before
        if not self._hash is None:
            return self._hash

        # Compute hash
        string = str(id(self.form)) + str(self.options)
        hexdigest = sha1(string).hexdigest()
        number = int(hexdigest, 16)
        
        return number

    def __eq__(self, other):
        "Check for equality"
        return hash(self) == hash(other)
    
    def signature(self):
        "Return unique string for form expression + options"
        
        # Check if we have computed the signature before
        if not self._signature is None:
            return self._signature
        
        # Compute signature
        if self.options["compiler"] == "ffc":
            self.form_data = analyze.analyze(self.form, simplify_form=True)
            element_signature = ";".join([element.signature() for element in self.form_data.elements])
        elif self.options["compiler"] == "ufl":
            self.form_data = FormData(self.form)
            element_signature = ";".join([element.__repr__() for element in self.form_data.elements])
        else:
            raise RuntimeError("Unknown compiler.")

        form_signature = str(self.form)
        swig_version = get_swig_version()
        options_signature = str(self.options)
        string = ";".join([form_signature, element_signature, swig_version, \
                           options_signature, FFC_VERSION])
        self._signature = "form_" + sha1(string).hexdigest()

        # Store form data and signature in form for later reuse
#        self._form_data = self.form_data

        return self._signature
