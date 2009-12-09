__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2008-09-04"
__copyright__ = "Copyright (C) 2008 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2009-12-09

from hashlib import sha1
from instant import get_swig_version

# FFC common modules
from ffc.common.constants import FFC_VERSION

# UFL modules
import ufl

class JITObject:
    """This class is a wrapper for a compiled object in the context of
    specific compiler options. A JITObject is identified either by its
    hash value or by its signature. The hash value is valid only in a
    single instance of an application (at runtime). The signature is
    persistent and may be used for caching modules on disk."""

    def __init__(self, form, options):
        "Create JITObject for given form and options"
        assert(isinstance(form, ufl.Form))

        # Store data
        self.form = form
        self.options = options
        self._hash = None
        self._signature = None

    def __hash__(self):
        "Return unique integer for form + options"

        # Compute hash if not computed before
        if self._hash is None:
            string = str(id(self.form)) + _options_signature(self.options)
            hexdigest = sha1(string).hexdigest()
            self._hash = int(hexdigest, 16)

        return self._hash

    def __eq__(self, other):
        "Check for equality"
        return hash(self) == hash(other)

    def signature(self):
        "Return unique string for form expression + options"

        # Check if we have computed the signature before
        if not self._signature is None:
            return self._signature

        # Compute form signature based on form stored in formdata
        form_signature = repr(self.form)

        # Build signature including form, options, FFC version and SWIG version
        options_signature = _options_signature(self.options)
        ffc_signature     = str(FFC_VERSION)
        swig_signature    = str(get_swig_version())
        signatures = [form_signature, options_signature, ffc_signature, swig_signature]
        string = ";".join(signatures)
        self._signature = "form_" + sha1(string).hexdigest()

        return self._signature

def _options_signature(options):
    "Return options signature (some options must be ignored)."
    options = options.copy()
    ignores = ["log_prefix"]
    for ignore in ignores:
        if ignore in options:
            del options[ignore]
    return str(options)
