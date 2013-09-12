# Copyright (C) 2008-2013 Anders Logg
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Martin Alnaes, 2013
#
# First added:  2008-09-04
# Last changed: 2013-01-25

# Python modules.
from hashlib import sha1

# Instant modules.
from instant import get_swig_version

# UFL modules.
import ufl

# FFC modules.
from constants import FFC_VERSION

# UFC modules.
import ufc_utils

# Compute signature of all ufc headers combined
ufc_signature = sha1(''.join(getattr(ufc_utils, header)
                             for header in
                             (k for k in vars(ufc_utils).keys()
                              if k.endswith("_header")))
                              ).hexdigest()

class JITObject:
    """This class is a wrapper for a compiled object in the context of
    specific compiler parameters. A JITObject is identified either by its
    hash value or by its signature. The hash value is valid only in a
    single instance of an application (at runtime). The signature is
    persistent and may be used for caching modules on disk."""

    def __init__(self, form, parameters):
        "Create JITObject for given form and parameters"
        assert(isinstance(form, ufl.Form))

        # Store data
        self.form = form
        self.parameters = parameters
        self._hash = None
        self._signature = None

    def __hash__(self):
        "Return unique integer for form + parameters"

        # Check if we have computed the hash before
        if not self._hash is None:
            return self._hash

        # Compute hash based on signature
        self._hash = int(self.signature(), 16)

        return self._hash

    def __eq__(self, other):
        "Check for equality"
        return hash(self) == hash(other)

    def signature(self):
        "Return unique string for form + parameters"

        # Check if we have computed the signature before
        if not self._signature is None:
            return self._signature

        # Get signature from assumed precomputed form_data
        form_signature = self.form.form_data().signature

        # Compute other relevant signatures
        parameters_signature = _parameters_signature(self.parameters)
        ffc_signature = str(FFC_VERSION)
        swig_signature = str(get_swig_version())
        cell_signature = str(self.form.form_data().cell)

        # Build common signature
        signatures = [form_signature,
                      parameters_signature,
                      ffc_signature,
                      cell_signature,
                      ufc_signature]
        string = ";".join(signatures)
        self._signature = sha1(string).hexdigest()

        # Uncomment for debugging
        #print "form_signature       =", form_signature
        #print "parameters_signature =", parameters_signature
        #print "ffc_signature        =", ffc_signature
        #print "cell_signature       =", cell_signature
        #print "signature            =", self._signature

        return self._signature

def _parameters_signature(parameters):
    "Return parameters signature (some parameters must be ignored)."
    parameters = parameters.copy()
    ignores = ["log_prefix"]
    for ignore in ignores:
        if ignore in parameters:
            del parameters[ignore]
    return str(parameters)
