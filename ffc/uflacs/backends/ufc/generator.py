# -*- coding: utf-8 -*-
# Copyright (C) 2015-2017 Martin Sandve Alnæs
#
# This file is part of UFLACS.
#
# UFLACS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFLACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>.

import re

from ffc.log import error, warning
import ffc.backends.ufc

from ffc.uflacs.language.format_lines import format_indented_lines
from ffc.uflacs.backends.ufc.templates import *

#__all__ = (["ufc_form", "ufc_dofmap", "ufc_finite_element", "ufc_integral"]
#           + ["ufc_%s_integral" % integral_type for integral_type in integral_types])


# These are all the integral types directly supported in ufc
from ffc.representation import ufc_integral_types


# These are the method names in ufc::form that are specialized for each integral type
integral_name_templates = (
    "max_%s_subdomain_id",
    "has_%s_integrals",
    "create_%s_integral",
    "create_default_%s_integral",
    )


class ufc_generator(object):
    """Common functionality for code generators producing ufc classes.

    The generate function is the driver for generating code for a class.
    It automatically extracts each template keyword %(foo)s and calls
    self.foo(...) to define the code snippet for that keyword.

    The standard arguments to every self.foo(...) function are:
    - L (the language module reference)
    - ir (the full ir dict)
    - parameters (the parameters dict, can be omitted)
    If the second argument is not named "ir", it must be
    a valid key in ir, and the value of ir[keyname] is
    passed instead of the full ir dict. Invalid keynames
    result in attempts at informative errors, meaning errors
    will often be caught early when making changes.
    """
    def __init__(self, basename):
        ufc_templates = ffc.backends.ufc.templates
        self._header_template = ufc_templates[basename + "_header"]
        self._implementation_template = ufc_templates[basename + "_implementation"]
        self._combined_template = ufc_templates[basename + "_combined"]
        self._jit_header_template = ufc_templates[basename + "_jit_header"]
        self._jit_implementation_template = ufc_templates[basename + "_jit_implementation"]

        r = re.compile(r"%\(([a-zA-Z0-9_]*)\)")
        self._header_keywords = set(r.findall(self._header_template))
        self._implementation_keywords = set(r.findall(self._implementation_template))
        self._combined_keywords = set(r.findall(self._combined_template))

        self._keywords = sorted(self._header_keywords | self._implementation_keywords)

        # Do some ufc interface template checking, to catch bugs
        # early when we change the ufc interface templates
        if set(self._keywords) != set(self._combined_keywords):
            a = set(self._header_keywords) - set(self._combined_keywords)
            b = set(self._implementation_keywords) - set(self._combined_keywords)
            c = set(self._combined_keywords) - set(self._keywords)
            msg = "Templates do not have matching sets of keywords:"
            if a:
                msg += "\n  Header template keywords '%s' are not in the combined template." % (sorted(a),)
            if b:
                msg += "\n  Implementation template keywords '%s' are not in the combined template." % (sorted(b),)
            if c:
                msg += "\n  Combined template keywords '%s' are not in the header or implementation templates." % (sorted(c),)
            error(msg)

    def generate_snippets(self, L, ir, parameters):
        "Generate code snippets for each keyword found in templates."
        snippets = {}
        for kw in self._keywords:
            handlerstr = "%s.%s" % (self.__class__.__name__, kw)

            # Check that attribute self.<keyword> is available
            if not hasattr(self, kw):
                error("Missing handler %s." % (handlerstr,))

            # Call self.<keyword>(*args) to get value to insert in snippets
            method = getattr(self, kw)
            vn = method.__code__.co_varnames[:method.__code__.co_argcount]
            file_line = "%s:%s" % (method.__code__.co_filename, method.__code__.co_firstlineno)

            #if handlerstr == "ufc_dofmap.create":
            #    import ipdb; ipdb.set_trace()

            # Always pass L
            assert vn[:2] == ("self", "L")
            vn = vn[2:]
            args = (L,)

            # Either pass full ir or extract ir value with keyword given by argument name
            if vn[0] == "ir":
                args += (ir,)
            elif vn[0] in ir:
                args += (ir[vn[0]],)
            else:
                error("Cannot find key '%s' in ir, argument to %s at %s." % (vn[0], handlerstr, file_line))
            vn = vn[1:]

            # Optionally pass parameters
            if vn == ("parameters",):
                args += (parameters,)
            elif vn:
                error("Invalid argument names %s to %s at %s." % (vn, handlerstr, file_line))

            # Call handler
            value = method(*args)


            if isinstance(value, list):
                value = L.StatementList(value)

            # Indent body and format to str
            if isinstance(value, L.CStatement):
                value = L.Indented(value.cs_format(precision=parameters["precision"]))
                value = format_indented_lines(value)
            elif not isinstance(value, str):
                error("Expecting code or string, not %s, returned from handler %s at %s." % (type(value), handlerstr, file_line))

            # Store formatted code in snippets dict
            snippets[kw] = value

        # Error checking (can detect some bugs early when changing the ufc interface)
        # Get all attributes of subclass class (skip "_foo")
        attrs = set(name for name in dir(self) if not (name.startswith("_") or name.startswith("generate")))
        # Get all attributes of this base class (skip "_foo" and "generate*")
        base_attrs = set(name for name in dir(ufc_generator) if not (name.startswith("_") or name.startswith("generate")))
        # The template keywords should not contain any names not among the class attributes
        missing = set(self._keywords) - attrs
        if missing:
            warning("*** Missing generator functions:\n%s" % ('\n'.join(map(str, sorted(missing))),))
        # The class attributes should not contain any names not among the template keywords
        # (this is strict, a useful check when changing ufc, but can be dropped)
        unused = attrs - set(self._keywords) - base_attrs
        if unused:
            warning("*** Unused generator functions:\n%s" % ('\n'.join(map(str, sorted(unused))),))

        # Return snippets, a dict of code strings
        return snippets

    def generate(self, L, ir, parameters=None, snippets=None, jit=False):
        "Return composition of templates with generated snippets."
        if snippets is None:
            snippets = self.generate_snippets(L, ir, parameters)
        if jit:
            ht = self._jit_header_template
            it = self._jit_implementation_template
        else:
            ht = self._header_template
            it = self._implementation_template
        h = ht % snippets
        cpp = it % snippets
        return h, cpp

    def classname(self, L, ir):
        "Return classname."
        return ir["classname"]

    def members(self, L, ir):
        "Return empty string. Override in classes that need members."
        if ir.get("members"):
            error("Missing generator function.")
        return ""

    def constructor(self, L, ir):
        "Return empty string. Override in classes that need constructor."
        if ir.get("constructor"):
            error("Missing generator function.")
        return L.NoOp()

    def constructor_arguments(self, L, ir):
        "Return empty string. Override in classes that need constructor."
        if ir.get("constructor_arguments"):
            error("Missing generator function.")
        return ""

    def initializer_list(self, L, ir):
        "Return empty string. Override in classes that need constructor."
        if ir.get("initializer_list"):
            error("Missing generator function.")
        return ""

    def destructor(self, L, ir):
        "Return empty string. Override in classes that need destructor."
        if ir.get("destructor"):
            error("Missing generator function.")
        return L.NoOp()

    def signature(self, L, ir):
        "Default implementation of returning signature string fetched from ir."
        sig = ir["signature"]
        return L.Return(L.LiteralString(sig))

    def create(self, L, ir):
        "Default implementation of creating a new object of the same type."
        classname = ir["classname"]
        return L.Return(L.New(classname))
