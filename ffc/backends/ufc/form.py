# -*- coding: utf-8 -*-
# Copyright (C) 2009-2017 Anders Logg and Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

# Note: Most of the code in this file is a direct translation from the
# old implementation in FFC

from ffc.backends.ufc import form_template as ufc_form
from ffc.backends.ufc.utils import (generate_return_new,
                                    generate_return_new_switch)
from ffc.representation import ufc_integral_types

# These are the method names in ufc_form that are specialized for each
# integral type
integral_name_templates = (
    "max_{}_subdomain_id",
    "has_{}_integrals",
    "create_{}_integral",
    "create_default_{}_integral",
)


def create_delegate(integral_type, declname, impl):
    def _delegate(self, L, ir, parameters):
        return impl(self, L, ir, parameters, integral_type, declname)

    _delegate.__doc__ = impl.__doc__ % {"declname": declname, "integral_type": integral_type}
    return _delegate


def add_ufc_form_integral_methods(cls):
    """Generate methods on the class decorated by this function,
    for each integral name template and for each integral type.

    This allows implementing e.g. create_###_integrals once in the
    decorated class as '_create_foo_integrals', and this function will
    expand that implementation into 'create_cell_integrals',
    'create_exterior_facet_integrals', etc.

    Name templates are taken from 'integral_name_templates' and 'ufc_integral_types'.
    """
    # The dummy name "foo" is chosen for familiarity for ffc developers
    dummy_integral_type = "foo"

    for template in integral_name_templates:
        implname = "_" + (template.format(dummy_integral_type))
        impl = getattr(cls, implname)
        for integral_type in ufc_integral_types:
            declname = template.format(integral_type)
            _delegate = create_delegate(integral_type, declname, impl)
            setattr(cls, declname, _delegate)
    return cls


@add_ufc_form_integral_methods
class UFCForm:
    """Each function maps to a keyword in the template.

    The exceptions are functions on the form
        def _*_foo_*(self, L, ir, parameters, integral_type, declname)
    which add_ufc_form_integral_methods will duplicate for foo = each integral type.
    """

    def num_coefficients(self, L, num_coefficients):
        return L.Return(num_coefficients)

    def rank(self, L, rank):
        return L.Return(rank)

    def original_coefficient_position(self, L, ir):
        i = L.Symbol("i")
        positions = ir["original_coefficient_position"]

        # Check argument
        msg = "Invalid original coefficient index."

        if positions:
            code = [L.If(L.GE(i, len(positions)), [L.Comment(msg), L.Return(-1)])]

            position = L.Symbol("position")
            code += [
                L.ArrayDecl("static const int64_t", position, len(positions), positions),
                L.Return(position[i]),
            ]
            return code
        else:
            code = [L.Comment(msg), L.Return(-1)]
        return code

    def create_coordinate_finite_element(self, L, ir):
        print('create coord')
        classnames = ir["create_coordinate_finite_element"]
        assert len(classnames) == 1
        return generate_return_new(L, classnames[0], factory=True)

    def coordinate_finite_element_declaration(self, L, ir):
        print('create coord decl')
        classname = ir["create_coordinate_finite_element"]
        code = "ufc_finite_element* create_{name}();\n".format(name=classname[0])
        return code

    def create_coordinate_dofmap(self, L, ir):
        classnames = ir["create_coordinate_dofmap"]
        assert len(classnames) == 1
        return generate_return_new(L, classnames[0], factory=True)

    def coordinate_dofmap_declaration(self, L, ir):
        classname = ir["create_coordinate_dofmap"]
        code = "ufc_dofmap* create_{name}();\n".format(name=classname[0])
        return code

    def create_coordinate_mapping(self, L, ir):
        classnames = ir["create_coordinate_mapping"]
        # list of length 1 until we support multiple domains
        assert len(classnames) == 1
        return generate_return_new(L, classnames[0], factory=True)

    def coordinate_mapping_declaration(self, L, ir):
        classname = ir["create_coordinate_mapping"]
        code = "ufc_coordinate_mapping* create_{name}();\n".format(name=classname[0])
        return code

    def create_finite_element(self, L, ir):
        i = L.Symbol("i")
        classnames = ir["create_finite_element"]
        return generate_return_new_switch(L, i, classnames, factory=True)

    def finite_element_declaration(self, L, ir):
        classnames = set(ir["create_finite_element"])
        code = ""
        for name in classnames:
            code += "ufc_finite_element* create_{name}();\n".format(name=name)
        return code

    def create_dofmap(self, L, ir):
        i = L.Symbol("i")
        classnames = ir["create_dofmap"]
        return generate_return_new_switch(L, i, classnames, factory=True)

    def dofmap_declaration(self, L, ir):
        classnames = set(ir["create_dofmap"])
        code = ""
        for name in classnames:
            code += "ufc_dofmap* create_{name}();\n".format(name=name)
        return code

    # This group of functions are repeated for each
    # foo_integral by add_ufc_form_integral_methods:

    def _max_foo_subdomain_id(self, L, ir, parameters, integral_type, declname):
        """Return implementation of ufc::form::%(declname)s()."""
        # e.g. max_subdomain_id = ir["max_cell_subdomain_id"]
        max_subdomain_id = ir[declname]
        return L.Return(int(max_subdomain_id))

    def _has_foo_integrals(self, L, ir, parameters, integral_type, declname):
        """Return implementation of ufc::form::%(declname)s()."""
        # e.g. has_integrals = ir["has_cell_integrals"]
        has_integrals = ir[declname]
        return L.Return(bool(has_integrals))

    def _create_foo_integral(self, L, ir, parameters, integral_type, declname):
        """Return implementation of ufc::form::%(declname)s()."""
        # e.g. subdomain_ids, classnames = ir["create_cell_integral"]
        subdomain_ids, classnames = ir[declname]
        subdomain_id = L.Symbol("subdomain_id")
        return generate_return_new_switch(L, subdomain_id, classnames, subdomain_ids, factory=True)

    def _create_default_foo_integral(self, L, ir, parameters, integral_type, declname):
        """Return implementation of ufc::form::%(declname)s()."""
        # e.g. classname = ir["create_default_cell_integral"]
        classname = ir[declname]
        if classname is None:
            return L.Return(L.Null())
        else:
            return generate_return_new(L, classname, factory=True)


def ufc_form_generator(ir, parameters):
    """Generate UFC code for a form"""

    factory_name = ir["classname"]

    d = {}
    d["factory_name"] = factory_name
    d["signature"] = "\"{}\"".format(ir["signature"])
    d["rank"] = ir["rank"]
    d["num_coefficients"] = ir["num_coefficients"]

    d["max_cell_subdomain_id"] = ir["max_cell_subdomain_id"]
    d["max_exterior_facet_subdomain_id"] = ir["max_exterior_facet_subdomain_id"]
    d["max_interior_facet_subdomain_id"] = ir["max_interior_facet_subdomain_id"]
    d["max_vertex_subdomain_id"] = ir["max_vertex_subdomain_id"]
    d["max_custom_subdomain_id"] = ir["max_custom_subdomain_id"]

    d["has_cell_integrals"] = "true" if ir["has_cell_integrals"] else "false"
    d["has_exterior_facet_integrals"] = "true" if ir["has_exterior_facet_integrals"] else "false"
    d["has_interior_facet_integrals"] = "true" if ir["has_interior_facet_integrals"] else "false"
    d["has_vertex_integrals"] = "true" if ir["has_vertex_integrals"] else "false"
    d["has_custom_integrals"] = "true" if ir["has_custom_integrals"] else "false"

    import ffc.uflacs.language.cnodes as L
    generator = UFCForm()

    statements = generator.original_coefficient_position(L, ir)
    assert isinstance(statements, list)
    d["original_coefficient_position"] = L.StatementList(statements)

    d["create_coordinate_finite_element"] = generator.create_coordinate_finite_element(L, ir)
    d["coordinate_finite_element_declaration"] = generator.coordinate_finite_element_declaration(L, ir)
    d["create_coordinate_dofmap"] = generator.create_coordinate_dofmap(L, ir)
    d["coordinate_dofmap_declaration"] = generator.coordinate_dofmap_declaration(L, ir)
    d["create_coordinate_mapping"] = generator.create_coordinate_mapping(L, ir)
    d["coordinate_mapping_declaration"] = generator.coordinate_mapping_declaration(L, ir)
    d["create_finite_element"] = generator.create_finite_element(L, ir)
    d["finite_element_declaration"] = generator.finite_element_declaration(L, ir)
    d["create_dofmap"] = generator.create_dofmap(L, ir)
    d["dofmap_declaration"] = generator.dofmap_declaration(L, ir)

    d["create_cell_integral"] = generator.create_cell_integral(L, ir, parameters)
    d["create_interior_facet_integral"] = generator.create_interior_facet_integral(
        L, ir, parameters)
    d["create_exterior_facet_integral"] = generator.create_exterior_facet_integral(
        L, ir, parameters)
    d["create_vertex_integral"] = generator.create_vertex_integral(L, ir, parameters)
    d["create_custom_integral"] = generator.create_custom_integral(L, ir, parameters)

    d["create_default_cell_integral"] = generator.create_default_cell_integral(L, ir, parameters)
    d["create_default_interior_facet_integral"] = generator.create_default_interior_facet_integral(
        L, ir, parameters)
    d["create_default_exterior_facet_integral"] = generator.create_default_exterior_facet_integral(
        L, ir, parameters)
    d["create_default_vertex_integral"] = generator.create_default_vertex_integral(
        L, ir, parameters)
    d["create_default_custom_integral"] = generator.create_default_custom_integral(
        L, ir, parameters)

    # Check that no keys are redundant or have been missed
    from string import Formatter
    fields = [fname for _, fname, _, _ in Formatter().parse(ufc_form.factory) if fname]
    assert set(fields) == set(d.keys()), "Mismatch between keys in template and in formattting dict"

    # Format implementation code
    implementation = ufc_form.factory.format_map(d)

    # Format declaration
    declaration = ufc_form.declaration.format(factory_name=factory_name)

    return declaration, implementation
