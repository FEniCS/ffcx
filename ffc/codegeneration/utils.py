# -*- coding: utf-8 -*-
# Copyright (C) 2015-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

# TODO: Move these to ffc.language utils?


def generate_return_new(L, classname):
    return L.Return(L.Call("create_" + classname))


def generate_return_new_switch(L, i, classnames, args=None):

    if isinstance(i, str):
        i = L.Symbol(i)

    def create(classname):
        return L.Call("create_" + classname)

    default = L.Return(L.Null())
    if classnames:
        cases = []
        if args is None:
            args = list(range(len(classnames)))
        for j, classname in zip(args, classnames):
            if classname:
                cases.append((j, L.Return(create(classname))))
        return L.Switch(i, cases, default=default)
    else:
        return default


def generate_return_literal_switch(L,
                                   i,
                                   values,
                                   default,
                                   literal_type,
                                   typename=None):
    # TODO: UFC functions of this type could be replaced with return vector<T>{values}.

    if isinstance(i, str):
        i = L.Symbol(i)
    return_default = L.Return(literal_type(default))

    if values and typename is not None:
        # Store values in static table and return from there
        V = L.Symbol("return_values")
        decl = L.ArrayDecl("static const %s" % typename, V, len(values),
                           [literal_type(k) for k in values])
        return L.StatementList(
            [decl,
             L.If(L.GE(i, len(values)), return_default),
             L.Return(V[i])])
    elif values:
        # Need typename to create static array, fallback to switch
        cases = [(j, L.Return(literal_type(k))) for j, k in enumerate(values)]
        return L.Switch(i, cases, default=return_default)
    else:
        # No values, just return default
        return return_default


def generate_return_int_switch(L, i, values, default):
    return generate_return_literal_switch(L, i, values, default, L.LiteralInt,
                                          "int")
