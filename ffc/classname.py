# -*- coding: utf-8 -*-
# Copyright (C) 2009-2017 Anders Logg
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""This module defines some basics for generating C code."""

def make_name(prefix, basename, signature):
    pre = prefix.lower() + "_" if prefix else ""
    sig = str(signature).lower()
    return "{}{}_{}".format(pre, basename, sig)


def make_integral_name(prefix, integral_type, form_id, subdomain_id):
    basename = "{}_integral_{}".format(integral_type, str(form_id).lower())
    return make_name(prefix, basename, subdomain_id)
