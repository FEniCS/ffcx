# Copyright (C) 2009-2018 Anders Logg, Martin Sandve Alnæs and Garth N. Wells
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# Note: Most of the code in this file is a direct translation from the
# old implementation in FFC

import logging
import pprint
import textwrap

from ffcx.codegeneration.C import file_template
from ffcx import __version__ as FFCX_VERSION
from ffcx.codegeneration import __version__ as UFC_VERSION


logger = logging.getLogger("ffcx")


def generator(options):
    """Generate UFC code for file output."""
    logger.info("Generating code for file")

    # Attributes
    d = {"ffcx_version": FFCX_VERSION, "ufcx_version": UFC_VERSION}
    d["options"] = textwrap.indent(pprint.pformat(options), "//  ")
    extra_c_includes = []
    if "_Complex" in options["scalar_type"]:
        extra_c_includes += ["complex.h"]
    d["extra_c_includes"] = "\n".join(
        f"#include <{header}>" for header in extra_c_includes
    )

    # Format declaration code
    code_pre = (
        file_template.declaration_pre.format_map(d),
        file_template.implementation_pre.format_map(d),
    )

    # Format implementation code
    code_post = (
        file_template.declaration_post.format_map(d),
        file_template.implementation_post.format_map(d),
    )

    return code_pre, code_post
