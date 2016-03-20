# Copyright (C) 2005-2015 Anders Logg
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

from ffc.log import INFO


# Comments from other places in code:
# FIXME: Document option -fconvert_exceptions_to_warnings
# FIXME: Remove option epsilon and just rely on precision?


FFC_PARAMETERS = {
  "format":                         "ufc",   # code generation format
  "representation":                 "auto",  # form representation / code
                                             # generation strategy
  "quadrature_rule":                "auto",  # quadrature rule used for
                                             # integration of element tensors
  "quadrature_degree":              -1,      # quadrature degree used for
                                             # computing integrals
  "precision":                      15,      # precision used when writing
                                             # numbers
  "epsilon":                        1e-14,   # machine precision, used for
                                             # dropping zero terms
  "split":                          False,   # split generated code into .h and
                                             # .cpp file
  "form_postfix":                   True,    # postfix form name with "Function",
                                             # "LinearForm" or BilinearForm
  "convert_exceptions_to_warnings": False,   # convert all exceptions to warning
                                             # in generated code
  "cache_dir":                      "",      # cache dir used by Instant
  "output_dir":                     ".",     # output directory for generated
                                             # code
  "cpp_optimize":                   True,    # optimization for the JIT compiler
  "cpp_optimize_flags":             "-O2",   # optimization flags for the JIT compiler
  "optimize":                       False,   # optimise the code generation
  "log_level":                      INFO+5,  # log level, displaying only
                                             # messages with level >= log_level
  "log_prefix":                     "",      # log prefix
  "error_control":                  False,   # with error control
}


def default_parameters():
    "Return (a copy of) the default parameter values for FFC."
    parameters = FFC_PARAMETERS.copy()

    # HACK
    import os
    r = os.environ.get("FFC_FORCE_REPRESENTATION")
    if r: parameters["representation"] = r

    return parameters


def default_jit_parameters():
    parameters = default_parameters()
    parameters["no-evaluate_basis_derivatives"] = True
    return parameters


def compilation_relevant_parameters(parameters):
    parameters = parameters.copy()
    ignores = ["log_prefix", "log_level", "cache_dir", "output_dir"]
    for ignore in ignores:
        assert ignore in FFC_PARAMETERS
        if ignore in parameters:
            del parameters[ignore]

    # HACK
    import os
    r = os.environ.get("FFC_FORCE_REPRESENTATION")
    if r: parameters["representation"] = r

    return parameters
