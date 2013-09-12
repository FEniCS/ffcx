# Copyright (C) 2005-2010 Anders Logg
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
# First added:  2005-05-20
# Last changed: 2005-05-20

from log import INFO

# Last changed: 2011-01-18

FFC_PARAMETERS = {
  "format":                         "ufc",   # code generation format
  "representation":                 "auto",  # form representation / code
                                             # generation strategy
  "quadrature_rule":                "auto",  # quadrature rule used for
                                             # integration of element tensors
  "quadrature_degree":              "auto",  # quadrature degree used for
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
  "restrict_keyword":               "",      # compiler specific "__restrict" or "__restrict__" keyword
  "log_level":                      INFO,    # log level, displaying only
                                             # messages with level >= log_level
  "log_prefix":                     "",      # log prefix
  "error_control":                  False,   # with error control
}

def default_parameters():
    "Return (a copy of) the default parameter values for FFC."
    return FFC_PARAMETERS.copy()
