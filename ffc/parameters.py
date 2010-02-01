__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2005-05-20"
__copyright__ = "Copyright (C) 2005-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from log import INFO

FFC_PARAMETERS = {"format":                         "ufc",  # code generation format
                  "representation":                 "auto", # form representation / code generation strategy
                  "quadrature_rule":                "auto", # quadrature rule used for integration of element tensors
                  "quadrature_degree":              "auto", # quadrature degree used for computing integrals
                  "precision":                      15,     # precision used when writing numbers
                  "epsilon":                        1e-14,  # machine precision, used for dropping zero terms
                  "split":                          False,  # split generated code into .h and .cpp file
                  "form_postfix":                   True,   # postfix form name with "Function", "LinearForm" or BilinearForm
                  "convert_exceptions_to_warnings": False,  # convert all exceptions to warning in generated code
                  "cache_dir":                      "",     # cache dir used by Instant
                  "output_dir":                     ".",    # output directory for generated code
                  "cpp_optimize":                   False,  # optimization for the JIT compiler
                  "optimize":                       False,  # optimise the quadrature code generation
                  "log_level":                      INFO,   # log level, displaying only messages with level >= log_level
                  "log_prefix":                     ""}     # log prefix

def default_parameters():
    "Return (a copy of) the default parameter values for FFC."
    return FFC_PARAMETERS.copy()
