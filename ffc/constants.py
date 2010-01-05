__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2005-05-20"
__copyright__ = "Copyright (C) 2005-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard, 2009
# Last changed: 2010-01-04

from log import INFO

FFC_VERSION = "0.7.1"

FFC_OPTIONS = {"representation":      "auto", # form representation / code generation strategy
               "format":              "ufc",  # code generation format
               "quadrature_degree":   "auto", # quadrature degree used for computing integrals
               "quadrature_rule":      None,  # quadrature rule used for integration of element tensors
               "precision":           "15",   # precision used when writing numbers
               "split":                False, # split generated code into .h and .cpp file
               "form_postfix":         True,  # postfix form name with "Function", "LinearForm" or BilinearForm
               "cache_dir":            None,  # cache dir used by Instant
               "output_dir":           ".",   # output directory for generated code
               "cpp optimize":         False, # optimization for the JIT compiler
               "optimize":             False, # optimise the quadrature code generation
               "log_level":            INFO,  # log level, displaying only messages with level >= log_level
               "log_prefix":           "",    # log prefix
               "epsilon":              1e-14} # machine precision, used for dropping zero terms
