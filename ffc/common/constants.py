__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2005-05-20 -- 2009-04-21"
__copyright__ = "Copyright (C) 2005-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

from log import INFO

FFC_VERSION = "0.6.2"

FFC_OPTIONS = {"debug_level":         INFO,   # debug level
               "representation":      "auto", # form representation / code generation strategy
               "format":              "ufc",  # code generation format
               "quadrature_order":    "auto", # quadrature order used for quadrature representation
               "precision":           "15",   # precision used when writing numbers
               "split_implementation": False, # split generated code into .h and .cpp file
               "form_postfix":         True,  # postfix form name with "Function", "LinearForm" or BilinearForm
               "cache_dir":            None,  # cache dir used by Instant
               "output_dir":           ".",   # output directory for generated code
               "cpp optimize":         False, # optimization for the JIT compiler
               "optimize":             False} # optimise the quadrature code generation
