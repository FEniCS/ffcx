__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2005-05-20 -- 2009-03-15"
__copyright__ = "Copyright (C) 2005-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

FFC_VERSION = "0.6.1"

FFC_DEBUG_LEVEL = 0

FFC_OPTIONS = {"representation": "tensor",
               "language": "ufc",
               "format": "ufc",
               "optimize": False,
               "cpp optimize": False,
               "blas": False,
               "precision": "15",
               "quadrature_points": False,
               "split_implementation": False,
               "form_postfix": True,
               "cache_dir": None,
               "output_dir": ".",
               "external_signature": None,
               "compiler": "ffc",
               "quadrature_order": "auto"}

# FIXME: New options to replace FFC options
UFL_OPTIONS = {"representation": "auto",      # form representation / code generation strategy
               "format": "ufc",               # code generation format
               "quadrature_order": "auto",    # quadrature order used for quadrature representation
               "precision": "15",             # precision used when writing numbers
               "split_implementation": False, # split generated code into .h and .cpp file
               "form_postfix": True,          # postfix form name with "Function", "LinearForm" or BilinearForm
               "cache_dir": None,             # cache dir used by Instant
               "output_dir": "."}             # output directory for generated code
