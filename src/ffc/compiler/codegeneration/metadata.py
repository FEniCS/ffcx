"Code generation for meta data"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-02-28 -- 2007-02-28"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

def generate_meta_data(form_data, format):
    "Set form meta data"

    code = {}

    # Set name of form
    code["name"] = form_data.name

    # Set number of arguments
    code["num_arguments"] = form_data.num_arguments

    # Set shape dimension (should be the same so pick first)
    code["shape_dimension"] = form_data.elements[0].shape_dimension()

    return code
