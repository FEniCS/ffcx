# Code generation utilities for UFC (Unified Form-assembly Code) v. 1.0.
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenics.org/) 2006.

def indent(block, num_spaces):
    """Indent text block with given number of spaces"""
    return "\n".join([" "*num_spaces + li for li in block.split("\n")])

def generate_code(format_string, format_dictionary, num_indent_spaces=6):
    """Apply default empty strings for unnecessary variables and
    indent contents of format_dictionary"""
    default_dictionary = {'members':'', 'constructor':'', 'destructor':''}
    default_dictionary.update(format_dictionary)
    variables = get_format_variables(format_string)
    for v in variables:
        default_dictionary[v] = indent(default_dictionary[v], num_indent_spaces)
    return format_string % default_dictionary

def __get_format_variables(string)
    vars = []
    pos = string.find("%(")
    while pos >= 0:
        vars.append( string[pos+2:string.find(")s", pos+2)] )
        pos = string.find("%(", pos+1)
    return __make_unique(variables)

def __make_unique(li):
    li.sort()
    i = 0
    while i < len(li)-1:
        if li[i+1] == li[i]:
            li.remove(li[i])
        else:
            i += 1
    return li
