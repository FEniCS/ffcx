"Raw output format."

__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-11-17 -- 2006-04-01"
__copyright__ = "Copyright (C) 2004-2006 Anders Logg"
__license__  = "GNU GPL Version 2"

# Specify formatting for code generation
format = { "sum": lambda l: " + ".join(l),
           "subtract": lambda l: " - ".join(l),
           "multiplication": lambda l: "*".join(l),
           "grouping": lambda s: "(%s)" % s,
           "determinant": None,
           "floating point": lambda a: "%.15e" % a,
           "constant": lambda j: None,
           "coefficient": lambda j, k: None,
           "coefficient table": lambda j, k: None,
           "transform": lambda j, k: None,
           "reference tensor" : lambda j, i, a: "(%d, %s, %s)" % (j, str(i), (str(a))),
           "geometry tensor": lambda j, a: None,
           "element tensor": lambda i, k: None,
           "tmp declaration": lambda j, k: None,
           "tmp access": lambda j, k: None }

def init(options):
    "Initialize code generation for raw format."
    return

def write(forms, options):
    "Generate code for raw format."
    print "Generating raw output"

    for j in range(len(forms)):

        # Write form
        output = ""
        output += __form(forms[j])
        
        # Write file
        if len(forms) > 1: extra = "-%d" % j
        else: extra = ""
        filename = "%s%s.raw" % (forms[j].name, extra)
        file = open(filename, "w")
        file.write(output)
        file.close()
        print "Output written to " + filename

    return

def __form(form):
    "Generate form in raw format."
    output = ""    

    # Interior contribution
    if form.AKi.terms:
        for a0 in form.AKi.a0:
            output += "interior %s %s\n" % (a0.name, a0.value)

    # Boundary contribution
    if form.AKb.terms:
        for g0 in form.AKi.a0:
            output += "boundary %s %s\n" % (a0.name, a0.value)

    return output
