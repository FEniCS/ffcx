"Raw output format."

__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2005-09-29 -- 2006-04-01"
__copyright__ = "Copyright (C) 2005-2006 Anders Logg"
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
           "reference tensor" : lambda j, i, a: None,
           "geometry tensor": lambda j, a: None,
           "element tensor": lambda i, k: None,
           "tmp declaration": lambda j, k: None,
           "tmp access": lambda j, k: None }

def init(options):
    "Initialize code generation for XML format."
    return

def write(forms, options):
    "Generate code for XML format."
    print "Generating XML output"

    for j in range(len(forms)):

        # Generate name of form
        if len(forms) > 1:
            name = "%s-%d" % (forms[j].name, j)
        else:
            name = forms[j].name

        # Open file
        filename = "%s.xml" % name
        file = open(filename, "w")        

        # Write form
        __header(file, forms[j])
        __form(file, forms[j], name)
        __footer(file, forms[j])

        # Close file
        file.close()
        print "Output written to " + filename

    return

def __header(file, form):
    "Generate header in XML format."
    file.write("""\
<?xml version="1.0" encoding="UTF-8"?>

<ffc xmlns:ffc="http://www.fenics.org/ffc/">
""")

def __footer(file, form):
    "Generate footer in XML format."
    file.write("""\
</ffc>
""")

def __form(file, form, name):
    "Generate form in XML format."
    file.write("  <form name=\"%s\">\n" % name)
    
    # Interior contribution
    if len(form.AKi.terms) > 0:
        file.write("    <interior>\n")
        for term in form.AKi.terms:

            # Extract data for term
            sig = term.signature()
            A0 = term.A0
            GK = term.GKs[0] # pick first, all should be the same
            rank_a = A0.rank
            rank_g = GK.rank
            size_a = len(A0.i.indices)
            size_g = len(GK.a.indices)

            # Write data for term
            file.write("""\
      <term signature=\"%s\" size=\"%s\">
        <geometrytensor rank=\"%d\" size=\"%d\"></geometrytensor>
        <referencetensor rank=\"%d\">
""" % (sig, size_a, rank_g, size_g, rank_a))
            iindices = A0.i.indices
            aindices = A0.a.indices or [[]]
            for i in iindices:
                for a in aindices:
                    index = " ".join([str(j) for j in (i + a)])
                    value = A0.A0[i + a]
                    file.write("          <entry index=\"%s\" value=\"%s\"/>\n" % (index, value))
            file.write("""\
        </referencetensor>
      </term>
""")
        file.write("    </interior>\n")

    file.write("  </form>\n")
