"Raw output format."

__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2005-09-29"
__copyright__ = "Copyright (c) 2005 Anders Logg"
__license__  = "GNU GPL Version 2"

format = { "sum": lambda l: " + ".join(l),
           "subtract": lambda l: " - ".join(l),
           "multiplication": lambda l: "*".join(l),
           "grouping": lambda s: "(%s)" % s,
           "determinant": "not defined",
           "floating point": lambda a: "%.15e" % a,
           "constant": lambda j: "not defined",
           "coefficient": lambda j, k: "not defined",
           "transform": lambda j, k: "not defined",
           "reference tensor" : lambda j, i, a: "not defined",
           "geometry tensor": lambda j, a: "not defined",
           "element tensor": lambda i, k: "not defined" }

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

        # Write form
        output  = __header(forms[j])
        output += __form(forms[j], name)
        output += __footer(forms[j])

        # Write file
        filename = "%s.xml" % name
        file = open(filename, "w")
        file.write(output)
        file.close()
        print "Output written to " + filename

    return

def __header(form):
    "Generate header in XML format."
    return """\
<?xml version="1.0" encoding="UTF-8"?>

<ffc xmlns:ffc="http://www.fenics.org/ffc/">
"""

def __footer(form):
    "Generate footer in XML format."
    return """\
</ffc>
"""

def __form(form, name):
    "Generate form in XML format."
    output = "  <form name=\"%s\">\n" % name
    
    # Interior contribution
    if len(form.AKi.terms) > 0:
        output += "    <interior>\n"
        for term in form.AKi.terms:

            # Extract data for term
            sig = term.signature()
            A0 = term.A0
            rank_a = A0.rank
            rank_g = term.GKs[0].rank

            # Write data for term
            output += "      <term signature=\"%s\">\n" % sig
            output += "        <geometrytensor rank=\"%d\"></geometrytensor>\n" % rank_g
            output += "        <referencetensor rank=\"%d\">\n" % rank_a
            iindices = A0.i.indices
            aindices = A0.a.indices or [[]]
            for i in iindices:
                for a in aindices:
                    index = " ".join([str(j) for j in (i + a)])
                    value = A0.A0[i + a]
                    output += "          <entry index=\"%s\" value=\"%s\"/>\n" % (index, value)
            output += "        </referencetensor>\n"
            output += "      </term>\n"
        output += "    </interior>\n"

    # Boundary contribution
    if form.AKb.terms:
        for g0 in form.AKi.a0:
            output += "    boundary %s %s\n" % (a0.name, a0.value)

    output += "  </form>\n"

    return output
