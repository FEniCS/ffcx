"LaTeX output format."

__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-10-14"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

format = { "multiplication": "",
           "determinant": "\\det F_K'",
           "floating point": lambda a: "%.2f" % a,
           "coefficient": lambda j, k: "c_{%d%d}" % (j, k),
           "transform": lambda j, k: "\\frac{\\partial X_{%d}}{\\partial x_{%d}}" % (j, k),
           "geometry tensor": lambda j, a: "G_{K,%d}^{%s}" % (j, "".join(["%d" % index for index in a])),
           "element tensor": lambda i: "A^K_{%s}" % "".join(["%d" % index for index in i]) }

def compile(form):
    "Generate code for LaTeX."
    print "Compiling form %s for LaTeX" % str(form)

    # Generate output
    output = ""
    output += __file_header()
    output += __form(form)
    output += __file_footer()

    # Write file
    filename = form.name + ".tex"
    file = open(filename, "w")
    file.write(output)
    file.close()
    print "Output written on " + filename

    return

def __file_header():
    "Generate file header for LaTeX."
    return """
\\documentclass[12pt]{article}

\\begin{document}
"""

def __file_footer():
    "Generate file footer for LaTeX."
    return """
\\end{document}
"""

def __form(form):
    "Generate form for LaTeX."
    output = ""    

    # Interior contribution
    if form.AKi.terms:
        output += """
\\subsection*{Interior contribution}

\\begin{equation}
  \\begin{array}{rcl}
"""
        for gK in form.AKi.gK:
            output += "    %s &=& %s \\\\\n" % (gK.name, gK.value)
        output += """  \\end{array}
\\end{equation}

\\begin{equation}
  \\begin{array}{rcl}
"""
        for aK in form.AKi.aK:
            output += "    %s &=& %s \\\\\n" % (aK.name, aK.value)
        output += """ \\end{array}
\\end{equation}
"""

    return output
