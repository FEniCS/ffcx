#!/usr/bin/env python

import hotshot, hotshot.stats, test.pystone

# Python modules
import sys
import getopt
import os.path

# UFL modules
from ufl.log import UFLException

# FFC modules
from ffc.common.log import info, set_level, DEBUG, ERROR
from ffc.common.constants import FFC_OPTIONS, FFC_VERSION

def error(msg):
    "Print error message (cannot use log system at top level)."
    print("\n".join(["*** FFC: " + line for line in msg.split("\n")]))

def info_version():
    "Print version number."
    info("""\
This is FFC, the FEniCS Form Compiler, version %s.
For further information, visit http://www.fenics.org/ffc/.
""" % FFC_VERSION)

def info_usage():
    "Print usage information."
    info_version()
    info("""Usage: ffc [OPTION]... input.form

For information about the FFC command-line interface, refer to
the FFC man page which may invoked by 'man ffc' (if installed).
""")

def main():
    "Main function."
    argv = sys.argv[1:]
    # Get command-line arguments
    try:
        opts, args = getopt.getopt(argv, \
        "hvdsl:r:f:Oo:q:", \
        ["help", "version", "debug", "silent", "language=", "representation=",
         "optimize", "output-directory=", "quadrature-rule="])
    except getopt.GetoptError:
        info_usage()
        error("Illegal command-line arguments.")
        return 1

    # Check for --help
    if ("-h", "") in opts or ("--help", "") in opts:
        info_usage()
        return 0

    # Check for --version
    if ("-v", "") in opts or ("--version", "") in opts:
        info_version()
        return 0

    # Check that we get at least one file
    if len(args) == 0:
        error("Missing file.")
        return 1

    # Parse command-line options
    options = FFC_OPTIONS.copy()
    for opt, arg in opts:
        if opt in ("-d", "--debug"):
             options["log_level"] = DEBUG
        elif opt in ("-s", "--silent"):
            options["log_level"] = ERROR
        elif opt in ("-l", "--language"):
            options["format"] = arg
        elif opt in ("-r", "--representation"):
            options["representation"] = arg
        elif opt in ("-q", "--quadrature-rule"):
            options["quadrature_rule"] = arg
        elif opt == "-f":
            if len(arg.split("=")) == 2:
                (key, value) = arg.split("=")
                options[key] = value
            elif len(arg.split("==")) == 1:
                key = arg.split("=")[0]
                options[arg] = True
            else:
                info_usage()
                return 1
        elif opt in ("-O", "--optimize"):
            options["optimize"] = True
        elif opt in ("-o", "--output-directory"):
            options["output_dir"] = arg

    # Set log level
    set_level(options["log_level"])

    # Print a nice message
    info_version()

    # Call parser and compiler for each file
    for filename in args:

        # Get filename suffix
        suffix = filename.split(".")[-1]

        # Check file suffix and parse file/generate module
        if suffix == "ufl":
            script = _make_script_ufl(filename, options)
        elif suffix == "form":
            error("Old style .form files are no longer supported. Use form2ufl to convert to UFL format.")
            return 1
        else:
            error("Expecting a UFL form file (.ufl).")
            return 1

        # Catch exceptions only when not in debug mode
        if options["log_level"] <= DEBUG:
            exec(compile(open(script).read(), script, 'exec'), {})
        else:
            try:
                exec(compile(open(script).read(), script, 'exec'), {})
            except UFLException as exception:
                info("")
                error(str(exception))
                error("To get more information about this error, rerun FFC with --debug.")
                return 1

    return 0

# New version for .ufl files
def _make_script_ufl(filename, options):
    "Create Python script from given .ufl file and return name of script"

    # FIXME: Use load_forms() from UFL here instead

    # Strip directories from filename for output file
    filename_out = filename.split(os.path.sep)[-1]

    # Get prefix of file name and generate Python script file name
    prefix = ".".join(filename_out.split(".")[:-1])
    script = prefix + ".py"
    info("Preprocessing form file: %s --> %s\n" % (filename, script))

    # Read input
    infile = open(filename, "r")
    input = infile.read()
    infile.close()

    # Generate output
    output = """\
from ufl import *
from ufl.log import set_level
from ffc.compiler.compiler import compile

# Set debug level
set_level(%d)

# Reserved variables for forms
(a, L, M) = (None, None, None)

# Reserved variable for element
element = None

%s
compile([a, L, M, element], \"%s\", %s, globals())
""" % (options["log_level"], input, prefix, options)

    # Write output
    outfile = open(script, "w")
    outfile.write(output)
    outfile.close()

    # Return script filename
    return script


if __name__ == "__main__":
#    sys.exit(main(sys.argv[1:]))

    name = "%s.prof" % sys.argv[-1].split(".")[0].lower()
#    print name
    prof = hotshot.Profile(name)
#    prof.runcall()
    prof.runcall(main)
    prof.close()
    stats = hotshot.stats.load(name)
    stats.strip_dirs()
    stats.sort_stats("time").print_stats(50)

#    stats.sort_stats('time', 'calls')
#    stats.print_stats(20)


