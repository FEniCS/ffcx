__author__ = "Anders Logg, Kristian B. Oelgaard and Marie E. Rognes"
__date__ = "2010-01-21"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from ffc.log import begin, end, info, info_red, info_green, info_blue
import os, sys, shutil, commands

# Parameters
output_directory = "output"

def generate_tests():
    "Generate form files for all test cases."

    begin("Generating test cases")

    # Create output directory if it does not exist
    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    # Copy demos files
    demo_files = [f for f in os.listdir("../../demo/") if f.endswith(".ufl")]
    for f in demo_files:
        shutil.copy("../../demo/%s" % f, output_directory)
    info_green("Copied %d demo files" % len(demo_files))

    end()

def generate_code():
    "Generate code for all test cases."

    begin("Generating code")

    # Get a list of all files
    form_files = [f for f in os.listdir(output_directory)]

    # Iterate over all files
    for f in form_files:

        # Generate code
        info("Generating code for %s" % f)
        status, output = commands.getstatusoutput("ffc %s/%s" % (output_directory, f))

        # Check status
        if status == 0:
            info_green("OK")
        else:
            info_red("Failed")

    end()

def main(args):
    "Run all regression tests."

    # Generate form files for all test cases
    generate_tests()

    # Generate code for all test cases
    generate_code()

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
