__author__ = "Anders Logg, Kristian B. Oelgaard and Marie E. Rognes"
__date__ = "2010-01-21"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from ffc.log import begin, end, info, info_red, info_green, info_blue
import os, sys, shutil, commands

# Parameters
output_directory = "output"

def generate_test_cases():
    "Generate form files for all test cases."

    begin("Generating test cases")

    # Copy demos files
    demo_files = [f for f in os.listdir("../../../demo/") if f.endswith(".ufl")]
    for f in demo_files:
        shutil.copy("../../../demo/%s" % f, ".")
    info_green("Copied %d demo files" % len(demo_files))

    # Generate form files for forms
    info("Not implemented")

    # Generate form files for elements
    info("Not implemented")

    end()

def generate_code():
    "Generate code for all test cases."

    begin("Generating code")

    # Get a list of all files
    form_files = [f for f in os.listdir(".") if f.endswith(".ufl")]

    # Iterate over all files
    for f in form_files:

        # Generate code
        status, output = commands.getstatusoutput("ffc %s" % f)

        # Check status
        if status == 0:
            info_green("%s OK" % f)
        else:
            info_red("%s failed" % f)

    end()

def validate_code():
    "Validate generated code against references."

    begin("Validating generated code")

    # Get a list of all files
    header_files = [f for f in os.listdir(".") if f.endswith(".h")]

    # Iterate over all files
    for f in header_files:

        # Get generated code
        generated_code = open(f).read()

        # Get reference code
        reference_file = "../reference/%s" % f
        if os.path.isfile(reference_file):
            reference_code = open(reference_file).read()
        else:
            info_blue("Missing reference for %s" % f)

    end()

def build_programs():
    "Build test programs for all test cases."

    begin("Building test programs")

    # Get a list of all files
    header_files = [f for f in os.listdir(".") if f.endswith(".h")]

    # Iterate over all files
    for f in header_files:
        info("Building test for %s" % f)

    end()

def validate_programs():
    "Validate generated programs against references."

    begin("Validating generated programs")

    info("Not implemented")

    end()

def main(args):
    "Run all regression tests."

    # Enter output directory
    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)
    os.chdir(output_directory)

    # Generate test cases
    generate_test_cases()

    # Generate and validate code
    #generate_code()
    validate_code()

    # Build and validate programs
    build_programs()
    validate_programs()

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
