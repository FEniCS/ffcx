__author__ = "Anders Logg, Kristian B. Oelgaard and Marie E. Rognes"
__date__ = "2010-01-21"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from ffc.log import begin, end, info, info_red, info_green, info_blue
import os, sys, shutil, commands

from ufctest import form_test_code, element_test_code

# Parameters
output_directory = "output"
logfile = None

def run_command(command):
    "Run command and collect errors in log file."
    (status, output) = commands.getstatusoutput(command)
    if status == 0:
        return True
    global logfile
    if logfile is None:
        logfile = open("../error.log", "w")
    logfile.write(output + "\n")
    return False

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
        ok = run_command("ffc -d %s" % f)

        # Check status
        if ok:
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
        reference_file = "../references/%s" % f
        if os.path.isfile(reference_file):
            reference_code = open(reference_file).read()
        else:
            info_blue("Missing reference for %s" % f)
            continue

        # Compare with reference
        if generated_code == reference_code:
            info_green("%s OK" % f)
        else:
            info_red("%s differs" % f)

    end()

def build_programs():
    "Build test programs for all test cases."

    begin("Building test programs")

    # Get a list of all files
    header_files = [f for f in os.listdir(".") if f.endswith(".h")]

    # Iterate over all files
    for f in header_files:

        # Generate test code
        prefix = f.split(".h")[0]
        num_forms = open(f).read().count("class %s_form_" % prefix.lower())
        test_file = open(prefix + ".cpp", "w")
        test_file.write(form_test_code(prefix, num_forms))
        test_file.close()

        # Compile test code
        command = "g++ `pkg-config --cflags ufc-1` -Wall -Werror -g -o %s %s.cpp" % (prefix, prefix)
        ok = run_command(command)

        # Check status
        if ok:
            info_green("%s OK" % prefix)
        else:
            info_red("%s failed" % prefix)

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
    generate_code()
    validate_code()

    # Build and validate programs
    build_programs()
    validate_programs()

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
