__author__ = "Anders Logg, Kristian B. Oelgaard and Marie E. Rognes"
__date__ = "2010-01-21"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

import os, sys

# Parameters
output_directory = "output"

def generate_tests():
    "Generate form files for all test cases."

    # Create output directory if it does not exist
    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    # Copy demos files
    os.system("cp ../../demo/*.ufl %s" % output_directory)

def main(args):
    "Run all regression tests."

    # Generate form files for all test cases
    generate_tests()

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

