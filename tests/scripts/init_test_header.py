#!/usr/bin/env python

import sys
sys.path.insert(0, "py")
from codegentestcase import create_initial_test_header

if __name__ == "__main__":
    arg, = sys.argv[1:]
    create_initial_test_header(arg)

