#!/usr/bin/env python

if __name__ == "__main__":
    import sys
    from codegentestcase import create_initial_test_header
    arg, = sys.argv[1:]
    create_initial_test_header(arg)
