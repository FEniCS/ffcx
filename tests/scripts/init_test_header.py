#!/usr/bin/env python

import sys, glob
sys.path.insert(0, "py")
from codegentestcase import create_initial_test_header#, update_testall_includes

def update_testall_includes():
    f = open('cpp/main_testall.h', 'w')
    for fn in glob.glob('cpp/test_*.h'):
        line = '#include "%s"\n' % fn.replace('cpp/','')
        print line
        f.write(line)
    f.close()

if __name__ == "__main__":
    args = sys.argv[1:]
    if args:
        arg, = args
        create_initial_test_header(arg)
    update_testall_includes()

