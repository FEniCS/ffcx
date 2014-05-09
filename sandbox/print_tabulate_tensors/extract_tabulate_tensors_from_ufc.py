#!/usr/bin/env python

import re
import sys
filenames = sys.argv[1:]

#tt = re.compile("^void (.*::)*tabulate_tensor")
tt = re.compile("void .*tabulate_tensor")

#bb = re.compile("^[ ]*{")
#be = re.compile("^[ ]*}")
bb = re.compile("{")
be = re.compile("}")

db = re.compile("#ifdef DEBUG")
de = re.compile("#endif // DEBUG")

for fn in filenames:
    print "========= tabulate_tensor in '%s':" % fn
    lines = open(fn).readlines()
    in_tabulate_tensor = False
    skip_lines = 0
    in_debug_code = False
    braces = 0
    for line in lines:

        if bb.search(line):
            braces += 1
        if be.search(line):
            braces -= 1

        if db.search(line):
            in_debug_code = True
        elif de.search(line):
            in_debug_code = False

        if in_tabulate_tensor:
            if skip_lines:
                skip_lines -= 1
                continue

            if not in_debug_code:
                print line,

            if braces == outside_tabulate_tensor_braces:
                in_tabulate_tensor = False
                break
        else:
            if tt.search(line):
                in_tabulate_tensor = True
                outside_tabulate_tensor_braces = braces
                skip_lines = 3

