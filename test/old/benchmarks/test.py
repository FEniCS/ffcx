#!/usr/bin/env python

"Test suite for FFC: Compile all forms in the demo directory using different representations"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@gmail.com)"
__date__ = "2009-03-10"
__copyright__ = "Copyright (C) 2009-2010 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-01-21

import sys
import getopt
from glob import glob
from os import chdir, path, pardir, system
import time
from commands import getstatusoutput

def main(argv):
    "Main function, handle arguments and run tests accordingly"

    # Get command-line arguments
    try:
        opts, args = getopt.getopt(argv, "hr:", \
        ["help", "representation="])
    except getopt.GetoptError:
        usage()
        return 2

    # Run tests for both representations and form types as default
    representations = ["tensor", "quadrature"]
    form_types = ["ufl"]

    # Get options
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            return 0
        elif opt in  ("-r", "--representation"):
            if arg in representations:
                representations = [arg]
            else:
                usage()
                return 2
        else:
            usage()
            return 2

    # CD to demo directory and get all possible forms
    chdir(path.join(pardir, pardir, "demo"))
    all_form_files = []
    for t in form_types:
        all_form_files += glob("*." + t)

    # Create sorted list of unique files
    all_form_files = list(set([path.splitext(f)[0] for f in all_form_files]))
    all_form_files.sort()

    # Strip file extension from user input and check if we have valid forms
    args = [path.splitext(a)[0] for a in args]
    if not args:
        args = all_form_files
    omitted = []
    form_files = []
    for a in args:
        if a in all_form_files:
            form_files.append(a)
        else:
            omitted.append(a)

    # Print test options
    print "\nThe following test options will be used"
    print "====================================================================\n"
    print "Representations:        %s" % ", ".join(representations)
    print "\nForm types:             %s" % ", ".join(["*." + t for t in form_types])
    print "\nFor the following forms:  "
    print "\n".join(form_files)
    if omitted:
        print "\n*** Warning: The following forms are omitted because they are not present in the demo directory"
        print "\n".join(["*** " + o for o in omitted])
    print "\n====================================================================\n"

    # Get values for padding
    max_f = max([len(f) for f in form_files])
    max_t = max([len(t) for t in form_types])
    max_r = max([len(r) for r in representations])

    # Create empty summary dictionary
    # Loop all form types and representation and get summary
    summary = {}
    for form_file in form_files:
        summary[form_file] = {}
        for form_type in form_types:
            for representation in representations:
                bench = "Benchmarking: " + form_file + " "*(max_f - len(form_file))
                t = ", type: " + form_type + " "*(max_t - len(form_type))
                r = ", representation: " + representation + " "*(max_r - len(representation))
                print bench + t + r
                summary[form_file][(form_type, representation)] = compile_forms(form_type, representation, form_file)

    # Print summary
    if summary:
        print_summary(summary, form_types, representations)

    return 0

def compile_forms(form_type, representation, file_name):
    "Compile the given form, and return the compile time and file size"
    compile_time = 0
    file_size = 0

    # We want to compare the parts of the compilation where the representations
    # differ. There are probably more functions that should be switched off.
    options = "-s -fno-evaluate_basis -fno-evaluate_basis_derivatives"

    form_file = file_name + "." + form_type
    start = time.time()
    (error, output) = getstatusoutput("python -s %s %s -r %s %s" % (path.join(pardir, "scripts", "ffc"), options, representation, form_file))
    if error == 0:
        compile_time =  time.time() - start
        file_size = path.getsize(file_name + ".h")/1024.0

    return (compile_time, file_size)

def print_summary(summary, form_types, representations):
    "Print nicely formatted summary"

    # Some text strings
    time_head = " time [s] "
    size_head = " size [KB] "
    time_format = "%.5f"
    size_format = "%.3f"
    failed = "---"

    # Text info
    max_name = max([len(f) for f in summary] + [len("Form files"), len("Best performance:")]) + 4
    form_file_header = center_text(max_name, "Form files")

    # Create column headers
    column_headers = ["(%s, %s)" %(t, r) for t in form_types for r in representations]
    max_column = max([len(c) for c in column_headers] + [len(size_head)*2])
    max_column += max_column%2
    column_headers = [center_text(max_column, c) for c in column_headers]

    # Create the top header, total width of table
    top_head = " "*len(form_file_header) + "|" + "|".join(column_headers) + "|"

    # Adjust the time and size headers to half of the column header
    time_head = center_text(max_column/2, time_head)
    size_head = center_text(max_column/2, size_head)
    time_size_headers = [time_head + size_head for t in form_types for r in representations]
    time_size_head = form_file_header + "|" + "|".join(time_size_headers) + "|"

    sorted_files = [f for f in summary]
    sorted_files.sort()
    columns = {}
    performance = {}
    for t in form_types:
        for r in representations:
            performance[(t, r)] = [0, 0]
    num_forms = 0
    # Create column entries
    for form_file, val in summary.iteritems():
        num_forms += 1
        row = left_text(max_name, form_file) + "|"
        min_t = None
        min_s = None
        # Get best time and file size
        for t in form_types:
            for r in representations:
                _time, _size = val[(t, r)]
                if min_t == None and _time:
                    min_t = _time
                if min_s == None and _size:
                    min_s = _size
                if _time:
                    min_t = min(min_t, _time)
                if _size:
                    min_s = min(min_s, _size)

        # Write columns
        for t in form_types:
            for r in representations:
                _time, _size = val[(t, r)]
                if not _time and not _size:
                    # Add column to row
                    row += right_text(max_column/2, failed) + right_text(max_column/2, failed) + "|"
                    continue
                if _time > min_t and min_t:
                    _time = time_format % (_time/min_t)
                else:
                    performance[(t,r)][0] += 1
                    _time = "*" + time_format % _time

                if _size > min_s and min_s:
                    _size = size_format % (_size/min_s)
                else:
                    performance[(t,r)][1] += 1
                    _size = "*" + size_format % _size
                # Add column to row
                row += right_text(max_column/2, _time) + right_text(max_column/2, _size) + "|"
        columns[form_file] = row

    # Start printing
    print ""
    print center_text(len(top_head), "*** SUMMARY ***")
    print "="*len(top_head)
    print ""
    print top_head
    print time_size_head
    print "-"*len(top_head)
    for f in sorted_files:
        print columns[f]
    p = left_text(max_name, "Performance:") + "|"
    for t in form_types:
        for r in representations:
            time, size = performance[(t, r)]
            p += right_text(max_column/2, "%d/%d" % (time, num_forms)) + right_text(max_column/2, "%d/%d" % (size, num_forms)) + "|"
    print p
    print ""
    print "Note:  '*' denotes the best measure, all other values are ratios e.g., file_size/(*file_size)"
    print ""


def center_text(width, text):
    space = width - len(text)
    return ' '*(space/2) + text + ' '*(space/2 + space%2)

def right_text(width, text):
    space = width - len(text)
    return ' '*space + text

def left_text(width, text):
    space = width - len(text)
    return text + ' '*space

def usage():
    "Display usage info."
    print """\nUsage: ./test.py [OPTION] [FILES (optional)], or: python test.py [OPTION] [FILES (optional)]

  -h, --help              display this info

  -r, --representation    specify the representation ('tensor' or 'quadrature') for
                          which to run the tests. E.g.,

                          ./test.py -r quadrature

                           - will test ALL forms in ../../demo using quadrature representation

                          ./test.py -r tensor Poisson.ufl

                           - will test ../../demo/Poisson.ufl using tensor representation

                          If no representation is specified the tests will be
                          run for both representations.
"""
    return


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))



