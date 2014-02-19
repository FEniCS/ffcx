#!/usr/bin/env python
import os, sys
import instant

def run_and_log(command, name, logfile):
    try:
        os.remove(logfile)
    except:
        pass

    if name == "Python": # Hack
        os.chdir("../test/")

    s, o = instant.get_status_output(command)
    if s != 0:
        print "%s tests FAILED with return code %d, see %s" % (name, s, logfile)
    else:
        print "%s tests PASSED" % (name,)

    if name == "Python": # Hack
        os.chdir("../tests/")

    try:
        f = open(logfile, 'w')
        f.write(o)
        f.close()
    except:
        print "Error while writing logfile %s" % (logfile,)

    return s

s1 = run_and_log("py.test", "Python", "python_tests.log")
s2 = run_and_log("make", "C++", "cpp_tests.log")

sys.exit(s1 or s2)
