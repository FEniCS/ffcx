#!/usr/bin/env python
import sys
import instant

s1, o1 = instant.get_status_output("python runpytests.py")
if s1 != 0:
    print "Python tests failed, see python_tests.log"

s2, o2 = instant.get_status_output("make")
if s2 != 0:
    print "C++ tests failed, see cpp_tests.log"

f = open('python_tests.log', 'w')
f.write(o1)
f.close()

f = open('cpp_tests.log', 'w')
f.write(o2)
f.close()

sys.exit(s1 or s2)
