#!/usr/bin/env python
import sys
import instant

s1, o1 = instant.get_status_output("python runpytests.py")
if s1 != 0:
    print "Python tests FAILED, see python_tests.log"
else:
    print "Python tests PASSED"
f = open('python_tests.log', 'w')
f.write(o1)
f.close()

s2, o2 = instant.get_status_output("make")
if s2 != 0:
    print "C++ tests FAILED, see cpp_tests.log"
else:
    print "C++ tests PASSED"
f = open('cpp_tests.log', 'w')
f.write(o2)
f.close()

sys.exit(s1 or s2)
