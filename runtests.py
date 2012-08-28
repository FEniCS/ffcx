#!/usr/bin/env python
import io, os, sys, unittest, glob

verbosity = 2
modulename = "uflacs"
testmodulename = "test_uflacs"
modulepath = "site-packages"

if __name__ == "__main__":
    # Import local uflacs and tests
    sys.path.insert(0, modulepath)
    mod = __import__(modulename)
    testmod = __import__(testmodulename)
    print "Running tests with %s version %s, last changed %s." % (
        modulename, mod.__version__, mod.__date__)

    # Running tests from all test_foo.py files
    pattern = os.path.join(modulepath, testmodulename, "test*.py")
    tests = sorted(os.path.basename(f).replace(".py", "") for f in glob.glob(pattern))

    # Setup unittest 
    loader = unittest.TestLoader()
    fullsuite = unittest.TestSuite()
    for case in tests:
        submodulename = '.'.join((testmodulename, case))
        casemodule = __import__(submodulename, fromlist=[testmodulename])
        print casemodule
        casesuite = loader.loadTestsFromModule(casemodule)
        fullsuite.addTests(casesuite)
    unittest.TextTestRunner(verbosity=verbosity).run(fullsuite)

