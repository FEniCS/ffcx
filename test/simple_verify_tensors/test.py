from cppcode import tabulate_tensor_code

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2009-03-15 -- 2009-03-15"
__copyright__ = "Copyright (C) 2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

from ufl.common import tstr
import sys, os, commands, pickle, numpy

#format = "form"
format = "ufl"

def tabulate_tensor(integral, header):
    "Generate code and tabulate tensor for integral."

    print "  Checking %s..." % integral

    # Generate and compile code
    options = {"n": 100, "N": 1000, "integral": integral, "header": header}
    code = tabulate_tensor_code % options
    open("tabulate_tensor.cpp", "w").write(code)
    (status, dummy) = commands.getstatusoutput("g++ -o tabulate_tensor tabulate_tensor.cpp")
    if not status is 0: return "GCC compilation failed"

    # Run code and get results
    (status, output) = commands.getstatusoutput("./tabulate_tensor")
    if not status is 0: return "Unable to tabulate tensor (segmentation fault?)"
    values = [float(value) for value in output.split(" ") if len(value) > 0]

    return numpy.array(values)

def get_integrals(form_file):
    "Extract all integral classes for form file."
    prefix = form_file.split("/")[-1].split(".")[0]
    form = prefix + "." + format
    header = prefix + ".h"
    integrals = commands.getoutput("grep integral %s | grep class" % header).split("\n")
    integrals = [integral.split(":")[0].split(" ")[-1] for integral in integrals if len(integral) > 0]
    return (integrals, form, header)

def to_dict(tuples):
    "Convert list of tuples to dictionary (dictionaries can't be pickled)."
    d = {}
    for key, value in tuples:
        d[key] = value
    return d

def main(args):
    "Call tabulate tensor for all integrals found in demo directory."

    # Change to temporary folder
    if not os.path.isdir("tmp"):
        os.mkdir("tmp")
    os.chdir("tmp")

    # Iterate over all form files
    form_files = commands.getoutput("ls ../../../demo/*.%s" % format).split("\n")
    values = []
    for form_file in form_files:

        # Compile form
        (integrals, form, header) = get_integrals(form_file)
        print "Compiling form %s..." % form
        (status, dummy) = commands.getstatusoutput("ffc %s" % form_file)

        # Tabulate tensors for all integrals
        print "  Found %d integrals" % len(integrals)
        for integral in integrals:
            if status == 0:
                values.append((integral, tabulate_tensor(integral, header)))
            else:
                values.append((integral, "FFC compilation failed"))

    # Load or update reference values
    if os.path.isfile("../reference.pickle"):
        reference = to_dict(pickle.load(open("../reference.pickle", "r")))
    else:
        print "Unable to find reference values, storing current values."
        pickle.dump(values, open("../reference.pickle", "w"))        
        return 0
        
    # Compare results
    print ""
    tol = 1e-12
    results = []
    for (integral, value) in values:
        if not isinstance(value, str):
            if integral in reference:
                e = max(abs(value - reference[integral]))
                if e < tol:
                    result = "OK  (diff = %g)" % e
                else:
                    result = "*** (diff = %g)" % e
            else:
                result = "missing reference"
        else:
            result = value
        results.append((integral, result))
    print tstr(results)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
