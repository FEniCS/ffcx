__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2009-03-15 -- 2009-03-15"
__copyright__ = "Copyright (C) 2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

from cppcode import tabulate_tensor_code, integral_types
from ufl.common import tstr
import sys, os, commands, pickle, numpy

# Temporary while testing new UFL compiler
#format = "form"
format = "ufl"

# Log file
logfile = None

# Set paths
if not "PATH" in os.environ: os.environ["PATH"] = ""
if not "PYTHONPATH" in os.environ: os.environ["PYTHONPATH"] = ""
os.environ["PATH"] = "../../../scripts:" + os.environ["PATH"]
os.environ["PYTHONPATH"] ="../../..:" + os.environ["PYTHONPATH"]

def tabulate_tensor(integral, integral_type, header):
    "Generate code and tabulate tensor for integral."

    print "  Checking %s..." % integral

    # Generate and compile code
    options = {"n": 100, "N": 1000, "integral": integral, "header": header}
    code = tabulate_tensor_code[integral_type] % options
    open("tabulate_tensor.cpp", "w").write(code)
    (ok, output) = run_command("g++ -o tabulate_tensor tabulate_tensor.cpp")
    if not ok: return "GCC compilation failed"

    # Run code and get results
    (ok, output) = run_command("./tabulate_tensor")
    if not ok: return "Unable to tabulate tensor (segmentation fault?)"
    values = [float(value) for value in output.split(" ") if len(value) > 0]

    return numpy.array(values)

def get_integrals(form_file):
    "Extract all integral classes for form file."
    prefix = form_file.split("/")[-1].split(".")[0]
    form = prefix + "." + format
    header = prefix + ".h"
    integrals = []
    for integral_type in integral_types:
        for integral in commands.getoutput("grep ufc::%s %s | grep class" % (integral_type, header)).split("\n"):
            integral = integral.split(":")[0].split(" ")[-1].strip()
            if not integral is "":
                integrals.append((integral, integral_type))
    return (integrals, form, header)

def to_dict(tuples):
    "Convert list of tuples to dictionary (dictionaries can't be pickled)."
    d = {}
    for key, value in tuples:
        d[key] = value
    return d

def run_command(command):
    "Run system command and collect any errors."
    global logfile
    (status, output) = commands.getstatusoutput(command)
    if not status is 0:
        if logfile is None:
            logfile = open("../error.log", "w")
        logfile.write(output + "\n")
    return (status == 0, output)

def check_results(values, reference):
    "Check results and print summary."

    print ""
    tol = 1e-12
    results = []
    ok = True
    for (integral, value) in values:
        if not isinstance(value, str):
            if integral in reference:
                e = max(abs(value - reference[integral]))
                if e < tol:
                    result = "OK" % e
                else:
                    result = "*** (diff = %g)" % e
                    ok = False
            else:
                result = "missing reference"
        else:
            result = value
            ok = False
        results.append((integral, result))
    print tstr(results, 100)

    if ok:
        print "\nAll tensors verified OK"
        return 0
    else:
        print "\nVerification failed."
        return 1
    
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
        (ok, output) = run_command("ffc %s %s" % (" ".join(args), form_file))

        # Tabulate tensors for all integrals
        print "  Found %d integrals" % len(integrals)
        for (integral, integral_type) in integrals:
            if ok:
                values.append((integral, tabulate_tensor(integral, integral_type, header)))
            else:
                values.append((integral, "FFC compilation failed"))

    # Load or update reference values
    if os.path.isfile("../reference.pickle"):
        reference = to_dict(pickle.load(open("../reference.pickle", "r")))
    else:
        print "Unable to find reference values, storing current values."
        pickle.dump(values, open("../reference.pickle", "w"))        
        return 0

    # Check results
    return check_results(values, reference)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
