__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2005-03-24"
__copyright__ = "Copyright (c) 2005 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
import sys

# FIXME: Get from width of terminal window
width = 80

class Progress:
    "A simple text-mode progress bar."

    def __init__(self, n):
        "Create progress bar for process consisting of n steps."
        self.n = n
        self.i = 0
        self.pos = 0
        #print "".join(["_" for j in range(width)])
        return

    def __iadd__(self, other):
        "Add increment to progress bar."
        self.i += other      
        newpos = int(float(self.i) / float(self.n) * float(width))
        if newpos > width:
            return
        if newpos > self.pos:
            sys.stdout.write("".join(["." for j in range(newpos - self.pos)]))
            sys.stdout.flush()
            self.pos = newpos
        if newpos == width:
            sys.stdout.write("\n")
            sys.stdout.flush();
        return self
