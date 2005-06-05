__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2005-05-16 -- 2005-05-17"
__copyright__ = "Copyright (c) 2005 Anders Logg"
__license__  = "GNU GPL Version 2"

# FIAT modules
from FIAT.dualbasis import *
from FIAT.shapes_new import *

# FFC modules
from declaration import *

# FIXME: Should not be DOLFIN-specific
format = { "point"     : lambda x : "points[%s]" % x,
           "affinemap" : lambda x : "map(%s)" % x,
           "component" : lambda i : "components[%s]" % i }

class PointMap:

    """A PointMap maps the coordinates of the degrees of freedom on
    the reference element to physical coordinates."""

    def __init__(self, dualbasis):
        "Create PointMap."

        # Get points (temporary until we can handle other types of nodes)
        points = dualbasis.pts

        # Map point to reference element (different in FFC and FIAT)
        newpoints = []
        for point in points:
            newpoints += [tuple([0.5*(x + 1.0) for x in point])]
        points = newpoints

        # Get the number of vector components
        num_components = dualbasis.num_reps

        self.declarations = []

        # Iterate over the dofs
        dof = 0
        for component in range(num_components):
            for point in points:

                # Coordinate
                x = (", ".join(["%.15e" % x for x in point]))
                name = format["point"](dof)
                value = format["affinemap"](x)
                self.declarations += [Declaration(name, value)]

                dof += 1

        # Iterate over the dofs
        dof = 0
        for component in range(num_components):
            for point in points:

                # Component
                name = format["component"](dof)
                value = "%d" % component
                self.declarations += [Declaration(name, value)]

                dof += 1

        #for declaration in self.declarations:
        #    print declaration.name + " = " + declaration.value
