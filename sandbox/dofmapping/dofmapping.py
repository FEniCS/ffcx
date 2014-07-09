
from six.moves import xrange

class EntitySet(object):
    def __init__(self, domainid, cellname, entitydim):
        self.domainid = domainid
        self.cellname = cellname
        self.entitydim = entitydim

    def hashdata(self):
        return (self.domainid, self.cellname, self.entitydim)

    def __hash__(self):
        return hash(self.hashdata())

    def __eq__(self, other):
        return isinstance(other, EntitySet) and self.hashdata() == other.hashdata()

    def __mul__(self, multiplicity):
        return DofSet({self: multiplicity})

    def __rmul__(self, multiplicity):
        return multiplicity*self

    def __add__(self, other):
        return 1*self + 1*other

    def __radd__(self, other):
        return other + self

    def __str__(self):
        return "entities of dim %d on a %s in domain %s" % (\
            self.entitydim, self.cellname, self.domainid)

class DofSet(object):
    def __init__(self, multiplicities=None):
        self.multiplicities = dict(multiplicities) if multiplicities else {}

    def __add__(self, other):
        dm = DofMap()
        dm += self
        dm += other
        return dm

    def __iadd__(self, other):
        sm = self.multiplicities
        om = other.multiplicities
        for es, ok in om.iteritems():
            sm[es] = sm.get(es,0) + ok
        return self

    def __mul__(self, factor):
        ds = DofSet(self.multiplicities)
        ds *= factor
        return ds

    def __rmul__(self, factor):
        return self*factor

    def __imul__(self, factor):
        for es, k in self.multiplicities.items():
            self.multiplicities[es] = k * factor
        return self
    
    def __str__(self):
        s = []
        for es, k in sorted(self.multiplicities.items(), key=lambda x: x[0].hashdata()):
            l = "%d %s" % (k, es)
            s.append(l)
        return '\n'.join(s)


def compute_dofranges(dofset, meshsizes, offset):
    dofrange = {}
    domainids = sorted(set([es.domainid for es in dofset.multiplicities]))
    for domainid in domainids:
        dofrange[domainid] = {}

        cellnames = sorted(set([es.cellname for es in dofset.multiplicities
                                if es.domainid == domainid]))

        for cellname in cellnames:
            dofrange[domainid][cellname] = {}

            entitydims = sorted(set([es.entitydim for es in dofset.multiplicities
                                 if (es.domainid, es.cellname) == (domainid, cellname)]))

            for entitydim in entitydims:
                dofrange[domainid][cellname][entitydim] = {}

                key = EntitySet(domainid, cellname, entitydim)
                multiplicity = dofset.multiplicities.get(key, 0)

                for entitydofnumber in xrange(multiplicity):
                    start = offset
                    count = meshsizes[key]
                    end = start + count
                    offset += count
                    dofrange[domainid][cellname][entitydim][entitydofnumber] = (start, end)

    return dofrange

def test():
    print str(EntitySet("A", "triangle", 2))
    print str(5*DofSet({EntitySet("A", "triangle", 2): 2}))

    es1 = EntitySet("A", "triangle", 1)
    es2 = EntitySet("A", "triangle", 0)
    dofset = DofSet({es1: 1, es2: 2})

    print str(dofset)

    meshsizes = { es1: 7, es2: 11 }

    offset = 0

    print compute_dofranges(dofset, meshsizes, offset)

