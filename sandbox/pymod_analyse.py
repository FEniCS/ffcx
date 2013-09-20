import shutil
import re
import os
from collections import defaultdict

modulename = r"uflacs"
skipmods = []#['utils']

r = re.compile(r'from +' + modulename + r'[.](.*) +import')

files = []
os.path.walk(modulename, lambda x, y, z: files.extend(os.path.join(y,fn) for fn in z if fn.endswith(".py")), None)

imports = defaultdict(dict)
pimports = {}
#for f in files:
for f in files[:]:
    mod = f.replace('uflacs/','').replace('.py', '').replace('/', '.')
    pmod = mod.split('.')[0]

    lines = [l for l in open(f).readlines() if "import" in l]

    matches = [r.search(l) for l in lines]
    groups = [m.groups() for m in matches if m]

    modules = sorted(set([g[0] for g in groups if len(g) == 1]))
    parent_modules = sorted(set([m.split('.')[0] for m in modules]))

    modules = [m for m in modules if not any(m.startswith(skipmod+'.') for skipmod in skipmods)]
    parent_modules = [m for m in parent_modules if m not in skipmods + [pmod]]

    crap = [g for g in groups if len(g) != 1]   
    if crap:
        print "Not sure what to do with this:"
        print crap

    imports[pmod][mod] = modules
    pimports[pmod] = sorted(set(parent_modules) | set(pimports.get(pmod,())))

# TODO: Make topological sorting of pmods
pmods = sorted(pimports.keys())

print '*'*80
for m in pmods:
    for k,v in imports[m].iteritems():
        if v:
            print '%s: %s' % (k, v)

print '*'*80
for m in pmods:
    v = pimports[m]
    if v:
        print '%s: %s' % (m, v)

