
from six.moves import iteritems, iterkeys
import shutil
import re
import os
from collections import defaultdict

modulename = r"uflacs"
#skipmods = ['utils', 'commands']
#skipmods = ['utils', 'commands', 'geometry', 'backends']
#skipmods = ['utils', 'commands', 'geometry', 'codeutils', 'backends']
skipmods = ['utils', 'commands']
#skipmods = []

r = re.compile(r'from +' + modulename + r'[.]([^ ]*) +import')

files = []
os.path.walk(modulename, lambda x, y, z: files.extend(os.path.join(y,fn) for fn in z if fn.endswith(".py")), None)

imports = defaultdict(dict)
pimports = {}
#for f in files:
for f in files[:]:
    mod = f.replace('uflacs/','').replace('.py', '').replace('/', '.').strip()
    pmod = mod.split('.')[0]
    if pmod in skipmods:
        continue

    lines = [l for l in open(f).readlines() if "import" in l]

    matches = [r.search(l) for l in lines]
    groups = [m.groups() for m in matches if m]

    modules = sorted({g[0] for g in groups if len(g) == 1})
    parent_modules = sorted({m.split('.')[0] for m in modules})

    modules = [m for m in modules if not any(m.startswith(skipmod+'.') for skipmod in skipmods)]
    parent_modules = [m for m in parent_modules if m not in skipmods + [pmod]]

    crap = [g for g in groups if len(g) != 1]
    if crap:
        print "Not sure what to do with this:"
        print crap

    imports[pmod][mod] = modules
    pimports[pmod] = sorted(set(parent_modules) | set(pimports.get(pmod,())))


# TODO: Make topological sorting of pmods for neater prints
pmods = sorted(iterkeys(pimports))

print
print '*'*80
print
for m in pmods:
    print "="*70, m
    l = max(len(k) for k in imports[m])
    fmt = '%s depends on'
    for k,v in iteritems(imports[m]):
        if v:
            print fmt % k
            print '\n'.join('    %s' % vv for vv in v)

print
print '*'*80
print
l = max(len(m) for m in pmods)
fmt = '%%%ds  depends on  %%s' % l
for m in pmods:
    v = pimports[m]
    if v:
        print fmt % (m, v)


def rename(x):
    return x.replace('.','__')

# Write submodule imports to graph
names = set()
lines = []
for m in pmods:
    lines.append("    subgraph {\n")
    for k,v in iteritems(imports[m]):
        if v:
            names.add(k)
        for t in v:
            names.add(t)
            lines.append("        %s -> %s;\n" % (rename(k), rename(t)))
    lines.append("    }\n")

for n in names:
    lines.append('    %s [label="%s"]\n' % (rename(n), n))

with open("deps.dot", "w") as df:
    df.write("digraph {\n")
    df.writelines(lines)
    df.write("}\n")


# Write coarse module dependencies to graph
names = set()
lines = []
for m in pmods:
    v = pimports[m]
    if v:
        names.add(m)
    for t in v:
        names.add(t)
        lines.append("        %s -> %s;\n" % (rename(m), rename(t)))

for n in names:
    lines.append('    %s [label="%s"]\n' % (rename(n), n))

with open("pdeps.dot", "w") as df:
    df.write("digraph {\n")
    df.writelines(lines)
    df.write("}\n")
