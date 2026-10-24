"""Build blocked symmetry prototypes for the known symmetric benchmark forms."""

# Generated C fragments retain their original source layout.
# ruff: noqa: E501

import json
import os
import re
import subprocess
from pathlib import Path

r = Path(os.environ.get("FFCX_STREAM_WORKDIR", "/tmp/ffcx-streaming-investigation"))
manifest = json.loads((r / "manifest.json").read_text())
flags = manifest["flags"]
for case, m in manifest["cases"].items():
    if m["size"] != m["dofs"] ** 2:
        continue
    source = (r / f"{case}-baseline.c").read_text()
    N = m["dofs"]
    n = N if case.startswith("poisson-") else N // m["dimension"]
    d = N // n
    if n < 6:
        continue
    updates = list(re.finditer(r"A\[(.*?)\] \+= (.*?);", source))
    flat = re.search(rf"for \(int i = 0; i < {N * N}; \+\+i\)\s*\{{", source)
    if flat:
        assert len(updates) == 1
        start = flat.start()
        brace = source.index("{", start)
        bodyupdates = [updates[0][2].replace("[i]", "[i*" + str(N) + "+j]")]
    else:
        assert len(updates) == d * d, (case, len(updates), d)
        loops = list(re.finditer(rf"for \(int i = 0; i < {n}; \+\+i\)\s*\{{", source))
        start = loops[-1].start()
        brace = source.index("{", start)
        bodyupdates = [u[2] for u in updates]
    end = brace + 1
    depth = 1
    while depth:
        depth += (source[end] == "{") - (source[end] == "}")
        end += 1
    for B in [2, 4, 8]:
        mode = f"blocksym{B}"
        D = B * d
        computes = "\n".join(
            f"tile[{d}*(i-ib)+{a}][{d}*(j-jb)+{b}]={bodyupdates[a * d + b]};"
            for a in range(d)
            for b in range(d)
        )
        body = f"""for(int ib=0;ib<{n};ib+={B}){{int ie=ib+{B}<{n}?ib+{B}:{n};
 for(int jb=ib;jb<{n};jb+={B}){{int je=jb+{B}<{n}?jb+{B}:{n};
 double tile[{D}][{D}];
 for(int i=ib;i<ie;++i)for(int j=jb;j<je;++j){{{computes}}}
 for(int ii=0;ii<{d}*(ie-ib);++ii)for(int jj=0;jj<{d}*(je-jb);++jj)A[({d}*ib+ii)*{N}+{d}*jb+jj]+=tile[ii][jj];
 if(ib!=jb)for(int jj=0;jj<{d}*(je-jb);++jj)for(int ii=0;ii<{d}*(ie-ib);++ii)A[({d}*jb+jj)*{N}+{d}*ib+ii]+=tile[ii][jj];
 }}
 }}"""
        transformed = source[:start] + body + source[end:]
        p = r / f"{case}-{mode}.c"
        p.write_text(transformed)
        m["modes"].append(mode)
        for cc in ["gcc-16", "clang-22"]:
            dest = r / cc
            obj = dest / f"{case}-{mode}.o"
            root = r / "baseline" / case.split("-")[0]
            subprocess.run([cc, *flags, f"-I{root}", "-c", str(p), "-o", str(obj)], check=True)
            subprocess.run(
                [
                    cc,
                    "-shared",
                    str(obj),
                    str(dest / f"{case}-driver.o"),
                    str(dest / f"{case}-verify.o"),
                    "-lm",
                    "-o",
                    str(dest / f"{case}-{mode}.so"),
                ],
                check=True,
            )
    print(case, "blocked symmetry", flush=True)
(r / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
