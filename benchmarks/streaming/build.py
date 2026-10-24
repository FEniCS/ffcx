"""Build experimental streamed-kernel variants and validation drivers."""

# Generated C fragments retain their original source layout.
# ruff: noqa: E501

import json
import os
import re
import subprocess
from pathlib import Path

import numpy as np

from benchmarks.elasticity import inputs

r = Path(os.environ.get("FFCX_STREAM_WORKDIR", "/tmp/ffcx-streaming-investigation"))
flags = ["-O3", "-march=native", "-ffp-contract=fast", "-ffast-math", "-fPIC"]
cases = {}


def balanced(s):
    """Balance sums in generated output contractions."""

    def replace(m):
        """Rewrite one matched output update."""
        expr = m[2]
        parts = []
        start = 0
        depth = 0
        for j, c in enumerate(expr):
            depth += (c in "([") - (c in ")]")
            if not depth and expr[j : j + 3] == " + ":
                parts.append(expr[start:j])
                start = j + 3
        parts.append(expr[start:])

        def tree(a):
            """Apply one step of the experimental transformation."""
            if len(a) == 1:
                return a[0]
            k = len(a) // 2
            return "(" + tree(a[:k]) + " + " + tree(a[k:]) + ")"

        return m[1] + tree(parts) + ";"

    return re.sub(r"(A\[[^\n]+?\] \+= )([^;\n]+);", replace, s)


for suite in ["poisson", "elasticity", "hyper"]:
    root = r / "baseline" / suite
    change = r / "scaling" / suite
    if suite == "hyper":
        entries = {
            f"{kind}{p}": {
                "degree": p,
                "dimension": 3,
                "dofs": len(np.load(root / f"nodes{p}.npy")) * 3,
                "kind": kind,
            }
            for p in [1, 2, 3]
            for kind in ["residual", "tangent"]
        }
    else:
        entries = json.loads((root / "manifest-ffcx.json").read_text())["cases"]
    for case, m in entries.items():
        name = suite + "-" + case
        if suite == "hyper":
            ws, xs = inputs(np.load(root / f"nodes{m['degree']}.npy"))
            size = m["dofs"] if m["kind"] == "residual" else m["dofs"] ** 2
            source = root / f"{case}.c"
            new = change / f"{case}.c"
            expected = np.zeros((32, size))
        else:
            data = np.load(root / f"{case}.npz")
            xs = data["ffcx_coordinates"]
            ws = np.zeros((32, 0))
            size = m["dofs"] ** 2
            source = root / f"ffcx_{case}.c"
            new = change / f"ffcx_{case}.c"
            expected = data["expected"].reshape(32, -1)
        np.savez(r / f"{name}.npz", coordinates=xs, coefficients=ws, expected=expected)
        text = source.read_text()
        fn = re.search(r"void (tabulate_tensor_\w+)\(", text)[1]
        nw = ws.shape[1]
        nx = xs.shape[1]
        wo = max(128, (nx * 8 + 63) // 64 * 64)
        ao = wo + (nw * 8 + 63) // 64 * 64
        stride = ao + (size * 8 + 63) // 64 * 64
        variants = {
            "baseline": text,
            "scaling": new.read_text(),
            "aligned": re.sub(
                r"static const (double|float) ", r"alignas(64) static const \1 ", text
            ),
            "balanced": balanced(text),
            "combined": balanced(
                re.sub(r"static const (double|float) ", r"alignas(64) static const \1 ", text)
            ),
            "fullscaling": (r / "fullscaling" / suite / source.name).read_text(),
            "aliases": (r / "aliases" / suite / source.name).read_text(),
        }
        start = text.index("{", text.index("void " + fn))
        close = start + 1
        depth = 1
        while depth:
            depth += (text[close] == "{") - (text[close] == "}")
            close += 1
        close -= 1
        if not re.search(r"for \(", text[start:close]) and size <= 144:
            variants["staged"] = (
                text[: start + 1]
                + f"\ndouble local_tensor[{size}] = {{0}};\n"
                + re.sub(r"\bA\[", "local_tensor[", text[start + 1 : close])
                + f"\nfor(int z=0;z<{size};++z)A[z]+=local_tensor[z];\n"
                + text[close:]
            )
        cases[name] = {
            **m,
            "size": size,
            "nx": nx,
            "nw": nw,
            "stride": stride,
            "useful_bytes_per_cell": nx // 3 * m["dimension"] * 8 + nw * 8 + size * 8,
            "modes": list(variants),
        }
        driver = f"""#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
void {fn}(double*,const double*,const double*,const double*,const int*,const uint8_t*,void*);
static unsigned char*records;static long cells;double checksum;
void prepare(const double*xs,const double*ws,long bytes){{if(records)return;cells=bytes/{stride};if(posix_memalign((void**)&records,64,cells*{stride}))abort();memset(records,0,cells*{stride});for(long k=0;k<cells;++k){{memcpy(records+k*{stride},xs+(k%32)*{nx},{nx}*8);memcpy(records+k*{stride}+{wo},ws+(k%32)*{nw},{nw}*8);}}}}
void release(void){{free(records);records=0;}}
void evaluate(double*A,const double*x,const double*w){{{fn}(A,w,NULL,x,NULL,NULL,NULL);}}
long cell_count(void){{return cells;}}
unsigned char*record_buffer(void){{return records;}}
double measure(long passes){{double check=0;struct timespec t0,t1;clock_gettime(CLOCK_MONOTONIC,&t0);for(long p=0;p<passes;++p)for(long k=0;k<cells;++k){{const double*x=(double*)(records+k*{stride});const double*w=(double*)(records+k*{stride}+{wo});double*A=(double*)(records+k*{stride}+{ao});memset(A,0,{size}*8);{fn}(A,w,NULL,x,NULL,NULL,NULL);check+=A[k%{size}];}}clock_gettime(CLOCK_MONOTONIC,&t1);checksum=check;return t1.tv_sec-t0.tv_sec+1e-9*(t1.tv_nsec-t0.tv_nsec);}}
"""
        verify = f"""#include <math.h>
unsigned char*record_buffer(void);long cell_count(void);
double verify(const double*ref){{unsigned char*r=record_buffer();double error=0;for(long k=0;k<cell_count();++k){{const double*A=(double*)(r+k*{stride}+{ao});for(long j=0;j<{size};++j){{if(!isfinite(A[j]))return INFINITY;double e=fabs(A[j]-ref[(k%32)*{size}+j]);if(e>error)error=e;}}}}return error;}}
"""
        (r / f"{name}-driver.c").write_text(driver)
        (r / f"{name}-verify.c").write_text(verify)
        for mode, t in variants.items():
            (r / f"{name}-{mode}.c").write_text(t)
        for cc in ["gcc-16", "clang-22"]:
            dest = r / cc
            dest.mkdir(exist_ok=True)
            d = dest / f"{name}-driver.o"
            v = dest / f"{name}-verify.o"
            if not d.exists():
                subprocess.run(
                    [cc, *flags, "-c", str(r / f"{name}-driver.c"), "-o", str(d)], check=True
                )
            if not v.exists():
                subprocess.run(
                    [
                        cc,
                        "-O3",
                        "-march=native",
                        "-fPIC",
                        "-c",
                        str(r / f"{name}-verify.c"),
                        "-o",
                        str(v),
                    ],
                    check=True,
                )
            for mode in variants:
                obj = dest / f"{name}-{mode}.o"
                if obj.exists():
                    continue
                subprocess.run(
                    [cc, *flags, f"-I{root}", "-c", str(r / f"{name}-{mode}.c"), "-o", str(obj)],
                    check=True,
                )
                subprocess.run(
                    [
                        cc,
                        "-shared",
                        str(obj),
                        str(d),
                        str(v),
                        "-lm",
                        "-o",
                        str(dest / f"{name}-{mode}.so"),
                    ],
                    check=True,
                )
        print(name, flush=True)
(r / "manifest.json").write_text(json.dumps({"flags": flags, "cases": cases}, indent=2) + "\n")
