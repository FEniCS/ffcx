"""Build symmetry prototypes for the known symmetric benchmark forms."""

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
    updates = list(re.finditer(r"A\[(.*?)\] \+= (.*?);", source))
    n = N if case.startswith("poisson-") else N // m["dimension"]
    d = N // n
    flat = re.search(rf"for \(int i = 0; i < {N * N}; \+\+i\)\s*\{{", source)
    if flat:
        assert len(updates) == 1
        u = updates[0]
        body = f"for(int row=0;row<{N};++row)for(int col=row;col<{N};++col){{int i=row*{N}+col;double value={u[2]};A[i]+=value;if(row!=col)A[col*{N}+row]+=value;}}"
        brace = source.index("{", flat.start())
        depth = 1
        end = brace + 1
        while depth:
            depth += (source[end] == "{") - (source[end] == "}")
            end += 1
        transformed = source[: flat.start()] + body + source[end:]
    elif "for (" not in source[source.index("void tabulate_tensor") :]:

        def replace(u):
            """Rewrite one matched output update."""
            index = int(eval(u[1], {"__builtins__": {}}, {}))
            row, col = divmod(index, N)
            if row > col:
                return ""
            if row == col:
                return u[0]
            return f"{{double value={u[2]};A[{index}]+=value;A[{col * N + row}]+=value;}}"

        transformed = re.sub(r"A\[(.*?)\] \+= (.*?);", replace, source)
    else:
        assert len(updates) == d * d, (case, len(updates), d)
        counter = [0]

        def replace(u):
            """Rewrite one matched output update."""
            a, b = divmod(counter[0], d)
            counter[0] += 1
            mirror = f"{N}*({d}*j+{b})+{d}*i+{a}"
            return f"{{if(i<j || {a}<={b}){{double value={u[2]};A[{u[1]}]+=value;if(i!=j || {a}!={b})A[{mirror}]+=value;}}}}"

        transformed = re.sub(r"A\[(.*?)\] \+= (.*?);", replace, source)
        transformed = transformed.replace(
            f"for (int j = 0; j < {n}; ++j)", f"for (int j = i; j < {n}; ++j)"
        )
    p = r / f"{case}-symmetric.c"
    p.write_text(transformed)
    m["modes"].append("symmetric")
    for cc in ["gcc-16", "clang-22"]:
        dest = r / cc
        obj = dest / f"{case}-symmetric.o"
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
                str(dest / f"{case}-symmetric.so"),
            ],
            check=True,
        )
    print(case, "symmetric", flush=True)
(r / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
