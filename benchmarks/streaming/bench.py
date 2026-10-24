"""Validate and time experimental kernels over DRAM-resident cell records."""

# Generated C fragments retain their original source layout.

import ctypes
import json
import os
import re
from pathlib import Path

import numpy as np

r = Path(os.environ.get("FFCX_STREAM_WORKDIR", "/tmp/ffcx-streaming-investigation"))
manifest = json.loads((r / "manifest.json").read_text())
os.sched_setaffinity(0, {2})
ptr = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
report = {}
output = r / os.environ.get("STREAM_OUTPUT", "screen.json")
for case, m in manifest["cases"].items():
    if os.environ.get("STREAM_FILTER") and not re.search(os.environ["STREAM_FILTER"], case):
        continue
    if os.environ.get("STREAM_MODES"):
        m["modes"] = [mode for mode in os.environ["STREAM_MODES"].split(",") if mode in m["modes"]]
    if len(m["modes"]) <= 1:
        continue
    data = np.load(r / f"{case}.npz")
    xs = np.ascontiguousarray(data["coordinates"])
    ws = np.ascontiguousarray(data["coefficients"])
    ref = np.ascontiguousarray(data["expected"])
    mods = {}
    samples = {}
    errors = {}
    for cc in ["gcc-16", "clang-22"]:
        baseline = None
        for mode in m["modes"]:
            label = cc + "-" + mode
            mod = ctypes.CDLL(str(r / cc / f"{case}-{mode}.so"))
            mod.evaluate.argtypes = [ptr, ptr, ptr]
            mod.prepare.argtypes = [ptr, ptr, ctypes.c_long]
            mod.measure.argtypes = [ctypes.c_long]
            mod.measure.restype = ctypes.c_double
            mod.verify.argtypes = [ptr]
            mod.verify.restype = ctypes.c_double
            mod.cell_count.restype = ctypes.c_long

            def evaluate(coords, nonzero=False):
                """Evaluate all validation cells with optional nonzero output."""
                values = []
                for k, x in enumerate(coords):
                    init = (
                        np.sin(np.arange(m["size"]) + k) * 0.23 if nonzero else np.zeros(m["size"])
                    )
                    a = np.ascontiguousarray(init.copy())
                    mod.evaluate(a, x, ws[k])
                    assert np.isfinite(a).all()
                    values.append(a - init)
                return np.array(values)

            actual = evaluate(xs)
            error = 0
            if mode == "baseline":
                baseline = mod
                if case.startswith("hyper-"):
                    ref = np.ascontiguousarray(actual)
                negxs = xs.copy().reshape(32, -1, 3)
                negxs[:, :, 0] *= -1
                negxs = negxs.reshape(xs.shape)
                negref = evaluate(negxs)
            np.testing.assert_allclose(actual, ref, rtol=2e-10, atol=2e-11)
            error = max(error, float(np.max(abs(actual - ref))))
            nonzero = evaluate(xs, True)
            np.testing.assert_allclose(nonzero, ref, rtol=2e-10, atol=2e-11)
            error = max(error, float(np.max(abs(nonzero - ref))))
            negative = evaluate(negxs)
            np.testing.assert_allclose(negative, negref, rtol=2e-10, atol=2e-11)
            error = max(error, float(np.max(abs(negative - negref))))
            mod.prepare(xs, ws, int(os.environ.get("STREAM_MIB", "256")) * 1024**2)
            mod.measure(1)
            mods[label] = mod
            samples[label] = []
            errors[label] = error
    labels = list(mods)
    for rep in range(int(os.environ.get("STREAM_REPS", "3"))):
        order = labels[rep % len(labels) :] + labels[: rep % len(labels)]
        for label in order:
            samples[label].append(mods[label].measure(1))
    rows = {}
    for label, mod in mods.items():
        # Each compiler's reference is numerically equivalent to the independent tensor.
        err = mod.verify(ref)
        assert err < 2e-10, (case, label, err)
        cells = mod.cell_count()
        ss = samples[label]
        rows[label] = {
            "samples_s": ss,
            "best_ns": min(ss) / cells * 1e9,
            "median_ns": float(np.median(ss)) / cells * 1e9,
            "useful_gbs": cells * m["useful_bytes_per_cell"] / min(ss) / 1e9,
            "max_error": max(err, errors[label]),
        }
        mod.release()
    report[case] = rows
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        case,
        {
            cc: {
                mode: round(rows[cc + "-baseline"]["best_ns"] / rows[cc + "-" + mode]["best_ns"], 3)
                for mode in m["modes"]
                if mode != "baseline"
            }
            for cc in ["gcc-16", "clang-22"]
        },
        flush=True,
    )
