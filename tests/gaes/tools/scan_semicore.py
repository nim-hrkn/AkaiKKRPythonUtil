"""scan RUN/*/out_dos.log of the 2019 HEA run for valence -> gap (DOS < th, contiguous) -> semicore (DOS > th again).

Usage (in fukushima_HEA_run_exprlattice/production_run/run0): python scan_semicore.py semicore_scan.json
See docs/core_levels_vs_dos.md and docs/ewidth_tuning_scheme.md.
"""
"""scan RUN/*/out_dos.log: valence -> gap (DOS<th contiguous) -> semicore (DOS>th again) going down in energy."""
import os, sys, json
import numpy as np
from multiprocessing import Pool
from pymatgen.core.periodic_table import Element

RUN = "RUN"
TH = [1e-3, 2e-2]
MIN_GAP = 0.10     # Ry, minimum contiguous low-DOS width
MIN_SEMI = 0.03    # Ry, minimum width of the DOS>th region below the gap

def read_total_dos(path):
    e, up, dn = [], [], []
    with open(path) as f:
        block = 0; on = False
        for line in f:
            s = line.strip()
            if s.startswith("total DOS"):
                block += 1; on = True; continue
            if on:
                if not s:
                    on = False; continue
                p = s.split()
                try:
                    x, y = float(p[0]), float(p[1])
                except (ValueError, IndexError):
                    on = False; continue
                if block == 1: e.append(x); up.append(y)
                else: dn.append(y)
    if not e:
        return None
    e = np.array(e); s = np.array(up) + (np.array(dn) if len(dn) == len(up) else 0.0)
    return e, s

def regions(e, mask):
    out = []; start = None
    for i, t in enumerate(mask):
        if t and start is None: start = i
        if not t and start is not None: out.append((e[start], e[i], start, i)); start = None
    if start is not None: out.append((e[start], e[-1], start, len(e) - 1))
    return out

def analyse(d):
    path = os.path.join(RUN, d, "out_dos.log")
    if not os.path.isfile(path): return None
    try:
        r = read_total_dos(path)
    except Exception:
        return None
    if r is None: return None
    e, s = r
    key = d.split(",")[0][4:]
    z = [int(key[i:i+2]) for i in range(0, len(key), 2)]
    elems = "".join(str(Element.from_Z(v)) for v in z)
    polytyp = [x for x in d.split(",") if x.startswith("polytyp_")][0][8:]
    ew = None
    try:
        for line in open(os.path.join(RUN, d, "inputcard_go")):
            if line.startswith("#"): continue
        lines = [l for l in open(os.path.join(RUN, d, "inputcard_go")) if not l.startswith("#")]
        ew = float(lines[2].split()[1])
    except Exception:
        pass
    res = {"dir": d, "key": key, "elements": elems, "z": z, "polytyp": polytyp, "ewidth_go": ew,
           "emin": float(e.min()), "emax": float(e.max())}
    for th in TH:
        low = regions(e, s < th)
        # gaps wide enough, below EF
        gaps = [g for g in low if g[1] - g[0] >= MIN_GAP and g[1] < 0]
        found = None
        for g in sorted(gaps, key=lambda g: -g[1]):   # highest gap first
            # DOS>th region below the gap
            below = e < g[0]
            high = regions(e, (s >= th) & below)
            high = [h for h in high if h[1] - h[0] >= MIN_SEMI]
            if high:
                h = max(high, key=lambda h: h[1])
                found = {"gap": (round(float(g[0]), 4), round(float(g[1]), 4)),
                         "semicore": (round(float(h[0]), 4), round(float(h[1]), 4)),
                         "semicore_max_dos": round(float(s[h[2]:h[3] + 1].max()), 4),
                         "ewidth_in_gap": bool(ew is not None and g[0] < -ew < g[1])}
                break
        res["th_%g" % th] = found
    return res

if __name__ == "__main__":
    dirs = sorted(os.listdir(RUN))
    with Pool(20) as p:
        out = [r for r in p.map(analyse, dirs, chunksize=50) if r]
    json.dump(out, open(sys.argv[1], "w"))
    n = len(out)
    for th in TH:
        k = "th_%g" % th
        hits = [r for r in out if r[k]]
        print("th", th, ": semicore-below-gap found in", len(hits), "of", n)
