"""DOS of the last judgement of every orbital-rule GAES run (RUN_orb*), with the rule's level, the bounds band and the regions."""
import json, glob, os, sys
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pyakaikkr import AkaikkrJob
from pyakaikkr.gaes import dos_curves_from_outputs, levels_from_go
from pyakaikkr.gaes.ewidth import choose_ewidth2
# newest run of each tag wins (RUN_orb2 < ... < RUN_orb6)
cases = {}
for k in sorted(glob.glob("RUN_orb*/*/key_*.json")):
    tag = os.path.basename(os.path.dirname(k)); cases[tag] = k
order = sys.argv[1:] or sorted(cases)
n = len(order); ncol = 2; nrow = (n + 1) // 2
fig, axes = plt.subplots(nrow, ncol, figsize=(16, 3.3 * nrow), sharex=True); axes = axes.ravel()
for ax, tag in zip(axes, order):
    d = json.load(open(cases[tag])); rules = d["parameters"].get("orbitals", []); pr = d["parameters"]
    jf = d["judgements"][-1]; dd = list(jf["directories"].values())[0]
    e, c, _ = dos_curves_from_outputs([(AkaikkrJob(dd), "out_dos.log")]); lv = levels_from_go(dd + "/out_go.log")
    lo, hi = jf["orbital_bounds"]
    dec = choose_ewidth2(e, c, jf["ewidth"], dosth=pr["dosth"], dosth2=pr["dosth2"], eth=pr["eth"], ediff=pr["ediff"],
                         margin=pr["margin"], dosth2_relax=pr["dosth2_relax"], min_ewidth=lo, max_ewidth=hi)
    ax.plot(e, c[0], color="k", lw=1.1)
    for g in dec.coarse: ax.axvspan(g.e1, g.e2, color="tab:green", alpha=0.10)
    for g in dec.fine: ax.axvspan(g.e1, g.e2, color="tab:blue", alpha=0.15)
    ax.axvspan(-hi if hi is not None else e[0], -lo if lo is not None else 0, facecolor="none", edgecolor="0.4", hatch="///", lw=0, alpha=0.25)
    ax.axvline(-jf["ewidth"], color="red", ls="-.", lw=1.3)
    ax.axvline(0, color="gray", lw=0.6, ls="--")
    for r in rules:
        key = r.split("=")[0]; v = lv.get(key)
        if v is not None and v.e > e[0]:
            ax.axvline(v.e, color="tab:blue", lw=1.5, ls="--" if v.star else "-")
            ax.text(v.e, 100, "%s%s %.2f" % (key, "*" if v.star else "", v.e), rotation=90, fontsize=8, color="tab:blue", ha="right", va="top")
    for key, v in lv.items():
        if v.e > e[0] and key not in [r.split("=")[0] for r in rules]:
            ax.axvline(v.e, color="0.6", lw=0.8, ls="--" if v.star else "-")
    ax.set_yscale("log"); ax.set_ylim(1e-4, 3e2); ax.set_xlim(-2.9, 0.9)
    ax.set_title("%s: %s, ewidth %.4f, conv=%s, bounds [%s, %s], tried %s" % (
        tag, d["status"], jf["ewidth"], d["converged"].get("fcc"), lo and round(lo, 3), hi and round(hi, 3),
        [round(t, 3) for t in d["ewidth_tried"]]), fontsize=8, loc="left")
for ax in axes[len(order):]: ax.axis("off")
for ax in axes[-ncol:]: ax.set_xlabel("E - EF (Ry)")
fig.suptitle("GAES with orbital rules (X-Mn-Fe-Co fcc): last judgement DOS; blue line = the rule's level (dashed = *), hatched = [min, max] ewidth, red = -ewidth", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.985)); out = os.environ.get("OUT", "orbital_dos.png"); fig.savefig(out, dpi=120); print("saved", out, len(order), "panels")
