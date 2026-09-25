"""DOS of the last judgement of every orbital-rule GAES run (RUN_orb*), with the rule's level, the bounds band and the regions."""
import json, glob, os, sys, textwrap
import matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pyakaikkr import AkaikkrJob
from pyakaikkr.gaes import dos_curves_from_outputs, levels_from_go
from pyakaikkr.gaes.ewidth import choose_ewidth2
from pyakaikkr.gaes.plot import draw_gaes_dos, legend_text
# newest run of each tag wins (RUN_orb2 < ... < RUN_orb6); GLOB overrides the search pattern
# (e.g. GLOB="RUN/*/key_*.json" for plain GAES runs without orbital rules)
cases = {}
for k in sorted(glob.glob(os.environ.get("GLOB", "RUN_orb*/*/key_*.json"))):
    tag = os.path.basename(os.path.dirname(k)); cases[tag] = k
order = sys.argv[1:] or sorted(cases)
n = len(order); ncol = 2; nrow = (n + 1) // 2
fig, axes = plt.subplots(nrow, ncol, figsize=(16, 3.9 * nrow), sharex=True); axes = axes.ravel()
for ax, tag in zip(axes, order):
    d = json.load(open(cases[tag])); rules = d["parameters"].get("orbitals", []); pr = d["parameters"]
    if not d["judgements"]:
        ax.set_title("%s: %s %s" % (tag, d["status"], d["message"][:80]), fontsize=8, loc="left"); continue
    jf = d["judgements"][-1]; dd = list(jf["directories"].values())[0]
    e, c, _ = dos_curves_from_outputs([(AkaikkrJob(dd), "out_dos.log")]); lv = levels_from_go(dd + "/out_go.log")
    lo, hi = jf.get("orbital_bounds") or (pr.get("min_ewidth"), pr.get("max_ewidth"))
    dec = choose_ewidth2(e, c, jf["ewidth"], dosth=pr["dosth"], dosth2=pr["dosth2"], eth=pr["eth"], ediff=pr["ediff"],
                         margin=pr["margin"], dosth2_relax=pr["dosth2_relax"], min_ewidth=lo, max_ewidth=hi)
    # highlight the rule orbitals, or (without rules) every level inside the window
    hl = [r.split("=")[0] for r in rules] or [k for k, v in lv.items() if v.e > e[0]]
    draw_gaes_dos(ax, e, c[0], ewidth=jf["ewidth"], decision=dec, bounds=(lo, hi), levels=lv, highlight=hl, xlim=(-2.9, 0.9))
    why = ("\n" + "\n".join(textwrap.fill("why: " + r, 120) for r in dec.reasons)) if dec.flag == "fail" else ""
    conv = ",".join("%s" % v for v in d["converged"].values()) if d["converged"] else None
    ax.set_title("%s: %s, ewidth %.4f, conv=%s, bounds [%s, %s], tried %s%s" % (
        tag, d["status"], jf["ewidth"], conv, lo and round(lo, 3), hi and round(hi, 3),
        [round(t, 3) for t in d["ewidth_tried"]], why), fontsize=8, loc="left")
for ax in axes[len(order):]: ax.axis("off")
for ax in axes[-ncol:]: ax.set_xlabel("E - EF (Ry)")
fig.suptitle(os.environ.get("TITLE", "GAES with orbital rules (X-Mn-Fe-Co fcc)") + ": last judgement DOS. " + legend_text(), fontsize=9)
fig.tight_layout(rect=(0, 0, 1, 0.985)); out = os.environ.get("OUT", "orbital_dos.png"); fig.savefig(out, dpi=120); print("saved", out, len(order), "panels")
