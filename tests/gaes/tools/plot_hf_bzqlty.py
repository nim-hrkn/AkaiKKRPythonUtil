"""AlSiSnHf bcc (cpa2021v01, lmxtyp=2, ewidth 1.15): bzqlty 10 / 14 / 18."""
import re, os
import matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pyakaikkr import AkaikkrJob
from pyakaikkr.gaes import dos_curves_from_outputs, levels_from_go
from pyakaikkr.gaes.ewidth import choose_ewidth2
from pyakaikkr.gaes.plot import draw_gaes_dos, legend_text
RUNS = [("bzqlty 10", "lmx2_ew_1.15"), ("bzqlty 14", "lmx2_ew_1.15_bz14"), ("bzqlty 18", "lmx2_ew_1.15_bz18")]
def info(d):
    txt = open(d + "/out_go.log").read()
    itr = re.findall(r"itr=\s*(\d+)", txt); err = re.findall(r"rms err=\s*(-?\d+\.\d+)", txt)
    core = re.search(r"\*\*\* type-\S+\s+Hf.*?core charge in the muffin-tin sphere =\s*(-?\d+\.\d+)", txt, re.S)
    return (int(itr[-1]) if itr else None), (float(err[-1]) if err else None), (float(core.group(1)) if core else None), [float(x) for x in err]
fig, axes = plt.subplots(len(RUNS), 2, figsize=(16, 4.4 * len(RUNS)))
for row, (label, d) in enumerate(RUNS):
    ax, axh = axes[row]
    if not os.path.isfile(d + "/out_dos.log"):
        ax.set_title("%s: no output yet" % label, loc="left"); continue
    e, c, _ = dos_curves_from_outputs([(AkaikkrJob(d), "out_dos.log")]); lv = levels_from_go(d + "/out_go.log")
    n, err, core, hist = info(d); conv = err is not None and err < -5.9
    dec = choose_ewidth2(e, c, 1.15, dosth=2e-2, dosth2=1e-3, eth=0.3, ediff=0.2, min_ewidth=1.0, max_ewidth=2.0)
    draw_gaes_dos(ax, e, c[0], ewidth=1.15, decision=dec, bounds=(1.0, 2.0), levels=lv, highlight=["Hf4f", "Hf5p", "Sn4d"], xlim=(-2.4, 0.9))
    roles = ", ".join("%s: %s (%.2f)" % (k, "VALENCE *" if lv[k].star else "CORE", lv[k].e) for k in ("Hf4f",) if k in lv)
    ax.set_title("AlSiSnHf bcc (cpa2021v01) lmxtyp=2, ewidth 1.15, %s: %s, itr %s, log10 rms %s\nHf core charge in MT %s; %s; gap judgement: %s" % (
        label, "CONVERGED" if conv else "NOT CONVERGED", n, err, core, roles, dec.flag), fontsize=8, loc="left", color="k" if conv else "darkred")
    axh.plot(range(1, len(hist) + 1), hist, color="k", lw=1.0); axh.set_xlabel("iteration"); axh.set_ylabel("log10 rms err")
    axh.set_title("%s: SCF history" % label, fontsize=8, loc="left"); axh.axhline(-6, color="gray", lw=0.6, ls="--")
axes[-1][0].set_xlabel("E - EF (Ry)")
fig.suptitle("AlSiSnHf bcc, ewidth 1.15, lmxtyp=2: k-point density (bzqlty). " + legend_text(), fontsize=9)
fig.tight_layout(rect=(0, 0, 1, 0.975)); fig.savefig("hf_bz_dos.png", dpi=120); print("saved hf_bz_dos.png")
