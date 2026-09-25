"""DOS of go runs where a small change of ewidth flips the SCF convergence (spec section 12 note)."""
import re, glob, os
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pyakaikkr import AkaikkrJob
from pyakaikkr.gaes import dos_curves_from_outputs, levels_from_go
from pyakaikkr.gaes.ewidth import choose_ewidth2
H = "/tmp/claude-1000/-home-kino-kino-Claude-Project-AKAIKKR/cf1c7f21-9c6c-416d-bd70-1577df1e05cf/scratchpad/ewscan_AlSiSnHf"
T = "/home/kino/kino/Claude/Project/AKAIKKR/AkaiKKRPythonUtil/tests/gaes/RUN_Hf_stopped/Al0p25Si0p25Sn0p25Hf0p25_bcc"
PAIRS = [
    ("AlSiSnHf bcc (cpa2021v01), Hf 4f", ["Hf4f", "Hf5p"], [(H + "/ew_0.9", 0.9), (T + "/key_13145072,ew_000-1.2000,ed_000-1e-04,polytyp_bcc,pm_000-5e-03", 1.2)]),
    ("InMnFeCo fcc, In4d=core", ["In4d"], [(glob.glob("RUN_orb4/InMnFeCo_In4d-core/key_*ew_000-0.6410*")[0], 0.641), (glob.glob("RUN_orb5/InMnFeCo_In4d-core/key_*ew_000-1.2000*")[0], 1.2)]),
    ("TlMnFeCo fcc, Tl5d=core", ["Tl5d"], [(glob.glob("RUN_orb4/TlMnFeCo_Tl5d-core/key_*ew_000-0.3860*")[0], 0.386), (glob.glob("RUN_orb5/TlMnFeCo_Tl5d-core/key_*ew_000-1.2000*")[0], 1.2)]),
]
def info(d):
    txt = open(d + "/out_go.log").read()
    itr = re.findall(r"itr=\s*(\d+)", txt); err = re.findall(r"rms err=\s*(-?\d+\.\d+)", txt)
    return (int(itr[-1]) if itr else None), (float(err[-1]) if err else None), ("cpu time" in txt or "sbrtime" in txt)
fig, axes = plt.subplots(len(PAIRS), 2, figsize=(16, 3.4 * len(PAIRS)), sharex=True)
for row, (name, keys, runs) in enumerate(PAIRS):
    for col, (d, ew) in enumerate(runs):
        ax = axes[row, col]
        e, c, _ = dos_curves_from_outputs([(AkaikkrJob(d), "out_dos.log")]); lv = levels_from_go(d + "/out_go.log")
        n, err, _ = info(d); conv = err is not None and err < -5.9
        ax.plot(e, c[0], color="k", lw=1.1)
        # gap regions by the Method 2 rules (dosth 2e-2 / dosth2 1e-3, eth 0.3, ediff 0.2), no ewidth bounds
        dec = choose_ewidth2(e, c, ew, dosth=2e-2, dosth2=1e-3, eth=0.3, ediff=0.2)
        for g in dec.coarse: ax.axvspan(g.e1, g.e2, color="tab:green", alpha=0.10)
        for g in dec.fine: ax.axvspan(g.e1, g.e2, color="tab:blue", alpha=0.15)
        flag = dec.flag
        ax.axvline(-ew, color="red", ls="-.", lw=1.4); ax.axvline(0, color="gray", lw=0.6, ls="--")
        for k in keys:
            v = lv.get(k)
            if v is not None:
                ax.axvline(v.e, color="tab:blue", lw=1.5, ls="--" if v.star else "-")
                ax.text(v.e, 100, "%s%s %.2f" % (k, "*" if v.star else "", v.e), rotation=90, fontsize=8, color="tab:blue", ha="right", va="top")
        for k, v in lv.items():
            if v.e > e[0] and k not in keys: ax.axvline(v.e, color="0.6", lw=0.8, ls="--" if v.star else "-")
        ax.set_yscale("log"); ax.set_ylim(1e-4, 3e2); ax.set_xlim(-2.9, 0.9)
        ax.set_title("%s: ewidth %.4f, %s (itr %s, log10 rms %s), gap judgement: %s" % (
            name, ew, "converged" if conv else "NOT converged", n, err, flag), fontsize=8, loc="left",
            color="k" if conv else "darkred")
for ax in axes[-1]: ax.set_xlabel("E - EF (Ry)")
fig.suptitle("SCF convergence flips with ewidth: red = -ewidth (contour bottom), blue = the level concerned (dashed = * valence), gray = other core levels; green / blue shade = coarse / fine gap regions (Method 2)", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.98)); fig.savefig("conv_pairs_dos.png", dpi=120); print("saved conv_pairs_dos.png")
