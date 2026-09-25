"""AlSiSnHf bcc (cpa2021v01) with lmxtyp=3 vs 2: total DOS, Hf f-PDOS, core levels, contour bottom."""
import re, os, sys
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pyakaikkr import AkaikkrJob
from pyakaikkr.gaes import dos_curves_from_outputs, levels_from_go, pdos_curves_from_output
from pyakaikkr.gaes.ewidth import choose_ewidth2
from pyakaikkr.gaes.plot import draw_gaes_dos, legend_text
C = "/tmp/claude-1000/-home-kino-kino-Claude-Project-AKAIKKR/cf1c7f21-9c6c-416d-bd70-1577df1e05cf/scratchpad/ewscan_AlSiSnHf"
T = "/home/kino/kino/Claude/Project/AKAIKKR/AkaiKKRPythonUtil/tests/gaes/RUN_Hf_stopped/Al0p25Si0p25Sn0p25Hf0p25_bcc/key_13145072,ew_000-1.2000,ed_000-1e-04,polytyp_bcc,pm_000-5e-03"
LMX2 = {0.9: C + "/ew_0.9", 1.2: T}
EWS = [0.9, 1.2, 1.5]
def info(d):
    txt = open(d + "/out_go.log").read()
    itr = re.findall(r"itr=\s*(\d+)", txt); err = re.findall(r"rms err=\s*(-?\d+\.\d+)", txt)
    tot = re.search(r"\*\*\* type-\S+\s+Hf.*?total charge=\s*(-?\d+\.\d+)", txt, re.S)
    return (int(itr[-1]) if itr else None), (float(err[-1]) if err else None), (float(tot.group(1)) if tot else None)
fig, axes = plt.subplots(len(EWS), 2, figsize=(16, 4.6 * len(EWS)), sharex=True)
for row, ew in enumerate(EWS):
    for col, (label, d) in enumerate((("lmxtyp=2", LMX2.get(ew)), ("lmxtyp=3", "ew_%s" % ew))):
        ax = axes[row, col]
        if d is None or not os.path.isfile(d + "/out_dos.log"):
            ax.set_title("%s ewidth %.1f: no dos output" % (label, ew), fontsize=9, loc="left"); continue
        e, c, _ = dos_curves_from_outputs([(AkaikkrJob(d), "out_dos.log")]); lv = levels_from_go(d + "/out_go.log")
        n, err, tot = info(d); conv = err is not None and err < -5.9
        dec = choose_ewidth2(e, c, ew, dosth=2e-2, dosth2=1e-3, eth=0.3, ediff=0.2, min_ewidth=1.0, max_ewidth=2.0)
        draw_gaes_dos(ax, e, c[0], ewidth=ew, decision=dec, bounds=(1.0, 2.0), levels=lv, highlight=["Hf4f", "Hf5p", "Sn4d"], xlim=(-2.4, 0.9))
        try:
            pe, pc, labels = pdos_curves_from_output(AkaikkrJob(d), "out_dos.log", l_sum=False)
            # PDOS columns are s, p, d[, f], total: the f column exists only with lmxtyp=3
            for lab, curve in zip(labels, pc):
                if "Hf" in lab and curve.ndim == 2:
                    if curve.shape[1] == 5:
                        ax.plot(pe, curve[:, 3] * 0.25, color="tab:purple", lw=1.0, label="Hf f-PDOS x conc")
                    else:
                        ax.text(0.99, 0.05, "no f channel (lmxtyp=2)", transform=ax.transAxes, fontsize=8, ha="right", color="tab:purple")
        except Exception as ex:
            ax.text(0.02, 0.05, "pdos: %s" % ex, transform=ax.transAxes, fontsize=7)
        roles = ", ".join("%s: %s (%.2f)" % (k, "VALENCE *" if lv[k].star else "CORE", lv[k].e) for k in ("Hf4f",) if k in lv)
        ax.set_title("AlSiSnHf bcc (cpa2021v01) %s, ewidth %.1f: %s, itr %s, log10 rms %s\nHf total charge %s (Z=72); %s; gap judgement: %s%s" % (
            label, ew, "CONVERGED" if conv else "NOT CONVERGED", n, err, tot, roles, dec.flag,
            " -> %.4f" % dec.ewidth if dec.flag == "new" else ""), fontsize=8, loc="left", color="k" if conv else "darkred")
        if ax.get_legend_handles_labels()[0]: ax.legend(fontsize=7, loc="lower right")
for ax in axes[-1]: ax.set_xlabel("E - EF (Ry)")
fig.suptitle("AlSiSnHf bcc: f channel absent (lmxtyp=2, left) vs present (lmxtyp=3, right). " + legend_text(), fontsize=9)
fig.tight_layout(rect=(0, 0, 1, 0.975)); fig.savefig("hf_lmx3_dos.png", dpi=120); print("saved hf_lmx3_dos.png")
