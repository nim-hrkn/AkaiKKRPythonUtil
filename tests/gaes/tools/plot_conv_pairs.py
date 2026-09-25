"""DOS of go runs where a small change of ewidth flips the SCF convergence (spec section 12 note)."""
import re, glob
import matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pyakaikkr import AkaikkrJob
from pyakaikkr.gaes import dos_curves_from_outputs, levels_from_go
from pyakaikkr.gaes.ewidth import choose_ewidth2
from pyakaikkr.gaes.plot import draw_gaes_dos, legend_text
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
fig, axes = plt.subplots(len(PAIRS), 2, figsize=(16, 3.8 * len(PAIRS)), sharex=True)
for row, (name, keys, runs) in enumerate(PAIRS):
    for col, (d, ew) in enumerate(runs):
        ax = axes[row, col]
        e, c, _ = dos_curves_from_outputs([(AkaikkrJob(d), "out_dos.log")]); lv = levels_from_go(d + "/out_go.log")
        n, err, _ = info(d); conv = err is not None and err < -5.9
        # gap regions by the Method 2 rules (dosth 2e-2 / dosth2 1e-3, eth 0.3, ediff 0.2), no ewidth bounds
        dec = choose_ewidth2(e, c, ew, dosth=2e-2, dosth2=1e-3, eth=0.3, ediff=0.2)
        flag = dec.flag
        draw_gaes_dos(ax, e, c[0], ewidth=ew, decision=dec, levels=lv, highlight=keys, xlim=(-2.9, 0.9))
        # how specx treated the orbital concerned in this go: '*' = valence (inside the contour), else core
        roles = ", ".join("%s: %s (%.2f Ry)" % (k, "VALENCE *" if lv[k].star else "CORE", lv[k].e) for k in keys if k in lv)
        ax.set_title("%s: ewidth %.4f, %s (itr %s, log10 rms %s), gap judgement: %s\n%s" % (
            name, ew, "converged" if conv else "NOT converged", n, err, flag, roles), fontsize=8, loc="left",
            color="k" if conv else "darkred")
for ax in axes[-1]: ax.set_xlabel("E - EF (Ry)")
fig.suptitle("SCF convergence flips with ewidth (second title line: how specx treated the orbital in that go, CORE or VALENCE *). " + legend_text(), fontsize=9)
fig.tight_layout(rect=(0, 0, 1, 0.98)); fig.savefig("conv_pairs_dos.png", dpi=120); print("saved conv_pairs_dos.png")
