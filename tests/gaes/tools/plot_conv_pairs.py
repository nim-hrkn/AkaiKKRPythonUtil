"""DOS of go runs where a small change of ewidth flips the SCF convergence (spec section 12 note)."""
import re, glob, os, textwrap
import matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pyakaikkr import AkaikkrJob
from pyakaikkr.gaes import dos_curves_from_outputs, levels_from_go
from pyakaikkr.gaes.ewidth import choose_ewidth2
from pyakaikkr.gaes.plot import draw_gaes_dos, legend_text
H = "/tmp/claude-1000/-home-kino-kino-Claude-Project-AKAIKKR/cf1c7f21-9c6c-416d-bd70-1577df1e05cf/scratchpad/ewscan_AlSiSnHf"
T = "/home/kino/kino/Claude/Project/AKAIKKR/AkaiKKRPythonUtil/tests/gaes/RUN_Hf_stopped/Al0p25Si0p25Sn0p25Hf0p25_bcc"
import json
def bounds_of(run_dir, ew):
    """[min_ewidth, max_ewidth] of the GAES judgement that used this go (key_*.json of the parent RUN directory)."""
    for k in glob.glob(os.path.join(os.path.dirname(run_dir), "key_*.json")):
        for j in json.load(open(k))["judgements"]:
            if abs(j["ewidth"] - ew) < 1e-6 and j.get("orbital_bounds"):
                return tuple(j["orbital_bounds"])
    return None
PAIRS = [
    # the Hf runs had no orbital rule: the default range [1.0, 2.0] applies
    ("AlSiSnHf bcc (cpa2021v01), Hf 4f", ["Hf4f", "Hf5p"], [(H + "/ew_0.9", 0.9, (1.0, 2.0)), (T + "/key_13145072,ew_000-1.2000,ed_000-1e-04,polytyp_bcc,pm_000-5e-03", 1.2, (1.0, 2.0))]),
    # the 0.641 / 0.386 starts came from the level table (step-0 rule of the first version, since replaced:
    # the first go now always uses ewidth_init); they were never chosen from a DOS
    ("InMnFeCo fcc, In4d=core", ["In4d"], [(d, ew, bounds_of(d, ew)) for d, ew in ((glob.glob("RUN_orb4/InMnFeCo_In4d-core/key_*ew_000-0.6410*")[0], 0.641), (glob.glob("RUN_orb5/InMnFeCo_In4d-core/key_*ew_000-1.2000*")[0], 1.2))]),
    ("TlMnFeCo fcc, Tl5d=core", ["Tl5d"], [(d, ew, bounds_of(d, ew)) for d, ew in ((glob.glob("RUN_orb4/TlMnFeCo_Tl5d-core/key_*ew_000-0.3860*")[0], 0.386), (glob.glob("RUN_orb5/TlMnFeCo_Tl5d-core/key_*ew_000-1.2000*")[0], 1.2))]),
]
START_NOTE = {0.9: "manual scan, not chosen by GAES", 0.641: "old step-0 start from the table, not from a DOS",
              0.386: "old step-0 start from the table, not from a DOS", 1.2: "ewidth_init"}
def info(d):
    txt = open(d + "/out_go.log").read()
    itr = re.findall(r"itr=\s*(\d+)", txt); err = re.findall(r"rms err=\s*(-?\d+\.\d+)", txt)
    return (int(itr[-1]) if itr else None), (float(err[-1]) if err else None), ("cpu time" in txt or "sbrtime" in txt)
fig, axes = plt.subplots(len(PAIRS), 2, figsize=(16, 5.0 * len(PAIRS)), sharex=True)
for row, (name, keys, runs) in enumerate(PAIRS):
    for col, (d, ew, bounds) in enumerate(runs):
        ax = axes[row, col]
        e, c, _ = dos_curves_from_outputs([(AkaikkrJob(d), "out_dos.log")]); lv = levels_from_go(d + "/out_go.log")
        n, err, _ = info(d); conv = err is not None and err < -5.9
        # gap regions by the Method 2 rules (dosth 2e-2 / dosth2 1e-3, eth 0.3, ediff 0.2), no ewidth bounds
        dec = choose_ewidth2(e, c, ew, dosth=2e-2, dosth2=1e-3, eth=0.3, ediff=0.2,
                             min_ewidth=bounds[0] if bounds else None, max_ewidth=bounds[1] if bounds else None)
        flag = dec.flag
        draw_gaes_dos(ax, e, c[0], ewidth=ew, decision=dec, bounds=bounds, levels=lv, highlight=keys, xlim=(-2.9, 0.9))
        # how specx treated the orbital concerned in this go: '*' = valence (inside the contour), else core
        roles = ", ".join("%s: %s (%.2f Ry)" % (k, "VALENCE *" if lv[k].star else "CORE", lv[k].e) for k in keys if k in lv)
        bl = "[%s, %s]" % tuple(None if b is None else round(b, 3) for b in bounds) if bounds else "none"
        why = ("\n" + "\n".join(textwrap.fill("why: " + r, 110) for r in dec.reasons)) if flag == "fail" else ""
        ax.set_title("%s: ewidth %.4f (%s)\n%s, itr %s, log10 rms %s\n%s\ngap judgement with [min, max] ewidth %s: %s%s%s" % (
            name, ew, START_NOTE.get(ew, ""), "CONVERGED" if conv else "NOT CONVERGED", n, err, roles, bl, flag,
            " -> %.4f" % dec.ewidth if flag == "new" else "", why), fontsize=8, loc="left",
            color="k" if conv else "darkred")
for ax in axes[-1]: ax.set_xlabel("E - EF (Ry)")
fig.suptitle("SCF convergence flips with ewidth. Second title line: how specx treated the orbital in that go (CORE, or VALENCE *).\n" + legend_text(), fontsize=9)
fig.tight_layout(rect=(0, 0, 1, 0.965)); fig.savefig("conv_pairs_dos.png", dpi=120); print("saved conv_pairs_dos.png")
