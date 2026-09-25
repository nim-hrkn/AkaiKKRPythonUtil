"""DOS of the final GAES solution of each X-Mn-Fe-Co system with the component core levels overlaid."""
import json, re, os, csv, glob, sys
import matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pyakaikkr import AkaikkrJob
from pyakaikkr.gaes import dos_curves_from_outputs
from pyakaikkr.gaes.ewidth import choose_ewidth2
from pyakaikkr.gaes.plot import draw_gaes_dos, legend_text
tableA = {}
for r in csv.DictReader(open("/home/kino/kino/Claude/Project/AKAIKKR/AkaiKKRPythonUtil/docs/data/converged_core_levels_2019.csv")):
    tableA.setdefault(r["element"], []).append((r["orbital"], float(r["E_minus_EF_median_Ry"])))
COL = {"Mn": "tab:pink", "Fe": "tab:orange", "Co": "gold"}
def own_levels(d):
    txt = open(d + "/out_go.log").read(); m = re.findall(r"\n\s+ef=\s*([-\d.]+)", txt)
    ef = float(m[-1]) if m else float("nan"); out = []
    for b in re.split(r"\*\*\* type-", txt)[1:]:
        mm = re.match(r"\S+\s+(\S+)\s+\(z=", b)
        if not mm: continue
        seen = set()
        for v, l, st in re.findall(r"(-?\d+\.\d+) Ry\((\d[spdf])\)(\*?)", b):
            e = float(v) - ef
            if e > -2.6 and (mm.group(1), l) not in seen: out.append((mm.group(1), l, e, bool(st))); seen.add((mm.group(1), l))
    return ef, out
summary = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "summary.json"))
names = [n for n in summary if summary[n].get("status") not in ("exception",)]
rows = []
n = len(names); ncol = 2; nrow = (n + 1) // 2
fig, axes = plt.subplots(nrow, ncol, figsize=(16, 3.4 * nrow), sharex=True); axes = axes.ravel()
for ax, name in zip(axes, names):
    x = name[:-6]; s = summary[name]
    js = glob.glob("RUN/%s_fcc/key_%s.json" % (name, name))
    if not js: ax.set_title(name + ": no result"); continue
    d = json.load(open(js[0]))
    if not d["judgements"]:
        ax.set_title("%s: %s %s" % (name, d["status"], d["message"][:60]), fontsize=9); continue
    jf = d["judgements"][-1]; dd = list(jf["directories"].values())[0]
    if not os.path.isfile(dd + "/out_dos.log"):
        ax.set_title("%s: %s (no dos)" % (name, d["status"]), fontsize=9); continue
    e, c, _ = dos_curves_from_outputs([(AkaikkrJob(dd), "out_dos.log")]); ef, lv = own_levels(dd)
    pr = d.get("parameters", {})
    dec = choose_ewidth2(e, c, jf["ewidth"], dosth=pr.get("dosth", 2e-2), dosth2=pr.get("dosth2", 1e-3), eth=pr.get("eth", 0.3),
                         ediff=pr.get("ediff", 0.2), margin=pr.get("margin", 0.01), dosth2_relax=pr.get("dosth2_relax", 2.0),
                         min_ewidth=pr.get("min_ewidth"), max_ewidth=pr.get("max_ewidth"))
    if dec.flag != jf["flag"] or [round(v, 4) for v in dec.candidates] != [round(v, 4) for v in jf.get("candidates", [])]:
        print("re-judged %s: stored %s %s -> now %s %s" % (name, jf["flag"], jf.get("candidates"), dec.flag, dec.candidates))
    # every level of the system: the X levels highlighted (blue), Mn / Fe / Co in their own colours below
    draw_gaes_dos(ax, e, c[0], ewidth=d["ewidth_final"], decision=dec, bounds=(pr.get("min_ewidth"), pr.get("max_ewidth")),
                  levels={el + l: [en, st] for el, l, en, st in lv if el == x}, highlight=[el + l for el, l, en, st in lv if el == x],
                  xlim=(-2.4, 0.8))
    for el, l, en, st in lv:
        if el != x:
            ax.axvline(en, color=COL.get(el, "k"), lw=1.3, ls="--" if st else "-", alpha=0.9)
            ax.text(en, 100, "%s %s%s %.2f" % (el, l, "*" if st else "", en), rotation=90, fontsize=7, color=COL.get(el, "k"), ha="right", va="top")
    for l, en in tableA.get(x, []): ax.axvline(en, color="tab:blue", lw=0.8, ls=":", alpha=0.7)
    conv = d["converged"].get("fcc")
    ax.set_title("%s fcc: %s, ewidth %.4f, conv=%s, E_F=%.3f, tried %s" % (name, d["status"], d["ewidth_final"] or 0, conv, ef,
                 [round(t, 3) for t in d["ewidth_tried"]]), fontsize=8, loc="left")
    rows.append((x, name, d["status"], d["ewidth_final"], conv, ef, [(el, l, round(en, 3), st) for el, l, en, st in lv if el == x],
                 d["gap_used"], [j["flag"] for j in d["judgements"]]))
for ax in axes[len(names):]: ax.axis("off")
for ax in axes[-ncol:]: ax.set_xlabel("E - EF (Ry)")
fig.suptitle("X-Mn-Fe-Co fcc HEAs: final GAES DOS (Method 2). " + legend_text() + "; dotted = 2019 table", fontsize=9)
fig.tight_layout(rect=(0, 0, 1, 0.985)); fig.savefig("hea_3d_dos.png", dpi=120); print("saved hea_3d_dos.png")
print("%-4s %-10s %-16s %-7s %-5s %-6s %s" % ("X", "system", "status", "ewidth", "conv", "EF", "X core levels (E-EF) | gap_used | flags"))
for r in rows: print("%-4s %-10s %-16s %-7s %-5s %-6.3f %s | %s | %s" % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8]))
