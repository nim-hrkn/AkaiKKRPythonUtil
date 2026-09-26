"""Plot atomic core levels (from atomic_levels_dsp.py) against Z: full range (symlog) and shallow range.

Usage: python plot_atomic_levels.py docs/data/atomic_core_levels_dsp.csv -o docs/data/atomic_core_levels_dsp.png
"""
import argparse
import csv
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ORDER = ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "6s", "6p"]
COLORS = {"s": "tab:blue", "p": "tab:orange", "d": "tab:green", "f": "tab:red"}
STYLE = {"1": "o", "2": "s", "3": "^", "4": "D", "5": "v", "6": "P"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("-o", "--output", default="atomic_core_levels_dsp.png")
    ap.add_argument("--ef", type=float, default=0.6, help="assumed E_F (Ry, absolute) for the right-hand E - E_F scale")
    ap.add_argument("--ewidth", type=float, default=1.2, help="ewidth of go for the contour-bottom guide line")
    a = ap.parse_args()
    ef = a.ef
    rows = [r for r in csv.DictReader(open(a.csv)) if r["orbital"]]
    by = defaultdict(list)
    for r in rows:
        by[r["orbital"]].append((int(r["Z"]), float(r["energy_Ry"]), r["star"] == "1"))
    fig, axes = plt.subplots(2, 1, figsize=(12, 10), gridspec_kw={"height_ratios": [1.2, 1]})
    for ax, (lo, hi, scale, title) in zip(axes, ((-3200, 0, "symlog", "all core levels (symlog scale)"),
                                                 (-6.0, 1.0, "linear", "shallow core levels, E > -6 Ry (semicore candidates)"))):
        for orb in ORDER:
            pts = [p for p in by.get(orb, []) if lo <= p[1] <= hi]
            if not pts:
                continue
            z = [p[0] for p in pts]
            e = [p[1] for p in pts]
            ax.plot(z, e, marker=STYLE[orb[0]], ms=4, lw=0.8, color=COLORS[orb[1]], label=orb, alpha=0.9)
            stars = [(p[0], p[1]) for p in pts if p[2]]
            if stars:
                ax.scatter([s[0] for s in stars], [s[1] for s in stars], s=90, facecolors="none", edgecolors="k",
                           linewidths=1.2, zorder=5, label="_nolegend_")
        if scale == "symlog":
            ax.set_yscale("symlog", linthresh=1.0)
        ax.set_ylim(lo, hi)
        ax.set_xlim(0, 84)
        ax.set_xlabel("Z")
        ax.set_ylabel("core level (Ry, absolute, initial atomic potential)")
        ax.set_title(title, loc="left", fontsize=10)
        ax.grid(True, lw=0.4, alpha=0.5)
        ax.axhline(ef, color="tab:red", lw=0.8, ls=":")
        ax.axhline(ef - a.ewidth, color="gray", lw=0.9, ls="--")
        ax.legend(ncol=8, fontsize=8, loc="lower left")
        if scale == "linear":
            # right-hand scale: E - E_F for the assumed E_F (typical E_F of the test materials: 0.5-0.7 Ry)
            ax2 = ax.secondary_yaxis("right", functions=(lambda e: e - ef, lambda x: x + ef))
            ax2.set_ylabel("$E - E_F$ (Ry), assuming $E_F$ = %.1f Ry" % ef)
        else:
            ax.text(83, ef - a.ewidth - 0.15, "$E_F$ - ewidth (E_F = %.1f, ewidth = %.1f)" % (ef, a.ewidth),
                    fontsize=8, color="gray", ha="right", va="top")
    axes[1].text(1, ef + 0.05, "assumed $E_F$ = %.1f Ry" % ef, fontsize=8, color="tab:red")
    axes[1].text(1, ef - a.ewidth - 0.25, "$E_F$ - ewidth = %.1f Ry (contour bottom for ewidth %.1f)" % (ef - a.ewidth, a.ewidth),
                 fontsize=8, color="gray")
    for z, sym in ((13, "Al"), (21, "Sc"), (31, "Ga"), (32, "Ge"), (39, "Y"), (40, "Zr"), (48, "Cd"), (49, "In"), (50, "Sn"),
                   (72, "Hf"), (80, "Hg"), (81, "Tl"), (82, "Pb"), (83, "Bi")):
        axes[1].axvline(z, color="lightgray", lw=0.6)
        axes[1].text(z, 0.9, sym, fontsize=7, ha="center", color="dimgray")
    fig.suptitle("Atomic core levels from specx go=dsp (fresh potential, fcc, a = experimental volume, nmag, pbe); "
                 "black circles: moved to valence (*); right scale: E - E_F for E_F = %.1f Ry" % ef, fontsize=10)
    fig.tight_layout()
    fig.savefig(a.output, dpi=140)
    print("saved", a.output)


if __name__ == "__main__":
    main()
