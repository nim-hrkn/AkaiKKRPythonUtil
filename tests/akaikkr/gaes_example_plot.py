# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""figure helper of the GAES examples: last judgement DOS of a KeyResult through pyakaikkr.plot."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from pyakaikkr import AkaikkrJob  # noqa: E402
from pyakaikkr.gaes import dos_curves_from_outputs, levels_from_go, per_atom  # noqa: E402
from pyakaikkr.gaes.ewidth import choose_ewidth2  # noqa: E402
from pyakaikkr.plot import plot_gaes_dos, gaes_legend_text  # noqa: E402


def draw_key_result(ax, res, highlight=()):
    """one panel: DOS of the last judgement of `res` (KeyResult) with regions, bounds, levels and reasons."""
    if not res.judgements:
        ax.set_title("{}: {} {}".format(res.key, res.status, res.message[:80]), fontsize=8, loc="left")
        return
    j = res.judgements[-1]
    d = list(j.directories.values())[0]
    pr = res.parameters
    e, c, _ = dos_curves_from_outputs([(AkaikkrJob(d), "out_dos.log")])
    lv = levels_from_go(d + "/out_go.log")
    lo, hi = j.orbital_bounds or (pr.get("min_ewidth"), pr.get("max_ewidth"))
    natm = j.natm if j.dos_unit == "per_atom" else None
    dec = choose_ewidth2(e, per_atom(c, natm), j.ewidth, dosth=pr["dosth"], dosth2=pr["dosth2"], eth=pr["eth"],
                         ediff=pr["ediff"], margin=pr["margin"], dosth2_relax=pr["dosth2_relax"], min_ewidth=lo, max_ewidth=hi)
    hl = list(highlight) or [k for k, v in lv.items() if v.e > e[0]]
    plot_gaes_dos(ax, e, c[0], ewidth=j.ewidth, decision=dec, bounds=(lo, hi), levels=lv, highlight=hl, natm=natm,
                  dosth=pr["dosth"], dosth2=pr["dosth2"], xlim=(-2.9, 0.9))
    why = "".join("\nwhy: " + r for r in dec.reasons) if dec.flag == "fail" else ""
    ax.set_title("%s: %s, ewidth %s, converged %s, bounds [%s, %s], tried %s%s" % (
        res.key, res.status, res.ewidth_final, res.converged, lo and round(lo, 3), hi and round(hi, 3),
        [round(t, 3) for t in res.ewidth_tried], why), fontsize=8, loc="left")
    ax.set_xlabel("E - EF (Ry)")


def plot_key_result(res, png, title="", highlight=()):
    fig, ax = plt.subplots(figsize=(10, 4.2))
    draw_key_result(ax, res, highlight)
    fig.suptitle(title + "\n" + gaes_legend_text(), fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(png, dpi=120)
    plt.close(fig)
    return png


def plot_key_results(results, png, title="", highlight=()):
    """results: {label: KeyResult}; one panel per entry side by side."""
    fig, axes = plt.subplots(1, len(results), figsize=(9 * len(results), 4.6), squeeze=False)
    for ax, (label, res) in zip(axes[0], results.items()):
        draw_key_result(ax, res, highlight)
        ax.set_title(label + "\n" + ax.get_title(loc="left"), fontsize=8, loc="left")
    fig.suptitle(title + "\n" + gaes_legend_text(), fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(png, dpi=120)
    plt.close(fig)
    return png
