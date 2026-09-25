# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""matplotlib helpers shared by every GAES DOS figure: the gap regions (green coarse / blue fine),
the [min_ewidth, max_ewidth] band (hatched), the contour bottom E_F - ewidth_go (red) and the
core levels of the components (solid core, dashed '*' valence).

    from pyakaikkr.gaes.plot import draw_gaes_dos
    draw_gaes_dos(ax, energy, dos, ewidth=1.2, decision=dec, bounds=(1.0, 2.0), levels=levels_from_go(out_go))
"""

COARSE_COLOR, COARSE_ALPHA = "tab:green", 0.10
FINE_COLOR, FINE_ALPHA = "tab:blue", 0.15
BOUNDS_COLOR = "0.4"
EWIDTH_COLOR = "red"
LEVEL_COLOR, LEVEL_COLOR_OTHER = "tab:blue", "0.6"


def shade_regions(ax, coarse=(), fine=()):
    """green bands for the coarse gap regions, blue for the fine sub-regions (GapRegion or (e1, e2))."""
    for g in coarse:
        e1, e2 = (g.e1, g.e2) if hasattr(g, "e1") else g[:2]
        ax.axvspan(e1, e2, color=COARSE_COLOR, alpha=COARSE_ALPHA, lw=0)
    for g in fine:
        e1, e2 = (g.e1, g.e2) if hasattr(g, "e1") else g[:2]
        ax.axvspan(e1, e2, color=FINE_COLOR, alpha=FINE_ALPHA, lw=0)


def shade_ewidth_bounds(ax, bounds, e_min=None, label=True):
    """hatched band E_F - max_ewidth .. E_F - min_ewidth (None = open side: e_min / E_F)."""
    if bounds is None:
        return
    lo, hi = bounds[0], bounds[1]
    if lo is None and hi is None:
        return
    left = -hi if hi is not None else (e_min if e_min is not None else ax.get_xlim()[0])
    right = -lo if lo is not None else 0.0
    ax.axvspan(left, right, facecolor="none", edgecolor=BOUNDS_COLOR, hatch="///", lw=0, alpha=0.25,
               label="[min, max] ewidth = [%s, %s]" % (lo, hi) if label else None)


def draw_ewidth(ax, ewidth, label=None, **kw):
    """the contour bottom E_F - ewidth_go (red dash-dot); ewidth may be None."""
    if ewidth is None:
        return
    kw = dict(color=EWIDTH_COLOR, ls="-.", lw=1.6, label=label, zorder=6) | kw   # on top of level lines at the same energy
    ax.axvline(-ewidth, **kw)


def draw_levels(ax, levels, highlight=(), e_min=None, y_text=100.0, fontsize=8, label_all=True):
    """core levels E - E_F from levels_from_go (dict key -> Level with .e / .star) or a dict key -> [e, star]:
    solid = core, dashed = '*' (valence). Keys in `highlight` are blue; the others gray. Every level inside
    the window is labelled (element, orbital, '*', E - E_F) unless label_all=False."""
    for key, v in (levels or {}).items():
        e, star = (v.e, v.star) if hasattr(v, "e") else (v[0], v[1])
        if e_min is not None and e <= e_min:
            continue
        hi = key in highlight
        color = LEVEL_COLOR if hi else LEVEL_COLOR_OTHER
        ax.axvline(e, color=color, lw=1.5 if hi else 0.9, ls="--" if star else "-", alpha=1.0 if hi else 0.9)
        if hi or label_all:
            ax.text(e, y_text, "%s%s %.2f" % (key, "*" if star else "", e), rotation=90, fontsize=fontsize if hi else fontsize - 1,
                    color=color, ha="right", va="top")


def draw_gaes_dos(ax, energy, dos, ewidth=None, decision=None, coarse=(), fine=(), bounds=None, levels=None,
                  highlight=(), log=True, ylim=(1e-4, 3e2), xlim=None):
    """one GAES DOS panel: curve + regions + bounds band + contour bottom + core levels + E_F line.

    decision: pyakaikkr.gaes.Decision (its coarse / fine are used unless coarse / fine are given).
    bounds: (min_ewidth, max_ewidth) or None. levels: levels_from_go(...) dict (or {key: [e, star]}).
    """
    ax.plot(energy, dos, color="k", lw=1.1)
    if decision is not None and not coarse and not fine:
        coarse, fine = decision.coarse, decision.fine
    shade_regions(ax, coarse, fine)
    e_min = float(min(energy))
    shade_ewidth_bounds(ax, bounds, e_min=e_min, label=False)
    draw_ewidth(ax, ewidth)
    ax.axvline(0, color="gray", lw=0.6, ls="--")
    draw_levels(ax, levels, highlight=highlight, e_min=e_min, y_text=ylim[1] / 3)
    if log:
        ax.set_yscale("log")
    ax.set_ylim(*ylim)
    if xlim is not None:
        ax.set_xlim(*xlim)


def legend_text():
    return ("black = DOS, green / blue = coarse / fine gap regions, hatched = [min, max] ewidth, "
            "red = -ewidth_go, blue lines = core levels (dashed = * valence), gray = other levels")


__all__ = ["draw_gaes_dos", "shade_regions", "shade_ewidth_bounds", "draw_ewidth", "draw_levels", "legend_text"]
