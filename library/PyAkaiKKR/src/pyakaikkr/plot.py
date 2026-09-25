# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""Array-based drawing functions shared by every AkaiKKR figure ("plot function A").

Every function here takes plain numpy arrays / dicts and a matplotlib Axes and draws on it.
Nothing here reads files or AiiDA nodes: pyakaikkr (DosPlotter, AwkPlotter, JijPlotter,
gaes.plot) reads out_*.log into arrays and calls these, and aiida-akaikkr passes the arrays of
its ArrayData / Dict nodes directly. Saving the figure is the caller's job.

Array conventions (the ones of AkaikkrJob.get_dos_as_list / get_pdos_as_list and of the
aiida-akaikkr parser):

    dos   [nspin, ne]                 total DOS, spin resolved (nspin 1 or 2)
    pdos  [nspin, ncomp, ne, nl]      per component and l (nl may be padded with NaN)
    awk   [nk, ne]                    A(w, k); kdist [nk], energy [ne], kcrt: indices of the labelled k
    jij   columns distance, J_ij(meV) of one (type, component) pair

`style` is a dict overriding DEFAULT_STYLE (colours and line widths); aiida-akaikkr passes its
own palette so the look of its figures does not change.
"""
import numpy as np

L_NAMES = ["s", "p", "d", "f", "g", "h", "i", "j", "k", "l", "m"]

DEFAULT_STYLE = {
    "series": ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown"],
    "ink": "black",           # text / main lines
    "ink2": "gray",           # secondary lines (E_F, zero, thresholds, other levels)
    "grid": "#d9d8d3",
    "ewidth": "tab:red",      # E_F - ewidth_go (contour bottom)
    "final": "tab:orange",    # final ewidth of a GAES run
    "coarse": "tab:green",    # coarse gap regions
    "fine": "tab:blue",       # fine gap sub-regions
    "bounds": "0.4",          # [min_ewidth, max_ewidth] band (hatched)
    "level": "tab:blue",      # highlighted core levels
    "level_other": "0.6",     # other core levels
    "linewidth": 1.2,
    "fill_alpha": 0.12,
    "coarse_alpha": 0.10,
    "fine_alpha": 0.15,
}

EWIDTH_LINE_LABEL = "$-$ewidth (go)"


def _style(style):
    s = dict(DEFAULT_STYLE)
    if style:
        s.update(style)
    return s


def _series(s, i):
    return s["series"][i % len(s["series"])]


# ---------------------------------------------------------------- common marks
def mark_ewidth_go(ax, energy, ewidth_go, style=None, label=EWIDTH_LINE_LABEL, **kw):
    """dash-dot vertical line at E - E_F = -|ewidth_go| (bottom of the go SCF contour).

    ewidth_go is the ewidth of the go run, never of the dos run (the dos ewidth only sets the
    energy mesh). Nothing is drawn for None. The x range is widened when the line lies below
    the mesh. Returns True when drawn.
    """
    if ewidth_go is None:
        return False
    s = _style(style)
    ebtm = -abs(float(ewidth_go))
    kw = dict(color=s["ewidth"], lw=1.3, ls="-.", label=label, zorder=6) | kw
    ax.axvline(ebtm, **kw)
    energy = np.asarray(energy, dtype=float)
    if energy.size:
        lo, hi = float(energy.min()), float(energy.max())
        if ebtm < lo:
            margin = 0.02 * (hi - ebtm)
            ax.set_xlim(ebtm - margin, hi + margin)
    return True


def mark_efermi(ax, style=None, **kw):
    """dashed vertical line at E - E_F = 0."""
    s = _style(style)
    ax.axvline(0.0, **(dict(color=s["ink2"], lw=0.8, ls="--") | kw))


def style_axes(ax, xlabel=None, ylabel=None, title=None, style=None, grid=True, title_loc="left", title_size=11):
    """labels, title, light grid and hidden top/right spines (the aiida-akaikkr look; harmless elsewhere)."""
    s = _style(style)
    if xlabel is not None:
        ax.set_xlabel(xlabel, color=s["ink"])
    if ylabel is not None:
        ax.set_ylabel(ylabel, color=s["ink"])
    if title is not None:
        ax.set_title(title, loc=title_loc, color=s["ink"], fontsize=title_size)
    if grid:
        ax.grid(True, color=s["grid"], linewidth=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(s["grid"])
    ax.tick_params(colors=s["ink2"])


# ---------------------------------------------------------------- DOS / PDOS
def plot_dos(ax, energy, dos, *, ewidth_go=None, ewidth_label=EWIDTH_LINE_LABEL, style=None,
             spin_labels=("up", "down"), mirror_down=True, fill=False, yscale=None, efermi=True):
    """total DOS on one Axes.

    Args:
        energy: [ne] E - E_F (Ry).
        dos: [nspin, ne] (or [ne] for one spin).
        ewidth_go: ewidth of the go run for the contour-bottom line (None: no line).
        mirror_down: draw the down spin negative (one panel); False draws both positive.
        fill: shade under the curve of the first spin.
        yscale: "log" / "linear" / None (leave as is). With "log" the down spin is not mirrored.
    Returns:
        list of Line2D of the spin curves.
    """
    s = _style(style)
    energy = np.asarray(energy, dtype=float)
    dos = np.asarray(dos, dtype=float)
    if dos.ndim == 1:
        dos = dos[None, :]
    if yscale == "log":
        mirror_down = False
    lines = []
    for ispin in range(dos.shape[0]):
        y = dos[ispin]
        sign = -1.0 if (ispin == 1 and mirror_down) else 1.0
        label = spin_labels[ispin] if dos.shape[0] > 1 and ispin < len(spin_labels) else None
        lines += ax.plot(energy, sign * y, color=_series(s, ispin), lw=s["linewidth"] + 0.4, label=label)
        if fill and ispin == 0:
            ax.fill_between(energy, sign * y, color=_series(s, ispin), alpha=s["fill_alpha"], lw=0)
    if dos.shape[0] > 1 and mirror_down:
        ax.axhline(0.0, color=s["ink2"], lw=0.6)
    if yscale:
        ax.set_yscale(yscale)
    if efermi:
        mark_efermi(ax, s)
    mark_ewidth_go(ax, energy, ewidth_go, s, label=ewidth_label)
    return lines


def plot_pdos(ax, energy, pdos, *, nl=None, ewidth_go=None, ewidth_label=EWIDTH_LINE_LABEL, style=None,
              l_names=L_NAMES, mirror_down=True, yscale=None, efermi=True, spin=None):
    """PDOS of one component on one Axes: one curve per l (down spin mirrored, or only `spin`).

    Args:
        energy: [ne].
        pdos: [nspin, ne, nl] (or [ne, nl] for one spin). Columns beyond `nl`, or all-NaN
            columns (padding of a type with a smaller mxl), are skipped.
        spin: None (both), 0 or 1 to draw one spin positive.
    Returns:
        list of Line2D (one per l of the first drawn spin).
    """
    s = _style(style)
    energy = np.asarray(energy, dtype=float)
    pdos = np.asarray(pdos, dtype=float)
    if pdos.ndim == 2:
        pdos = pdos[None, :, :]
    if yscale == "log":
        mirror_down = False
    nspin = pdos.shape[0]
    nl_eff = pdos.shape[2] if nl is None else min(nl, pdos.shape[2])
    spins = range(nspin) if spin is None else [spin]
    lines = []
    for il in range(nl_eff):
        if np.all(np.isnan(pdos[:, :, il])):
            continue
        color = _series(s, il)
        label = l_names[il] if il < len(l_names) else "l={}".format(il)
        for k, ispin in enumerate(spins):
            sign = -1.0 if (ispin == 1 and mirror_down and spin is None) else 1.0
            ln = ax.plot(energy, sign * pdos[ispin, :, il], color=color, lw=s["linewidth"] + 0.4,
                         label=label if k == 0 else None)
            if k == 0:
                lines += ln
    if nspin > 1 and spin is None and mirror_down:
        ax.axhline(0.0, color=s["ink2"], lw=0.6)
    if yscale:
        ax.set_yscale(yscale)
    if efermi:
        mark_efermi(ax, s)
    mark_ewidth_go(ax, energy, ewidth_go, s, label=ewidth_label)
    return lines


def component_names(type_of_site, long=True):
    """component names in PDOS order from AkaikkrJob.get_type_of_site (also stored in aiida results).

    long=True: "Rh 50.0% in HEA" style (component in type); False: the short name as written
    ("HEA_Rh_50.0%" -> "Rh").
    """
    names = []
    for t in type_of_site:
        for c in t["comp_shortname"]:
            comp = c[len(t["type"]) + 1:] if c.startswith(t["type"] + "_") else c
            if long:
                names.append("{} in {}".format(comp.replace("_", " "), t["type"]) if comp != c else c)
            else:
                names.append(comp.split("_")[0] if comp != c else c)
    return names


# ---------------------------------------------------------------- A(w, k)
def plot_awk(ax, kdist, energy, awk, kcrt, klabel=None, *, style=None, cmap="Blues", vmax_percentile=99.5,
             vmin=0.0, vmax=None, show_kgrid=True, show_ef=True, linecolor=None, linewidth=0.6, **pcolorargs):
    """Bloch spectral function A(w, k) as a colour mesh; returns the QuadMesh (for a colorbar).

    Args:
        kdist: [nk] cumulative k distance; energy: [ne]; awk: [nk, ne]; kcrt: indices in kdist of
        the labelled k points; klabel: their labels (None: no tick labels).
    """
    s = _style(style)
    kdist = np.asarray(kdist, dtype=float)
    energy = np.asarray(energy, dtype=float)
    awk = np.asarray(awk, dtype=float)
    kcrt = np.asarray(kcrt, dtype=int)
    if vmax is None and vmax_percentile is not None:
        vmax = float(np.percentile(awk, vmax_percentile))
    mesh = ax.pcolormesh(kdist, energy, awk.T, cmap=cmap, shading="auto", vmin=vmin, vmax=vmax, **pcolorargs)
    lc = linecolor or s["ink2"]
    ax.set_xticks(kdist[kcrt])
    if klabel is not None:
        ax.set_xticklabels(list(klabel)[:len(kcrt)])
    else:
        ax.get_xaxis().set_visible(False)
    if show_kgrid:
        for x in kdist[kcrt]:
            ax.axvline(x, color=lc, lw=linewidth)
    if show_ef:
        ax.axhline(0.0, color=lc, lw=linewidth + 0.2, ls="--")
    return mesh


# ---------------------------------------------------------------- J_ij
def jij_limits(distance, jij, pad=0.05):
    """(xlim, ylim) with a 5 % margin, shared by the panels of one J_ij figure."""
    distance = np.asarray(distance, dtype=float)
    jij = np.asarray(jij, dtype=float)
    dx = pad * ((distance.max() - distance.min()) or 1.0)
    dy = pad * ((jij.max() - jij.min()) or 1.0)
    return (distance.min() - dx, distance.max() + dx), (jij.min() - dy, jij.max() + dy)


def plot_jij(ax, distance, jij, *, label=None, style=None, marker="o", xlim=None, ylim=None, a=1.0):
    """J_ij(R) of one pair on one Axes (R in units of a, or scaled by `a`); returns the Line2D list."""
    s = _style(style)
    distance = np.asarray(distance, dtype=float) * a
    jij = np.asarray(jij, dtype=float)
    order = np.argsort(distance)
    lines = ax.plot(distance[order], jij[order], color=_series(s, 0), lw=s["linewidth"] + 0.2, marker=marker,
                    markersize=4, label=label)
    ax.axhline(0.0, color=s["ink2"], lw=0.8, ls="--")
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    return lines


# ---------------------------------------------------------------- GAES
def shade_regions(ax, coarse=(), fine=(), style=None):
    """green bands for the coarse gap regions, blue for the fine sub-regions (GapRegion or (e1, e2))."""
    s = _style(style)
    for g in coarse:
        e1, e2 = (g.e1, g.e2) if hasattr(g, "e1") else g[:2]
        ax.axvspan(e1, e2, color=s["coarse"], alpha=s["coarse_alpha"], lw=0)
    for g in fine:
        e1, e2 = (g.e1, g.e2) if hasattr(g, "e1") else g[:2]
        ax.axvspan(e1, e2, color=s["fine"], alpha=s["fine_alpha"], lw=0)


def shade_ewidth_bounds(ax, bounds, e_min=None, style=None, label=True):
    """hatched band E_F - max_ewidth .. E_F - min_ewidth (None = open side: e_min / E_F)."""
    if bounds is None:
        return
    s = _style(style)
    lo, hi = bounds[0], bounds[1]
    if lo is None and hi is None:
        return
    left = -hi if hi is not None else (e_min if e_min is not None else ax.get_xlim()[0])
    right = -lo if lo is not None else 0.0
    ax.axvspan(left, right, facecolor="none", edgecolor=s["bounds"], hatch="///", lw=0, alpha=0.25,
               label="[min, max] ewidth = [%s, %s]" % (lo, hi) if label else None)


def draw_ewidth(ax, ewidth, label=None, style=None, **kw):
    """the contour bottom E_F - ewidth_go (dash-dot); ewidth may be None."""
    if ewidth is None:
        return
    s = _style(style)
    ax.axvline(-ewidth, **(dict(color=s["ewidth"], ls="-.", lw=1.6, label=label, zorder=6) | kw))


def draw_levels(ax, levels, highlight=(), e_min=None, y_text=None, fontsize=8, label_all=True, style=None):
    """core levels E - E_F from gaes.levels_from_go (key -> Level with .e / .star) or {key: [e, star]}:
    solid = core, dashed = '*' (valence). Keys in `highlight` use the level colour; others gray.
    y_text: data y of the labels (None: top of the axes in axes coordinates)."""
    s = _style(style)
    for key, v in (levels or {}).items():
        e, star = (v.e, v.star) if hasattr(v, "e") else (v[0], v[1])
        if e_min is not None and e <= e_min:
            continue
        hi = key in highlight
        color = s["level"] if hi else s["level_other"]
        ax.axvline(e, color=color, lw=1.5 if hi else 0.9, ls="--" if star else "-", alpha=1.0 if hi else 0.9)
        if hi or label_all:
            text = "%s%s %.2f" % (key, "*" if star else "", e)
            if y_text is None:
                ax.text(e, 0.98, text, rotation=90, fontsize=fontsize if hi else fontsize - 1, color=color,
                        ha="right", va="top", transform=ax.get_xaxis_transform())
            else:
                ax.text(e, y_text, text, rotation=90, fontsize=fontsize if hi else fontsize - 1, color=color,
                        ha="right", va="top")


def draw_thresholds(ax, dosth=None, dosth2=None, natm=None, style=None):
    """horizontal lines of the gap thresholds in the displayed (per cell) unit: dosth x natm (dashed)
    and dosth2 x natm (dotted). natm None = 1."""
    s = _style(style)
    scale = natm or 1
    if dosth is not None:
        ax.axhline(dosth * scale, color=s["ink2"], lw=0.6, ls="--")
    if dosth2 is not None:
        ax.axhline(dosth2 * scale, color=s["ink2"], lw=0.6, ls=":")


def plot_gaes_dos(ax, energy, dos, *, ewidth=None, ewidth_label=None, final=None, decision=None, coarse=(), fine=(),
                  bounds=None, levels=None, highlight=(), natm=None, dosth=None, dosth2=None, style=None,
                  log=True, ylim=(1e-4, 3e2), xlim=None, dos_label=None, efermi=True):
    """one GAES DOS panel: curve + gap regions + bounds band + contour bottom (+ final ewidth) + core
    levels + thresholds + E_F line.

    decision: pyakaikkr.gaes.Decision (its coarse / fine are used unless coarse / fine are given).
    bounds: (min_ewidth, max_ewidth) or None. levels: gaes.levels_from_go(...) dict or {key: [e, star]}.
    The DOS is drawn in the unit given (per cell); thresholds are scaled by natm.
    """
    s = _style(style)
    energy = np.asarray(energy, dtype=float)
    dos = np.asarray(dos, dtype=float)
    if dos.ndim == 2:
        dos = dos.sum(axis=0)
    lines = ax.plot(energy, dos, color=s["ink"] if style is None else _series(s, 0), lw=s["linewidth"], label=dos_label)
    if decision is not None and not coarse and not fine:
        coarse, fine = decision.coarse, decision.fine
    shade_regions(ax, coarse, fine, s)
    e_min = float(energy.min()) if energy.size else None
    shade_ewidth_bounds(ax, bounds, e_min=e_min, style=s, label=style is not None)
    draw_ewidth(ax, ewidth, label=ewidth_label, style=s)
    if final is not None and (ewidth is None or abs(final - ewidth) > 1e-6):
        ax.axvline(-final, color=s["final"], lw=1.0, ls="--", label="$-$ewidth final (%.4f)" % final)
    if efermi:
        mark_efermi(ax, s)
    draw_levels(ax, levels, highlight=highlight, e_min=e_min, y_text=None if (style is not None or ylim is None) else ylim[1] / 3, style=s)
    draw_thresholds(ax, dosth, dosth2, natm, s)
    if log:
        ax.set_yscale("log")
    if ylim is not None:
        ax.set_ylim(*ylim)
    if xlim is not None:
        ax.set_xlim(*xlim)
    return lines


# ---------------------------------------------------------------- figure export (PNG + SVG)
def _rasterize_meshes(fig):
    """rasterize colour meshes / images (A(w,k)) so that an SVG stays small: the mesh becomes an embedded
    PNG while axes, text and lines remain vector."""
    from matplotlib.collections import QuadMesh
    from matplotlib.image import AxesImage
    n = 0
    for ax in fig.axes:
        for art in list(ax.collections) + list(ax.images):
            if isinstance(art, (QuadMesh, AxesImage)):
                art.set_rasterized(True)
                n += 1
    return n


def figure_to_svg(fig, rasterize_meshes=True, dpi=150):
    """the figure as an SVG string (XML declaration and DOCTYPE removed, ready to inline in HTML)."""
    import io
    if rasterize_meshes:
        _rasterize_meshes(fig)
    buf = io.StringIO()
    fig.savefig(buf, format="svg", dpi=dpi, bbox_inches="tight", metadata={"Date": None})
    text = buf.getvalue()
    i = text.find("<svg")
    return text[i:] if i >= 0 else text


def figure_to_png(fig, dpi=150):
    """the figure as PNG bytes."""
    import io
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    return buf.getvalue()


def data_uri(png_bytes):
    """PNG bytes as a data URI for an <img> tag."""
    import base64
    return "data:image/png;base64," + base64.b64encode(png_bytes).decode("ascii")


def save_figure(fig, path_noext, formats=("png", "svg"), rasterize_meshes=True, dpi=150):
    """save the figure as <path_noext>.png and .svg (meshes rasterized inside the SVG); returns {fmt: path}."""
    out = {}
    for fmt in formats:
        path = "{}.{}".format(path_noext, fmt)
        if fmt == "svg":
            with open(path, "w") as f:
                f.write(figure_to_svg(fig, rasterize_meshes=rasterize_meshes, dpi=dpi))
        else:
            fig.savefig(path, dpi=dpi)
        out[fmt] = path
    return out


def gaes_legend_text():
    return ("black = DOS, green / blue = coarse / fine gap regions, hatched = [min, max] ewidth, "
            "red = -ewidth_go, blue lines = core levels (dashed = * valence), gray = other levels")


__all__ = ["DEFAULT_STYLE", "L_NAMES", "EWIDTH_LINE_LABEL", "mark_ewidth_go", "mark_efermi", "style_axes",
           "plot_dos", "plot_pdos", "component_names", "plot_awk", "jij_limits", "plot_jij",
           "shade_regions", "shade_ewidth_bounds", "draw_ewidth", "draw_levels", "draw_thresholds",
           "plot_gaes_dos", "gaes_legend_text", "figure_to_svg", "figure_to_png", "data_uri", "save_figure"]
