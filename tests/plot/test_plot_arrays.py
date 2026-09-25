"""pyakaikkr.plot: array-based drawing shared with aiida-akaikkr (no specx, no files)."""
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from pyakaikkr.plot import (plot_dos, plot_pdos, plot_awk, plot_jij, jij_limits, plot_gaes_dos, mark_ewidth_go,
                            component_names, DEFAULT_STYLE)  # noqa: E402
from pyakaikkr.gaes.plot import draw_gaes_dos, legend_text  # noqa: E402


def _ax():
    fig, ax = plt.subplots()
    return fig, ax


def test_plot_dos_two_spins_mirrored_and_ewidth_line():
    e = np.linspace(-2.0, 1.0, 60)
    d = np.vstack([np.ones(60), 2 * np.ones(60)])
    fig, ax = _ax()
    lines = plot_dos(ax, e, d, ewidth_go=1.2)
    assert len(lines) == 2 and lines[1].get_ydata().max() < 0          # down spin mirrored
    xs = [l.get_xdata()[0] for l in ax.lines if len(set(l.get_xdata())) == 1]
    assert any(abs(x + 1.2) < 1e-9 for x in xs) and any(abs(x) < 1e-9 for x in xs)   # -ewidth and E_F lines
    fig, ax = _ax()
    lines = plot_dos(ax, e, d, yscale="log")                          # log: no mirroring
    assert lines[1].get_ydata().min() > 0 and ax.get_yscale() == "log"
    fig, ax = _ax()
    assert len(plot_dos(ax, e, d[0])) == 1                              # 1-D input = one spin
    plt.close("all")


def test_ewidth_line_widens_the_x_range_and_none_draws_nothing():
    e = np.linspace(-1.0, 0.5, 20)
    fig, ax = _ax()
    assert mark_ewidth_go(ax, e, 1.5) and ax.get_xlim()[0] < -1.5
    fig, ax = _ax()
    assert not mark_ewidth_go(ax, e, None) and not ax.lines
    plt.close("all")


def test_plot_pdos_skips_nan_padding_and_one_spin():
    e = np.linspace(-2.0, 1.0, 30)
    p = np.abs(np.random.default_rng(0).random((2, 30, 4)))
    p[:, :, 3] = np.nan                                                 # f column padded (mxl 2)
    fig, ax = _ax()
    lines = plot_pdos(ax, e, p)
    assert [l.get_label() for l in lines] == ["s", "p", "d"]
    fig, ax = _ax()
    assert len(plot_pdos(ax, e, p, spin=1, nl=2)) == 2 and all(l.get_ydata().min() >= 0 for l in ax.lines[:2])
    plt.close("all")


def test_component_names_from_type_of_site():
    tos = [{"type": "HEA", "comp_shortname": ["HEA_Rh_50.0%", "HEA_Pt_50.0%"]}, {"type": "Cu", "comp_shortname": ["Cu"]}]
    assert component_names(tos) == ["Rh 50.0% in HEA", "Pt 50.0% in HEA", "Cu"]
    assert component_names(tos, long=False) == ["Rh", "Pt", "Cu"]


def test_plot_awk_mesh_and_ticks():
    kdist = np.linspace(0, 3, 12)
    energy = np.linspace(-1, 1, 20)
    awk = np.abs(np.random.default_rng(1).random((12, 20)))
    fig, ax = _ax()
    mesh = plot_awk(ax, kdist, energy, awk, kcrt=[0, 5, 11], klabel=["G", "X", "L"])
    assert mesh is not None and [t.get_text() for t in ax.get_xticklabels()] == ["G", "X", "L"]
    plt.close("all")


def test_plot_jij_sorted_with_shared_limits():
    dist = np.array([1.0, 0.5, 1.5])
    jij = np.array([2.0, 5.0, -1.0])
    xlim, ylim = jij_limits(dist, jij)
    fig, ax = _ax()
    lines = plot_jij(ax, dist, jij, label="Fe-Fe", xlim=xlim, ylim=ylim)
    assert list(lines[0].get_xdata()) == [0.5, 1.0, 1.5] and ax.get_xlim() == xlim
    plt.close("all")


def test_plot_gaes_dos_elements_and_alias():
    e = np.linspace(-2.2, 0.7, 100)
    d = np.full(100, 5e-4)
    d[e > -0.5] = 5.0
    fig, ax = _ax()
    plot_gaes_dos(ax, e, [d, d], ewidth=1.2, final=1.4, coarse=[(-2.2, -0.5)], fine=[(-2.2, -0.7)],
                  bounds=(1.0, 2.0), levels={"Bi6s": [-0.9, True], "Bi5d": [-1.7, False]}, highlight=["Bi6s"],
                  natm=2, dosth=2e-2, dosth2=1e-3, style={"series": ["#2a78d6"]})
    assert ax.get_yscale() == "log"
    xs = [l.get_xdata()[0] for l in ax.lines if len(set(l.get_xdata())) == 1]
    for x in (-1.2, -1.4, -0.9, -1.7, 0.0):
        assert any(abs(v - x) < 1e-9 for v in xs), x
    ys = [l.get_ydata()[0] for l in ax.lines if len(set(l.get_ydata())) == 1 and len(set(l.get_xdata())) > 1]
    assert any(abs(y - 4e-2) < 1e-12 for y in ys) and any(abs(y - 2e-3) < 1e-12 for y in ys)   # thresholds x natm
    assert len(ax.patches) >= 3                                        # coarse, fine, bounds bands
    assert draw_gaes_dos is plot_gaes_dos and "coarse" in legend_text()
    plt.close("all")


def test_style_override_does_not_change_default():
    before = dict(DEFAULT_STYLE)
    e = np.linspace(-1, 1, 10)
    fig, ax = _ax()
    lines = plot_dos(ax, e, np.ones(10), style={"series": ["#123456"]})
    assert lines[0].get_color() == "#123456" and DEFAULT_STYLE == before
    plt.close("all")
