"""DOS/PDOS plots carry a vertical line at E-EF = -|ewidth| of the go run (not of the dos run)."""
import importlib
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from pyakaikkr import DosEXPlotter, DosPlotter, PDosEXPlotter
from pyakaikkr.DosPlotter import _EWIDTH_LINE_LABEL, _mark_ewidth_go, resolve_ewidth_go

CU_DIR = os.path.join(os.path.dirname(__file__), "..", "akaikkr", "Cu")
needs_cu = pytest.mark.skipif(
    not (os.path.isfile(os.path.join(CU_DIR, "out_go.log"))
         and os.path.isfile(os.path.join(CU_DIR, "out_dos.log"))),
    reason="tests/akaikkr/Cu/out_go.log and out_dos.log (made by testrun.py) are needed")


def _ewidth_lines(ax):
    return [l for l in ax.get_lines() if l.get_label() == _EWIDTH_LINE_LABEL]


def test_mark_line_inside_mesh():
    energy = np.linspace(-1.495, 0.495, 200)
    fig, ax = plt.subplots()
    assert _mark_ewidth_go(ax, energy, 1.0)
    lines = _ewidth_lines(ax)
    assert len(lines) == 1
    assert lines[0].get_xdata()[0] == pytest.approx(-1.0)
    plt.close(fig)


def test_mark_line_widens_xlim_when_below_mesh():
    energy = np.linspace(-0.995, 0.995, 200)  # akaikkr_cnd build: ref=0.5, ewidth_dos=2
    fig, ax = plt.subplots()
    ax.plot(energy, np.ones_like(energy))
    assert _mark_ewidth_go(ax, energy, 1.5)
    assert ax.get_xlim()[0] < -1.5
    assert _ewidth_lines(ax)[0].get_xdata()[0] == pytest.approx(-1.5)
    plt.close(fig)


def test_mark_line_absolute_value_and_none():
    energy = np.linspace(-1.5, 0.5, 10)
    fig, ax = plt.subplots()
    assert _mark_ewidth_go(ax, energy, -1.2)
    assert _ewidth_lines(ax)[0].get_xdata()[0] == pytest.approx(-1.2)
    assert not _mark_ewidth_go(ax, energy, None)
    assert len(_ewidth_lines(ax)) == 1
    plt.close(fig)


def test_resolve_priority_argument_first(tmp_path):
    # explicit ewidth_go wins even if out_go.log is missing, and no warning is issued
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert resolve_ewidth_go(str(tmp_path), ewidth_go=1.3) == 1.3


def test_resolve_missing_go_outfile_warns_and_returns_none(tmp_path):
    with pytest.warns(UserWarning, match="not found"):
        assert resolve_ewidth_go(str(tmp_path)) is None
    # explicitly switched off: no line, no warning
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert resolve_ewidth_go(str(tmp_path), read_go_outfile=False) is None


def test_dosplotter_make_with_ewidth_go(tmp_path):
    energy = list(np.linspace(-1.495, 0.495, 200))
    dos = list(np.abs(np.sin(np.linspace(0, 6, 200))) + 1e-3)
    path = DosPlotter(str(tmp_path)).make(energy, [dos, dos], ewidth_go=1.0)
    assert os.path.isfile(path)


@needs_cu
def test_resolve_reads_go_not_dos():
    # Cu test case: go ewidth=1.0, dos ewidth=2.0
    assert resolve_ewidth_go(CU_DIR) == pytest.approx(1.0)
    assert resolve_ewidth_go(CU_DIR, go_outfile="out_dos.log") == pytest.approx(2.0)


@needs_cu
def test_explotters_use_go_ewidth(tmp_path, monkeypatch):
    drawn = []

    def spy(ax, energy, ewidth_go):
        drawn.append(ewidth_go)
        return _mark_ewidth_go(ax, energy, ewidth_go)

    # pyakaikkr.DosPlotter as attribute is the class; patch the module object
    monkeypatch.setattr(importlib.import_module("pyakaikkr.DosPlotter"), "_mark_ewidth_go", spy)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        DosEXPlotter(CU_DIR, "out_dos.log", str(tmp_path)).make()
        PDosEXPlotter(CU_DIR, "out_dos.log", str(tmp_path)).make()
    assert drawn and all(e == pytest.approx(1.0) for e in drawn)
    assert os.path.isfile(os.path.join(str(tmp_path), "dos.png"))
    assert os.path.isfile(os.path.join(str(tmp_path), "pdos_0.png"))

    # argument overrides the file
    drawn.clear()
    DosEXPlotter(CU_DIR, "out_dos.log", str(tmp_path)).make(ewidth_go=0.8)
    assert drawn == [pytest.approx(0.8)]

    # no go file available -> warning, no line
    drawn.clear()
    with pytest.warns(UserWarning):
        DosEXPlotter(CU_DIR, "out_dos.log", str(tmp_path), go_outfile="no_such.log").make()
    assert drawn == [None]
