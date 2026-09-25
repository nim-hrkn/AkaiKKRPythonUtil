import numpy as np
import pytest

from pyakaikkr import GaesError
from pyakaikkr.gaes import gap_regions, check_ewidth, choose_ewidth, ewidth_candidates


def synthetic(gap=(-1.7, -0.7), floor=5e-4, semicore=-2.0):
    e = np.linspace(-2.25, 0.75, 201)
    d = np.full_like(e, 10.0)
    d[(e > gap[0]) & (e < gap[1])] = floor
    d[e < semicore] = 50.0
    return e, d


def test_single_gap_and_upper_edge_convention():
    e, d = synthetic()
    r = gap_regions(e, [d])
    assert len(r) == 1
    g = r[0]
    # e1: first low point; e2: the first point after the run (2019 convention)
    assert d[g.i1] < 1e-3 and d[g.i2] >= 1e-3 and d[g.i2 - 1] < 1e-3
    assert g.e1 == pytest.approx(e[g.i1]) and g.e2 == pytest.approx(e[g.i2])
    assert 0.9 < g.width < 1.1


def test_threshold_changes_width_and_no_derivative_needed():
    e, d = synthetic(floor=5e-4)
    # jagged DOS: single-point spikes inside the gap below the threshold and a one-point dip in the band
    d[80] = 9e-4
    d[150] = 5e-4   # one-point dip inside the valence band (would fool a derivative criterion)
    r1 = gap_regions(e, [d], dosth=1e-3)
    r2 = gap_regions(e, [d], dosth=6e-4)
    assert len(r1) == 2 and r1[0].width > 0.9   # the main gap plus the 1-point dip
    assert r1[1].width == pytest.approx(e[1] - e[0])
    # at 6e-4 the 9e-4 spike splits the main gap in two; the 1-point dip remains
    assert len(r2) == 3 and r2[0].i2 == 80 and r2[1].i1 == 81
    assert r2[0].width + r2[1].width < r1[0].width


def test_low_dos_to_the_end_of_mesh():
    e = np.linspace(-2.25, 0.75, 201)
    d = np.where(e < 0.0, 10.0, 1e-4)  # low DOS continues to the last point (2019 NameError case)
    r = gap_regions(e, [d])
    assert len(r) == 1 and r[0].i2 == len(e) - 1 and r[0].e2 == pytest.approx(e[-1])


def test_and_of_curves_and_mesh_mismatch():
    e, d1 = synthetic(gap=(-1.7, -0.7))
    _, d2 = synthetic(gap=(-1.4, -0.9))
    r = gap_regions(e, [d1, d2])
    assert len(r) == 1 and r[0].e1 == pytest.approx(-1.4, abs=0.02) and r[0].e2 == pytest.approx(-0.9, abs=0.02)
    with pytest.raises(GaesError):
        gap_regions(e, [d1[:-1]])
    with pytest.raises(GaesError):
        gap_regions(e, [])


def test_choose_ewidth_old_new_fail_and_min_ewidth():
    e, d = synthetic(gap=(-1.7, -0.7))
    r = gap_regions(e, [d])
    assert check_ewidth(r, 1.2) is not None
    assert choose_ewidth(r, 1.2)[0] == "old"
    flag, ew, cands = choose_ewidth(r, 0.5)          # -0.5 above the gap -> new
    assert flag == "new" and ew == pytest.approx(-(r[0].e2 - 0.2 - 0.01)) and cands == [ew]
    flag, ew, cands = choose_ewidth(r, 0.8)          # -0.8 inside but closer than ediff to the upper edge -> new
    assert flag == "new"
    assert choose_ewidth(r, None)[0] == "new"        # no current value is never "old"
    narrow = gap_regions(e, [synthetic(gap=(-1.2, -1.0))[1]])
    assert choose_ewidth(narrow, 1.1)[0] == "fail"   # width 0.2 < eth 0.3
    assert ewidth_candidates(r, min_ewidth=1.0) == []
    assert choose_ewidth(r, 0.5, min_ewidth=1.0)[0] == "fail"
