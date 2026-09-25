"""Method 2 (two thresholds) of the ewidth decision."""
import glob
import os

import numpy as np
import pytest

from pyakaikkr import AkaikkrJob
from pyakaikkr.gaes import choose_ewidth2, decide, dos_curves_from_outputs
from gaes_env import RUN_DIR


def synthetic(floor=8e-4, tail=3e-3):
    """valence band above -1.0, semicore peak at -1.75 with a tail, gap floor `floor`."""
    e = np.linspace(-2.2425, 0.7425, 200)
    d = np.full_like(e, floor)
    d[e > -1.05] = 5.0
    d += 150.0 * np.exp(-((e + 1.75) / 0.03) ** 2) + tail * np.exp(-((e + 1.75) / 0.25) ** 2)
    d[e < -2.0] = 20.0
    return e, d


def test_method2_rejects_semicore_tail_that_method1_accepts():
    e, d = synthetic()
    m1 = decide(1, e, [d], 1.6, dosth=2e-2)
    m2 = decide(2, e, [d], 1.6, dosth=2e-2, dosth2=1e-3)
    assert m1.flag == "old"                     # tail at -1.6 is below 2e-2
    assert m2.flag == "new" and 1.1 < m2.ewidth < 1.35
    assert m2.coarse and m2.fine and not m2.relaxed
    f = m2.fine[0]
    assert f.e1 > -1.6                          # -1.6 is outside the fine sub-region
    top = m2.coarse[0].e2 - 0.20
    assert -m2.ewidth == pytest.approx(min(f.e2, top) - 0.01)
    # a value inside the fine region and below top is "old"
    assert decide(2, e, [d], m2.ewidth, dosth=2e-2).flag == "old"


def test_method2_relaxes_dosth2_when_floor_is_above_it():
    e, d = synthetic(floor=1.2e-3)
    m2 = choose_ewidth2(e, [d], 1.6, dosth=2e-2, dosth2=1e-3)
    assert m2.flag == "new" and m2.relaxed and m2.dosth2_used == pytest.approx(2e-3)
    m2n = choose_ewidth2(e, [d], 1.6, dosth=2e-2, dosth2=1e-3, dosth2_relax=1.0)
    assert m2n.flag == "fail"


def test_method2_bottom_region_needs_verified_low_dos_below():
    e, d = synthetic()
    d[e < -2.0] = 1e-4                          # low DOS continues to the mesh bottom (-2.2425)
    m2 = choose_ewidth2(e, [d], 0.5, dosth=2e-2)
    assert m2.coarse[-1].i1 == 0
    # the real gap gives the first candidate; the bottom region [-2.24, -2.0] has 0.24 < eth of verified
    # low DOS below any position in it, so it yields nothing and is flagged
    assert 1.1 < m2.candidates[0] < 1.35 and m2.window_limited
    assert choose_ewidth2(e, [d], 2.1, dosth=2e-2).flag == "new"        # -2.1 in the bottom region is not "old"
    only_bottom = np.where(e < -1.9, 1e-4, 5.0)
    m = choose_ewidth2(e, [only_bottom], 2.1, dosth=2e-2)
    assert m.flag == "fail" and m.window_limited
    # a bottom region with enough verified low DOS below is usable: gap from the bottom up to -1.4
    wide = np.where(e < -1.4, 1e-4, 5.0)
    m = choose_ewidth2(e, [wide], 1.2, dosth=2e-2, max_ewidth=2.0)
    assert m.flag == "new" and m.ewidth == pytest.approx(1.61, abs=0.02) and not m.window_limited
    assert choose_ewidth2(e, [wide], 1.7, dosth=2e-2, max_ewidth=2.0).flag == "old"
    assert choose_ewidth2(e, [wide], 2.0, dosth=2e-2, max_ewidth=2.0).flag == "new"   # -2.0 has only 0.24 Ry verified below
    assert choose_ewidth2(e, [np.full_like(e, 5.0)], 1.2).flag == "fail"


@pytest.mark.skipif(not glob.glob(os.path.join(RUN_DIR, "*Bi*_fcc", "key_*,ew_000-1.6000,*", "out_dos.log")),
                    reason="Bi test runs (ewidth 1.6) not present")
def test_method2_on_bi_runs():
    expected = {"Al0p25Si0p25Rh0p25Bi0p25_fcc": 1.2825, "Al0p25Sc0p25Ni0p25Bi0p25_fcc": 1.2525,
                "Al0p25Si0p25Ge0p25Bi0p25_fcc": 1.2375}
    for name, ew in expected.items():
        dirs = glob.glob(os.path.join(RUN_DIR, name, "key_*,ew_000-1.6000,ed_000-1e-04,*,pm_000-5e-03"))
        if not dirs:
            continue
        e, c, _ = dos_curves_from_outputs([(AkaikkrJob(dirs[0]), "out_dos.log")])
        m2 = decide(2, e, c, 1.6, dosth=2e-2, dosth2=1e-3)
        assert m2.flag == "new" and m2.ewidth == pytest.approx(ew, abs=1e-3), name
        assert decide(1, e, c, 1.6, dosth=2e-2).flag == "old"


def test_candidate_is_clipped_to_the_bounds_inside_the_gap():
    """Pm-Mn-Fe-Co case: fine gap [-1.18, -0.77] straddles -1.0; the shallowest candidate 0.78 is below
    min_ewidth=1.0, but ewidth=1.0 still lies inside the gap, so it must be proposed (not 'fail')."""
    e = np.linspace(-2.2425, 0.7425, 200)
    d = np.full_like(e, 8e-4)
    d[e > -0.50] = 5.0                              # valence band bottom at -0.50
    d[(e > -1.70) & (e < -1.30)] = 5.0              # X 5p band
    d[(e > -1.30) & (e < -1.18)] = 2e-3             # its upper tail (above dosth2, below dosth)
    d[(e > -0.77) & (e < -0.50)] = 1.5e-3           # valence tail
    m = choose_ewidth2(e, [d], 1.2, dosth=2e-2, dosth2=1e-3, min_ewidth=1.0, max_ewidth=2.0)
    assert m.flag == "new" and m.ewidth == pytest.approx(1.0)
    m2 = choose_ewidth2(e, [d], 1.2, dosth=2e-2, dosth2=1e-3, min_ewidth=0.5, max_ewidth=2.0)
    assert m2.flag == "new" and m2.ewidth == pytest.approx(0.7825, abs=0.02)
    # a bound outside the upper gap: the upper gap gives nothing, the deep region below the X band
    # (touching the window bottom, but with > eth of verified low DOS) gives the candidate instead
    m3 = choose_ewidth2(e, [d], 1.2, dosth=2e-2, dosth2=1e-3, min_ewidth=1.25, max_ewidth=2.0)
    assert m3.flag == "new" and m3.ewidth == pytest.approx(1.91, abs=0.02)
    # and ewidth = 1.0 is then judged "old"
    assert choose_ewidth2(e, [d], 1.0, dosth=2e-2, dosth2=1e-3, min_ewidth=1.0, max_ewidth=2.0).flag == "old"
    # Method 1: same clipping through ewidth_candidates
    from pyakaikkr.gaes import decide
    m1 = decide(1, e, [d], 1.2, dosth=2e-2, min_ewidth=1.0, max_ewidth=2.0)
    assert m1.flag == "old"   # at 2e-2 the coarse gap [-1.3, -0.5] contains -1.2
