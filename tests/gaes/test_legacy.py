"""checks against the 2019 RUN data (skipped when the directory is absent)."""
import os

import pytest

from pyakaikkr import AkaikkrJob
from pyakaikkr.gaes import gap_regions, choose_ewidth, dos_curves_from_outputs, KkrRunner, Layout, collect_legacy
from kkr_env import legacy_dir, needs_run0, RUN0_RUN


@needs_run0
def test_gap_regions_reproduce_str_out_13142122():
    pairs = [(AkaikkrJob(legacy_dir("13142122", p)), "out_dos.log") for p in ("bcc", "fcc")]
    energy, curves, _ = dos_curves_from_outputs(pairs)
    assert energy.min() == pytest.approx(-2.2425) and energy.max() == pytest.approx(0.7425)
    regions = gap_regions(energy, curves, dosth=2e-2)     # 2019: 1e-2 on the spin average
    assert [g.as_tuple() for g in regions] == [pytest.approx((-2.2425, -2.1075)), pytest.approx((-2.0025, -0.7275))]
    flag, ew, cands = choose_ewidth(regions, 1.2)
    assert flag == "old" and ew == 1.2 and cands == [pytest.approx(0.9375)]


@needs_run0
def test_gap_regions_reproduce_str_out_13487580():
    pairs = [(AkaikkrJob(legacy_dir("13487580", p)), "out_dos.log") for p in ("bcc", "fcc")]
    energy, curves, _ = dos_curves_from_outputs(pairs)
    regions = gap_regions(energy, curves, dosth=2e-2)
    assert [g.as_tuple() for g in regions] == [pytest.approx((-2.2425, -0.7425))]
    assert choose_ewidth(regions, 1.2)[2] == [pytest.approx(0.9525)]


@needs_run0
def test_runner_result_of_legacy_directory():
    d = legacy_dir("13487580", "fcc")
    r = KkrRunner("", str(d), {"go": "go"}, compat=True)
    row = r.result(dosth=2e-2, save_csv=False)
    assert row["converged"] is True
    assert row["total_energy_Ry"] == pytest.approx(-21075.683771993)
    assert row["a_bohr"] == pytest.approx(8.01512)
    assert row["Tc_K"] == pytest.approx(0.0)
    assert all(abs(c["spin"]) < 1e-6 for c in row["component_moment"])
    assert row["ewidth"] == pytest.approx(1.2) and row["ewidth_dos"] == pytest.approx(3.0)
    assert row["low_dos_regions"] == [pytest.approx((-2.2425, -0.7425))]


@needs_run0
def test_layout_v1_finds_legacy_dirs():
    L = Layout(RUN0_RUN, version=1)
    pts = L.find(key="13487580")
    assert {p.polytyp for p in pts} == {"bcc", "fcc"}
