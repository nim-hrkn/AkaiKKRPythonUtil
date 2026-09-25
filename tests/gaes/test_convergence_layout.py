import numpy as np
import pytest

from pyakaikkr import GaesError
from pyakaikkr.gaes import is_converging, Layout, RunPoint


def test_is_converging():
    err = list(np.linspace(-1.0, -5.0, 150))
    mom = [2.0] * 150
    assert is_converging(err, mom)
    rng = np.random.default_rng(0)
    assert not is_converging(list(rng.normal(-3, 0.5, 150)), list(rng.normal(2, 0.5, 150)))
    assert is_converging(list(np.linspace(-3.0, -3.001, 150)), list(rng.normal(2, 1, 150)))  # mae_err tiny
    assert not is_converging([-1.0, -2.0], [1.0, 1.0])  # fewer than 3 points


def test_layout_round_trip_and_legacy_parse():
    p = RunPoint("13142122", 1, 1.2, 0, 1e-4, "fcc", 2, 5e-3)
    L2 = Layout("RUN", version=2)
    name = L2.dirname(p)
    assert name == "key_13142122,ew_001-1.2000,ed_000-1e-04,polytyp_fcc,pm_002-5e-03"
    assert L2.parse(name) == p
    L1 = Layout("RUN", version=1)
    assert L1.dirname(p) == "key_13142122,ew_001,ed_000,polytyp_fcc,pm_002"
    q = L1.parse("key_13142122,ew_000,ed_000,polytyp_bcc,pm_000")
    assert (q.key, q.iew, q.ied, q.polytyp, q.ipm) == ("13142122", 0, 0, "bcc", 0)
    assert q.ewidth is None
    with pytest.raises(GaesError):
        L1.parse("nonsense")
    assert Layout("RUN", version=2, sep=";").dirname(p).count(";") == 4
    with pytest.raises(ValueError):
        Layout("RUN", version=3)


def test_layout_find(tmp_path):
    L = Layout(str(tmp_path), version=2)
    for p in (RunPoint("a", 0, 1.2, 0, 1e-4, "fcc", 0, 5e-3), RunPoint("a", 1, 0.9, 0, 1e-4, "bcc", 0, 5e-3)):
        (tmp_path / L.dirname(p)).mkdir()
    (tmp_path / "other").mkdir()
    assert len(L.find()) == 2
    assert [p.polytyp for p in L.find(iew=1)] == ["bcc"]
