import io
import os

import pytest

from pyakaikkr import AkaikkrJob, GaesError
from pyakaikkr.gaes import (SiteComposition, make_single_site_param, check_type_name, TYPE_NAME_MAX_LEN,
                            heakey_to_composition, composition_to_heakey)
from kkr_env import legacy_dir, needs_run0


def test_from_type_name_and_back():
    c = SiteComposition.from_type_name("Rh0.5Pt0.5_1")
    assert c.elements == ("Rh", "Pt") and c.fractions == (0.5, 0.5) and c.z == (45, 78) and c.conc == (50, 50)
    assert c.type_name(suffix="_1") == "Rh0.5Pt0.5_1"
    b = SiteComposition.from_type_name("B0.975Vc0.025")
    assert b.z == (5, 0) and b.conc == (97.5, 2.5) and b.type_name() == "B0.975Vc0.025"
    assert SiteComposition.from_type_name("Fe").type_name() == "Fe"
    assert SiteComposition.from_type_name("Fe_0").conc == (100,)
    e = SiteComposition.from_elements("AlSiScTi")
    assert e.is_equiatomic and e.conc == (25, 25, 25, 25) and e.type_name() == "Al0.25Si0.25Sc0.25Ti0.25"
    assert SiteComposition.from_elements(["Al", "Si", "Sc", "Ti", "V"]).conc[0] == 20
    assert SiteComposition.from_dict({"La": 0.999, "Ge": 0.001}).type_name() == "La0.999Ge0.001"
    assert e.key() == "Al0p25Si0p25Sc0p25Ti0p25"


def test_type_name_limits():
    assert TYPE_NAME_MAX_LEN == 40
    check_type_name("x" * 40)
    with pytest.raises(GaesError):
        check_type_name("x" * 41)
    with pytest.raises(GaesError):
        check_type_name("a b")
    with pytest.raises(GaesError):
        check_type_name("a,b")
    eight = SiteComposition.from_elements("AlSiScTiVCrMnFe")
    with pytest.raises(GaesError):
        eight.type_name(suffix="_1")            # 57 characters
    assert len(eight.type_name(ndigits=2, max_len=10 ** 6)) > 40   # still too long with 2 digits
    five = SiteComposition.from_elements("AlSiScTiV")
    assert five.type_name(suffix="_1") == "Al0.2Si0.2Sc0.2Ti0.2V0.2_1"
    with pytest.raises(GaesError):
        SiteComposition.from_type_name("Rh0.5Pt")   # mixed with/without fractions
    with pytest.raises(GaesError):
        SiteComposition.from_type_name("Xx0.5Pt0.5")
    with pytest.raises(GaesError):
        SiteComposition(("Rh", "Pt"), (0.6, 0.6))


def test_heakey():
    c = heakey_to_composition("13142122")
    assert c.elements == ("Al", "Si", "Sc", "Ti")
    assert composition_to_heakey(c) == "13142122"
    assert heakey_to_composition(1314212223).elements == ("Al", "Si", "Sc", "Ti", "V")
    with pytest.raises(GaesError):
        heakey_to_composition("1314212")
    with pytest.raises(GaesError):
        composition_to_heakey(SiteComposition.from_type_name("Rh0.4Pt0.6"))


def _tokens(text):
    out = []
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        row = []
        for t in line.split():
            t = t.rstrip("abc") if t[:1].isdigit() and t[-1:] in "abc" else t
            try:
                row.append(float(t))
            except ValueError:
                row.append(t)
        out.append(row)
    return out


@needs_run0
def test_single_site_param_matches_2019_inputcard():
    d = legacy_dir("13487580", "fcc")
    ref = open(os.path.join(d, "inputcard_go")).read()
    comp = heakey_to_composition("13487580")
    dic = make_single_site_param(comp, "fcc", type_name="HEA")
    buf = io.StringIO()
    AkaikkrJob(str(d)).make_inputcard(dic, buf)
    assert _tokens(buf.getvalue()) == _tokens(ref)


def test_single_site_param_lattice_and_option():
    comp = SiteComposition.from_type_name("Rh0.5Pt0.5")
    assert make_single_site_param(comp, "fcc")["a"] == 1000000
    assert make_single_site_param(comp, "fcc", lattice="mjw")["a"] == 0
    assert make_single_site_param(comp, "fcc", lattice=7.2)["a"] == 7.2
    d = make_single_site_param(comp, "bcc", option={"number_emesh": 5})
    assert d["option"] == {"mse": "5"} and d["type"] == ["Rh0.5Pt0.5"] and d["conc"] == [[50, 50]]
    with pytest.raises(GaesError):
        make_single_site_param(comp, "fcc", lattice=-1)
