"""per-orbital valence / core rules -> ewidth bounds (docs/ewidth_tuning_scheme.md section 15)."""
import os

import numpy as np
import pytest

from pyakaikkr import AkaikkrJob
from pyakaikkr.Error import GaesError
from pyakaikkr.gaes import (parse_orbital_rules, levels_from_go, levels_from_tables, bounds_from_rules, check_rules,
                            initial_ewidth, choose_ewidth, choose_ewidth2, decide, Gaes)
from pyakaikkr.gaes.gap import GapRegion

from kkr_env import DATA_DIR as DATA
RB_GO = os.path.join(DATA, "out_go_RbMnFeCo_fcc_ew1.2.log")   # Rb-Mn-Fe-Co fcc, ewidth 1.2, E_F 0.245 / 0.280


def test_parse_rules_and_aliases():
    rules = parse_orbital_rules(["Rb4p=occupied", "Se 4s:core", "Bi6s valence", "In4d=UNOCCUPIED"])
    assert [str(r) for r in rules] == ["Rb4p=valence", "Se4s=core", "Bi6s=valence", "In4d=core"]
    assert parse_orbital_rules(["Rb4p=core", "Rb4p=core"]) == parse_orbital_rules(["Rb4p=core"])
    with pytest.raises(GaesError):
        parse_orbital_rules(["Rb4p=core", "Rb4p=valence"])
    with pytest.raises(GaesError):
        parse_orbital_rules(["Rb4p=filled"])
    with pytest.raises(GaesError):
        parse_orbital_rules(["rb4p=core"])


def test_core_levels_by_component_and_levels_from_go():
    recs = AkaikkrJob(DATA).get_core_levels_by_component(os.path.basename(RB_GO))
    rb = [r for r in recs if r["element"] == "Rb" and r["orbital"] == "4p"]
    assert {r["spin"] for r in rb} == {"up", "down"} and all(r["star"] for r in rb)
    assert rb[0]["level_Ry"] == pytest.approx(-0.6911086) and rb[0]["ef_Ry"] == pytest.approx(0.2448142)
    assert rb[0]["e_minus_ef_Ry"] == pytest.approx(-0.9359, abs=1e-3)
    assert not [r for r in recs if r["element"] == "Rb" and r["orbital"] == "4s" and r["star"]]
    lv = levels_from_go(RB_GO)
    assert lv["Rb4p"].star is True and lv["Rb4p"].e == pytest.approx(-0.937, abs=2e-3)
    assert lv["Rb4p"].e_min <= lv["Rb4p"].e <= lv["Rb4p"].e_max
    assert lv["Rb4s"].star is False and lv["Rb4s"].e == pytest.approx(-2.10, abs=1e-2)
    assert "Mn3d" not in lv     # not core-configured: never in the core level table


def test_bounds_from_rules_on_the_rb_output():
    lv = levels_from_go(RB_GO)
    b = bounds_from_rules(parse_orbital_rules(["Rb4p=valence"]), lv, ediff=0.2)
    assert b.min_ewidth == pytest.approx(0.938 + 0.2, abs=2e-3) and b.max_ewidth is None
    b = bounds_from_rules(parse_orbital_rules(["Rb4p=core"]), lv, ediff=0.2)
    assert b.min_ewidth is None and b.max_ewidth == pytest.approx(0.936 - 0.2, abs=2e-3)
    b = bounds_from_rules(parse_orbital_rules(["Rb4p=valence", "Rb4s=core"]), lv, ediff=0.2)
    assert b.min_ewidth == pytest.approx(1.138, abs=2e-3) and b.max_ewidth == pytest.approx(1.899, abs=2e-3)
    # explicit user bounds are intersected; defaults are replaced
    b = bounds_from_rules(parse_orbital_rules(["Rb4p=core"]), lv, 0.2, user_min=0.5, default_min=1.0, default_max=2.0)
    assert b.as_list() == pytest.approx([0.5, 0.736], abs=2e-3)
    b = bounds_from_rules([], lv, 0.2, default_min=1.0, default_max=2.0)
    assert b.as_list() == [1.0, 2.0]
    with pytest.raises(GaesError):     # Rb 4p core needs ewidth <= 0.74, Rb 4s valence needs >= 2.3
        bounds_from_rules(parse_orbital_rules(["Rb4p=core", "Rb4s=valence"]), lv, 0.2)
    with pytest.raises(GaesError):     # Mn 3d is not core-configured: cannot be made core
        bounds_from_rules(parse_orbital_rules(["Mn3d=core"]), lv, 0.2)
    # a rule that gives no bound leaves the defaults in force
    assert bounds_from_rules(parse_orbital_rules(["Mn3d=valence"]), lv, 0.2, default_min=1.0).min_ewidth == 1.0


def test_check_rules_reports_mismatch():
    lv = levels_from_go(RB_GO)
    assert check_rules(parse_orbital_rules(["Rb4p=valence", "Rb4s=core"]), lv) == []
    bad = check_rules(parse_orbital_rules(["Rb4p=core"]), lv)
    assert len(bad) == 1 and bad[0]["rule"] == "Rb4p=core" and bad[0]["star"] is True


def test_levels_from_tables_prefers_the_converged_table():
    lv = levels_from_tables(["Bi", "Rb"], ef_assumed=0.6)
    assert lv["Bi6s"].source == "table_converged" and lv["Bi6s"].e == pytest.approx(-0.882)
    assert lv["Bi6s"].e_min == pytest.approx(-0.945) and lv["Bi6s"].e_max == pytest.approx(-0.8)
    assert lv["Rb4p"].source == "table_atomic" and lv["Rb4p"].e == pytest.approx(-0.862 - 0.6, abs=1e-3)


def test_initial_ewidth_is_only_raised():
    b = bounds_from_rules(parse_orbital_rules(["Rb4p=core"]), levels_from_go(RB_GO), 0.2)
    assert initial_ewidth(1.2, b) == 1.2            # a core rule never lowers the first ewidth (the first go's DOS does)
    assert initial_ewidth(0.6, b) == 0.6
    b = bounds_from_rules(parse_orbital_rules(["Rb4p=valence", "Rb4s=core"]), levels_from_go(RB_GO), 0.2)
    assert initial_ewidth(0.8, b) == pytest.approx(1.138 + 0.2, abs=2e-3)   # raised to min + ediff
    assert initial_ewidth(1.5, b) == 1.5
    b = bounds_from_rules(parse_orbital_rules(["Rb4s=valence"]), levels_from_go(RB_GO), 0.2)
    assert initial_ewidth(1.2, b) == pytest.approx(2.30 + 0.2, abs=2e-3)


def test_old_is_refused_outside_the_bounds():
    """the defect found while writing section 15: an ewidth anchored in a gap was 'old' even when it
    violated min_ewidth / max_ewidth, so an orbital rule could never move it."""
    r = [GapRegion(-1.7, -0.7, 10, 40)]
    assert choose_ewidth(r, 1.2)[0] == "old"
    flag, ew, _ = choose_ewidth(r, 1.2, max_ewidth=1.0)
    assert flag == "new" and ew == pytest.approx(0.91)                   # -(e2 - ediff - margin), inside the bound
    assert choose_ewidth(r, 1.2, min_ewidth=1.5)[:2] == ("new", 1.5)
    assert choose_ewidth(r, 1.2, min_ewidth=1.8)[0] == "fail"           # -1.8 outside the gap
    e = np.linspace(-2.2425, 0.7425, 200)
    d = np.full_like(e, 5e-4)
    d[e > -0.50] = 5.0
    d[(e > -1.70) & (e < -1.30)] = 5.0
    assert choose_ewidth2(e, [d], 1.2, dosth=2e-2).flag == "old"
    m = choose_ewidth2(e, [d], 1.2, dosth=2e-2, max_ewidth=0.9)
    assert m.flag == "new" and m.ewidth == pytest.approx(0.70, abs=0.02)     # top of the gap [-1.3, -0.5] minus ediff
    m = choose_ewidth2(e, [d], 1.2, dosth=2e-2, min_ewidth=1.9)
    assert m.flag == "new" and m.ewidth == pytest.approx(1.91, abs=0.02)     # deep region below the X band
    assert decide(1, e, [d], 1.2, dosth=2e-2, max_ewidth=0.9).flag == "new"


def test_gaes_keeps_user_bounds_and_records_orbitals():
    g = Gaes("specx")
    assert (g.min_ewidth, g.max_ewidth, g._user_min, g._user_max) == (1.0, 2.0, None, None)
    g = Gaes("specx", min_ewidth=0.5, orbitals=["Bi6s=unoccupied"])
    assert g.parameters()["orbitals"] == ["Bi6s=core"] and g._user_min == 0.5
    b = g._bounds(levels_from_tables(["Bi"]))
    assert b.min_ewidth == 0.5 and b.max_ewidth == pytest.approx(0.8 - 0.2)      # shallowest table value


def test_step0_uses_only_the_converged_table():
    from pyakaikkr.gaes import levels_for_step0
    lv = levels_for_step0(["Bi", "Rb"])
    assert "Bi6s" in lv and "Bi5d" in lv and "Rb4p" not in lv      # Rb has no converged entry: no step-0 bound
    g = Gaes("specx", orbitals=["Rb4p=valence"])
    assert g._bounds(levels_for_step0(["Rb"]), strict=False).as_list() == [1.0, 2.0]     # no level -> the defaults stay
    g = Gaes("specx", orbitals=["Rb4p=core"])
    assert g._bounds(levels_for_step0(["Rb"]), strict=False).as_list() == [1.0, 2.0]     # step 0: missing level is no error
    with pytest.raises(GaesError):
        g._bounds(levels_for_step0(["Rb"]))                                            # after a go it is
