"""tests of the lattice mapping and the calculator input (no specx needed)."""
import io
import os

import numpy as np
import pytest
from ase.build import bulk
from ase.calculators.calculator import InputError, all_changes

from pyakaikkr import AkaikkrJob
from pyakaikkr.ase import AkaiKKR, atoms_to_kkr_lattice, atoms_to_kkr_param, set_occupancy

BOHR = 0.529177


def test_lattice_aux():
    cu = bulk("Cu", "fcc", a=3.615)
    lat = atoms_to_kkr_lattice(cu)
    assert lat["brvtyp"] == "aux"
    a = np.linalg.norm(cu.cell[0])
    assert lat["a"] == pytest.approx(a / BOHR)
    r = np.array([lat["r1"], lat["r2"], lat["r3"]])
    assert np.allclose(r * a, np.array(cu.cell))
    assert np.linalg.norm(r[0]) == pytest.approx(1.0)


def test_left_handed_cell_rejected():
    cu = bulk("Cu", "fcc", a=3.615)
    cu.set_cell(-np.array(cu.cell), scale_atoms=True)
    with pytest.raises(InputError):
        atoms_to_kkr_lattice(cu)


def test_nonperiodic_rejected():
    cu = bulk("Cu", "fcc", a=3.615)
    cu.pbc = [True, True, False]
    with pytest.raises(InputError):
        atoms_to_kkr_lattice(cu)


def test_atmicx_are_fractional():
    co = bulk("Co", "hcp", a=2.507, c=4.07)
    param, toa = atoms_to_kkr_param(co)
    assert param["natm"] == 2
    frac = co.get_scaled_positions(wrap=True)
    for row, f in zip(param["atmicx"], frac):
        assert [float(x[:-1]) for x in row[:3]] == pytest.approx(list(f))
        assert [x[-1] for x in row[:3]] == ["a", "b", "c"]
        assert row[3] == toa[0]


def test_inputcard_text(tmp_path):
    nife = set_occupancy(bulk("Ni", "fcc", a=3.571), {0: {"Fe": 0.1, "Ni": 0.9}})
    calc = AkaiKKR(directory=str(tmp_path), command="true", sdftyp="pbeasa", bzqlty=8)
    dic, toa = calc.make_input_dict(nife)
    dic["record"] = "init"
    f = io.StringIO()
    AkaikkrJob(str(tmp_path)).make_inputcard(dic, f)
    text = f.getvalue()
    assert "aux" in text
    assert "pbeasa" in text
    assert "26 10.0" in text and "28 90.0" in text
    assert "Fe0.1Ni0.9_0 2" in text


def test_spc_rejected(tmp_path):
    calc = AkaiKKR(directory=str(tmp_path), command="true", go="spc31")
    with pytest.raises(InputError):
        calc.make_input_dict(bulk("Cu", "fcc", a=3.615))


def test_check_state_sees_occupancy(tmp_path):
    a1 = set_occupancy(bulk("Ni", "fcc", a=3.571), {0: {"Fe": 0.1, "Ni": 0.9}})
    a2 = set_occupancy(bulk("Ni", "fcc", a=3.571), {0: {"Fe": 0.2, "Ni": 0.8}})
    calc = AkaiKKR(directory=str(tmp_path), command="true")
    calc.atoms = a1.copy()
    assert calc.check_state(a1.copy()) == []
    assert calc.check_state(a2) == all_changes


def test_record_reuse(tmp_path):
    calc = AkaiKKR(directory=str(tmp_path), command="true")
    cu = bulk("Cu", "fcc", a=3.615)
    dic, _ = calc.make_input_dict(cu)
    assert calc._decide_record(dic) == "init"
    calc._last_dic = dict(dic)
    assert calc._decide_record(dic) == "init"          # no potential file yet
    open(os.path.join(str(tmp_path), "pot.dat"), "w").close()
    dic2, _ = calc.make_input_dict(bulk("Cu", "fcc", a=3.70))
    assert calc._decide_record(dic2) == "2nd"          # same types, different volume
    dic3, _ = calc.make_input_dict(set_occupancy(cu, {0: {"Cu": 0.9, "Ni": 0.1}}))
    assert calc._decide_record(dic3) == "init"         # different composition
    calc.set(record="2nd")
    assert calc._decide_record(dic3) == "2nd"          # explicit


def test_atoms_grouped_by_type():
    # B, Fe in the Atoms; the inputcard lists Fe (lower electronegativity) first
    atoms = bulk("Fe", "bcc", a=2.87, cubic=True)
    atoms.set_chemical_symbols(["B", "Fe"])
    param, toa = atoms_to_kkr_param(atoms)
    assert toa == ["B_1", "Fe_0"]
    assert [row[3] for row in param["atmicx"]] == ["Fe_0", "B_1"]
    assert param["atom_order"] == [1, 0]
