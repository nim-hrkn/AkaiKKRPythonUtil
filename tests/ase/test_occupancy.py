"""tests of the occupancy -> type/ncmp/anclr/conc mapping (no specx needed)."""
import numpy as np
import pytest
from ase.build import bulk
from ase.calculators.calculator import InputError

from pyakaikkr.ase import set_occupancy, atoms_to_kkr_types, get_kind_occupancies, to_primitive


def test_pure_element():
    types, toa = atoms_to_kkr_types(bulk("Cu", "fcc", a=3.615))
    assert types["ntyp"] == 1
    assert types["ncmp"] == [1]
    assert types["anclr"] == [[29]]
    assert types["conc"] == [[100.0]]
    assert types["type"] == ["Cu_0"]
    assert toa == ["Cu_0"]


def test_binary_alloy():
    atoms = set_occupancy(bulk("Ni", "fcc", a=3.571), {0: {"Fe": 0.1, "Ni": 0.9}})
    assert atoms.info["occupancy"] == {"0": {"Fe": 0.1, "Ni": 0.9}}
    assert list(atoms.arrays["spacegroup_kinds"]) == [0]
    types, _ = atoms_to_kkr_types(atoms)
    assert types["ncmp"] == [2]
    assert types["anclr"] == [[26, 28]]        # Z ascending
    assert types["conc"] == [[10.0, 90.0]]
    assert types["type"] == ["Fe0.1Ni0.9_0"]


def test_quaternary_alloy_sum_is_100():
    atoms = set_occupancy(bulk("Fe", "bcc", a=2.87),
                          {0: {"Al": 0.25, "Mn": 0.25, "Fe": 0.25, "Co": 0.25}})
    types, _ = atoms_to_kkr_types(atoms)
    assert types["anclr"] == [[13, 25, 26, 27]]
    assert sum(types["conc"][0]) == 100.0


def test_vacancy():
    atoms = set_occupancy(bulk("Ni", "fcc", a=3.571), {0: {"Ni": 0.95}})
    types, _ = atoms_to_kkr_types(atoms, Vc="Og")
    assert types["ncmp"] == [2]
    assert types["anclr"] == [[28, 0]]
    assert types["conc"][0] == pytest.approx([95.0, 5.0])
    assert types["type"] == ["Ni0.95Vc0.05_0"]


@pytest.mark.parametrize("occ", [{"Ni": 1.2}, {"Fe": 0.6, "Ni": 0.6}, {"Ni": 0.0}])
def test_invalid_occupancy(occ):
    # set_occupancy validates, and so does atoms_to_kkr_types on hand-made info
    with pytest.raises(InputError):
        set_occupancy(bulk("Ni", "fcc", a=3.571), {0: occ})
    atoms = bulk("Ni", "fcc", a=3.571)
    atoms.info["occupancy"] = {"0": occ}
    with pytest.raises(InputError):
        atoms_to_kkr_types(atoms)


def test_symbol_not_in_occupancy():
    atoms = bulk("Ni", "fcc", a=3.571)
    atoms.info["occupancy"] = {"0": {"Fe": 0.5, "Co": 0.5}}
    with pytest.raises(InputError):
        get_kind_occupancies(atoms)


def test_missing_kind_key():
    atoms = bulk("Cu", "fcc", a=3.615, cubic=True)
    atoms.info["occupancy"] = {"0": {"Cu": 1.0}}
    atoms.new_array("spacegroup_kinds", np.array([0, 0, 0, 1]))
    with pytest.raises(InputError):
        get_kind_occupancies(atoms)


def test_type_splitting_by_symmetry():
    sc = bulk("Cu", "fcc", a=3.615, cubic=True) * (2, 1, 1)
    sc.positions[0, 2] += 0.05
    t_sym, _ = atoms_to_kkr_types(sc, type_mode="symmetry")
    t_kind, _ = atoms_to_kkr_types(sc, type_mode="kind")
    t_atom, _ = atoms_to_kkr_types(sc, type_mode="atom")
    assert t_sym["ntyp"] > 1
    assert t_kind["ntyp"] == 1
    assert t_atom["ntyp"] == len(sc)


def test_symmetry_does_not_merge_different_kinds():
    # B2-like CsCl cell: two bcc sites with different occupancies
    atoms = bulk("Fe", "bcc", a=2.87, cubic=True)
    atoms = set_occupancy(atoms, {1: {"Rh": 0.5, "Pt": 0.5}})
    types, toa = atoms_to_kkr_types(atoms)
    assert types["ntyp"] == 2
    assert types["anclr"] == [[26], [45, 78]]
    assert toa == ["Fe_0", "Rh0.5Pt0.5_1"]


def test_to_primitive_keeps_occupancy():
    conv = bulk("Ni", "fcc", a=3.571, cubic=True)
    conv = set_occupancy(conv, {i: {"Fe": 0.1, "Ni": 0.9} for i in range(4)})
    prim = to_primitive(conv)
    assert len(prim) == 1
    types, _ = atoms_to_kkr_types(prim)
    assert types["anclr"] == [[26, 28]]
    assert types["conc"] == [[10.0, 90.0]]
    assert abs(prim.get_volume() - conv.get_volume() / 4) < 1e-8


def test_type_params_override():
    atoms = set_occupancy(bulk("Ni", "fcc", a=3.571), {0: {"Fe": 0.1, "Ni": 0.9}})
    types, _ = atoms_to_kkr_types(atoms, type_params={"Fe": {"mxl": 3}})
    assert types["mxl"] == [3]
    types, _ = atoms_to_kkr_types(atoms, type_params={"Fe0.1Ni0.9_0": {"rmt": 0.4}})
    assert types["rmt"] == [0.4]


def test_type_order_by_electronegativity():
    # CsCl-type FeB-like cell: B (X=2.04) must come after Fe (X=1.83) whatever the atom order
    atoms = bulk("Fe", "bcc", a=2.87, cubic=True)
    atoms.set_chemical_symbols(["B", "Fe"])
    types, toa = atoms_to_kkr_types(atoms)
    assert types["type"] == ["Fe_0", "B_1"]
    assert toa == ["B_1", "Fe_0"]
    types, toa = atoms_to_kkr_types(atoms, type_order="appearance")
    assert types["type"] == ["B_0", "Fe_1"]
