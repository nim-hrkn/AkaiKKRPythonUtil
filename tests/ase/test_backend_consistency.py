"""the cif (pymatgen) backend and the ase backend must give the same structure.

needs specx: set AKAIKKR_PROGRAM_PATH.
"""
import os

import numpy as np
import pytest

from conftest import needs_specx, specx_path, STRUCTURE_DIR

CIFS = ["Cu-Fm3m", "Fe-Im3m", "Co_P63mmc", "Ni-Fm3m", "NiFe-Fm3m", "AlMnFeCo-Im3m",
        "FeRh0.5Pt0.5", "FeB1.95-P6mmm", "GaAsVc-F43m", "Co2MnSi-Fm3m", "SmCo5_P6mmm",
        "LiH-Fm-3m"]

BOHR = 0.529177


def _lattice_bohr(param):
    """primitive lattice vectors in bohr from a kkr param dict."""
    from pyakaikkr.Cif2Kkr import _TranslationKKR
    a = param["a"]
    if param["brvtyp"] == "aux":
        return np.array([param["r1"], param["r2"], param["r3"]]) * a
    conv = np.diag([a, a * param["b/a"], a * param["c/a"]])
    if param["brvtyp"] in ("hcp",):
        gamma = np.deg2rad(param["gamma"])
        conv = np.array([[a, 0, 0], [a * np.cos(gamma), a * np.sin(gamma), 0], [0, 0, a * param["c/a"]]])
    return np.dot(_TranslationKKR.getMatrix(param["brvtyp"]), conv)


def _atom_types_by_position(param, lattice):
    """{type name: sorted cartesian positions} built from atmicx."""
    out = {}
    for row in param["atmicx"]:
        frac = np.array([float(x[:-1]) for x in row[:3]])
        out.setdefault(row[3], []).append(np.dot(frac, lattice))
    return out


@needs_specx
@pytest.mark.parametrize("cif", CIFS)
def test_same_structure(cif, tmp_path):
    from akaikkr_testscript.testrun_class import get_kkr_struc_from_cif
    from akaikkr_testscript.asestruc import get_kkr_struc_from_ase

    path = os.path.join(STRUCTURE_DIR, cif + ".cif")
    p_cif = get_kkr_struc_from_cif(path, specx_path(), False, use_bravais=True,
                                   remove_temperaryfiles=True, Vc="Og",
                                   directory=str(tmp_path / "cif"))
    p_ase = get_kkr_struc_from_ase(path, specx_path(), False, remove_temperaryfiles=True,
                                   Vc="Og", directory=str(tmp_path / "ase"))
    assert p_cif is not None

    assert p_ase["ntyp"] == p_cif["ntyp"]
    assert p_ase["natm"] == p_cif["natm"]

    # composition of each type: compare as multisets of (anclr, conc)
    def comps(p):
        return sorted(sorted((int(zz), round(float(cc), 6)) for zz, cc in zip(z, c))
                      for z, c in zip(p["anclr"], p["conc"]))
    assert comps(p_ase) == comps(p_cif)

    # the two primitive lattices span the same lattice: integer unimodular transform
    L_cif = _lattice_bohr(p_cif)
    L_ase = _lattice_bohr(p_ase)
    assert abs(abs(np.linalg.det(L_ase)) - abs(np.linalg.det(L_cif))) < 1e-6 * abs(np.linalg.det(L_cif))
    M = np.dot(L_ase, np.linalg.inv(L_cif))
    assert np.allclose(M, np.round(M), atol=1e-5)

    # each type of the ase backend corresponds to one type of the cif backend
    # (same composition and same set of positions modulo lattice translations)
    T_cif = _atom_types_by_position(p_cif, L_cif)
    T_ase = _atom_types_by_position(p_ase, L_ase)
    inv = np.linalg.inv(L_cif)
    matched = set()
    for ta, pos_a in T_ase.items():
        ia = p_ase["type"].index(ta)
        found = None
        for tc, pos_c in T_cif.items():
            ic = p_cif["type"].index(tc)
            ca = sorted((int(z), round(float(c), 6)) for z, c in zip(p_ase["anclr"][ia], p_ase["conc"][ia]))
            cc = sorted((int(z), round(float(c), 6)) for z, c in zip(p_cif["anclr"][ic], p_cif["conc"][ic]))
            if ca != cc:
                continue
            if len(pos_a) != len(pos_c):
                continue
            ok = True
            for pa in pos_a:
                hit = False
                for pc in pos_c:
                    d = np.dot(pa - pc, inv)
                    if np.allclose(d, np.round(d), atol=1e-4):
                        hit = True
                        break
                ok = ok and hit
            if ok:
                found = tc
                break
        assert found is not None, f"type {ta} of the ase backend has no counterpart"
        assert found not in matched
        matched.add(found)
        # the type order must agree too (automatic rmt depends on the order)
        assert p_cif["type"].index(found) == ia, f"type order differs: {p_ase['type']} vs {p_cif['type']}"
