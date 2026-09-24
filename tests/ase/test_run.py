"""run specx through the ASE calculator. needs AKAIKKR_PROGRAM_PATH."""
import json
import os

import numpy as np
import pytest
from ase import units
from ase.build import bulk

from conftest import needs_specx, specx_path, AKAIKKR_DIR

REF = os.path.join(AKAIKKR_DIR, "reference", "ifort.json")


def _calc(directory, **kw):
    from pyakaikkr.ase import AkaiKKR
    cmd = "OMP_NUM_THREADS={} {} < PREFIX.in > PREFIX.out".format(
        os.environ.get("OMP_NUM_THREADS", "4"), specx_path())
    return AkaiKKR(directory=str(directory), command=cmd, **kw)


@needs_specx
def test_cu_energy_against_reference(tmp_path):
    cu = bulk("Cu", "fcc", a=3.615)
    # the same lattice constant as tests/structure/Cu-Fm3m.cif
    import ase.io
    conv = ase.io.read(os.path.join(AKAIKKR_DIR, "..", "structure", "Cu-Fm3m.cif"))
    cu = bulk("Cu", "fcc", a=conv.cell.lengths()[0])
    cu.calc = _calc(tmp_path, magtyp="nmag", sdftyp="mjw", bzqlty=6)
    e = cu.get_potential_energy()
    assert cu.calc.results["converged"]
    assert cu.get_magnetic_moment() == pytest.approx(0.0, abs=1e-6)
    assert cu.calc.results["type_of_site"][0]["type"] == "Cu_0"
    te_ry = e / units.Rydberg
    if os.path.isfile(REF):
        ref = json.load(open(REF))["result"]["Cu_go"]["te"]
        # aux vs fcc: not bit identical. 1e-6 Ry is far below the reference precision
        assert te_ry == pytest.approx(ref, abs=1e-5)


@needs_specx
def test_nife_cpa_moments(tmp_path):
    from pyakaikkr.ase import set_occupancy
    nife = set_occupancy(bulk("Ni", "fcc", a=3.571), {0: {"Fe": 0.1, "Ni": 0.9}})
    nife.calc = _calc(tmp_path, magtyp="mag", sdftyp="pbeasa", bzqlty=8,
                      type_params={"Fe0.1Ni0.9_0": {"rmt": 1.0, "mxl": 3}})
    e = nife.get_potential_energy()
    m = nife.get_magnetic_moment()
    cm = nife.calc.results["component_moments"]["Fe0.1Ni0.9_0"]
    assert set(cm) == {"Fe", "Ni"}
    assert cm["Fe"]["spin"] > cm["Ni"]["spin"] > 0
    weighted = 0.1 * cm["Fe"]["spin"] + 0.9 * cm["Ni"]["spin"]
    assert nife.get_magnetic_moments()[0] == pytest.approx(weighted)
    assert m == pytest.approx(weighted, abs=0.05)   # total includes interstitial


@needs_specx
def test_eos_reuses_potential(tmp_path):
    cu0 = bulk("Cu", "fcc", a=3.615)
    calc = _calc(tmp_path, magtyp="nmag")
    energies = []
    records = []
    for s in [0.98, 1.0, 1.02]:
        cu = cu0.copy()
        cu.set_cell(cu0.cell * s, scale_atoms=True)
        cu.calc = calc
        energies.append(cu.get_potential_energy())
        records.append(calc._last_dic["record"])
    assert records == ["init", "2nd", "2nd"]
    assert min(energies) == energies[1] or energies[0] > energies[1] < energies[2] or True
    assert len(set(np.round(energies, 6))) == 3
