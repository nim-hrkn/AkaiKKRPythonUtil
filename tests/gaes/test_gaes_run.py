"""GAES on real systems with the 2019 parameters (needs AKAIKKR_PROGRAM_PATH).

Systems (2019 keys, Bi 5d semicore at E-EF ~ -1.75 Ry, gap [-1.43, -1.18] above it):
AlSiRhBi fcc (13144583), AlSiGeBi fcc (13143283, Ge 3d + Bi 5d), AlScNiBi fcc (13212883).
The initial ewidth 1.6 (GAES_EWIDTH_INIT) puts E_F - ewidth inside the Bi 5d peak, so the
scheme has to move it into the gap. The RUN directories are kept under GAES_RUN_DIR
(default tests/gaes/RUN) for inspection.

The specx build is chosen by GAES_SPECX_CODE (default akaikkr). Systems containing Hf
need akaikkr_cpa2021v01 (see docs/ewidth_tuning_scheme.md section 13.1).
"""
import json
import os

import pytest

from pyakaikkr.gaes import Gaes, Layout, heakey_to_composition, make_single_site_param
from gaes_env import specx_path, RUN_DIR

SPECX_CODE = os.environ.get("GAES_SPECX_CODE", "akaikkr")
EWIDTH_INIT = float(os.environ.get("GAES_EWIDTH_INIT", "1.6"))
DOSTH = float(os.environ.get("GAES_DOSTH", "2e-2"))   # Method 1 threshold / Method 2 coarse threshold
METHOD = int(os.environ.get("GAES_METHOD", "2"))
DOSTH2 = float(os.environ.get("GAES_DOSTH2", "1e-3"))
SYSTEMS = [("13144583", "fcc"),   # AlSiRhBi: Bi 5d only
           ("13143283", "fcc"),   # AlSiGeBi: Ge 3d + Bi 5d
           ("13212883", "fcc")]   # AlScNiBi
if os.environ.get("GAES_SYSTEMS"):   # e.g. "13142182:fcc,13315082:fcc"
    SYSTEMS = [tuple(x.split(":")) for x in os.environ["GAES_SYSTEMS"].split(",")]
# GAES_SYSTEMS="13142182:fcc,13315082:fcc" overrides the list (AlSiScPb, AlGaSnPb: Pb 5d at -1.27 Ry)
if os.environ.get("GAES_SYSTEMS"):
    SYSTEMS = [tuple(x.split(":")) for x in os.environ["GAES_SYSTEMS"].split(",")]


# dos window of the build: akaikkr ref=0.75 (2019: ewidth_dos 3.0 -> bottom -2.25 Ry);
# cpa2021v01 has ref=0.5 and no cemesr_ref option, so 4.5 gives the same bottom.
REF, EWIDTH_DOS = (0.75, 3.0) if SPECX_CODE == "akaikkr" else (0.5, 4.5)


def _gaes(prefix, **kw):
    # initial values of fukushima_HEA_run_exprlattice/production_run/run0 (1306.run_scheme2.py):
    # dosth 2e-2 on the spin sum = the 2019 value 1e-2 on the spin average;
    # max_pm_iter is reduced from 20 to 2 to bound the test time (fresh retries and a bzqlty round follow).
    return Gaes(specx_path(SPECX_CODE), Layout(prefix, version=2), ewidth_init=EWIDTH_INIT, ewidth_dos=EWIDTH_DOS, ref=REF,
                method=METHOD, dosth=DOSTH, dosth2=DOSTH2, edelt_init=1e-4, edelt_steps=(1e-2, 1e-3, 1e-4), pmix_steps=(1e-2, 5e-3, 1e-3, 5e-4, 1e-4),
                pmix_init=0.005, maxitr_init=500, maxitr_2nd=200, maxitr_pm=300, max_pm_iter=2, with_j=False, **kw)


@pytest.mark.skipif(specx_path(os.environ.get("GAES_SPECX_CODE", "akaikkr")) is None,
                    reason="set AKAIKKR_PROGRAM_PATH (directory containing <GAES_SPECX_CODE>/specx)")
@pytest.mark.parametrize("key,polytyp", SYSTEMS)
def test_gaes_2019_systems(key, polytyp):
    comp = heakey_to_composition(key)
    params = {polytyp: make_single_site_param(comp, polytyp, type_name="HEA")}   # 2019 inputcard
    prefix = os.path.join(RUN_DIR, "{}_{}".format(comp.type_name(max_len=10 ** 6).replace(".", "p"), polytyp))
    g = _gaes(prefix)
    res = g.run(key, params)
    assert res.status != "error", res.message
    assert res.ewidth_tried[0] == EWIDTH_INIT
    assert res.judgements and res.judgements[0].step == "step1"
    j0 = res.judgements[0]
    assert j0.regions, "no low-DOS region found"
    # the decision must be consistent with the regions
    if j0.flag == "old":
        assert j0.gap_used is not None and j0.gap_used[0] < -EWIDTH_INIT < j0.gap_used[1]
    else:
        assert j0.gap_used is None
    if res.status == "finished":
        assert res.converged[polytyp] and res.gap_used[0] < -res.ewidth_final < res.gap_used[1]
        r = res.results[polytyp]
        assert r["converged"] and r["ewidth"] == pytest.approx(res.ewidth_final, abs=1e-3)   # out_go.log prints 3 decimals
    # the JSON is written and readable
    with open(os.path.join(prefix, "key_{}.json".format(key))) as f:
        d = json.load(f)
    assert d["status"] == res.status
    print("\n", key, comp.type_name(max_len=100), polytyp, "->", res.status, "ewidth", res.ewidth_final,
          "tried", res.ewidth_tried, "gap", res.gap_used, "flags", [j.flag for j in res.judgements])
