"""begin_option round trip through specx (needs AKAIKKR_PROGRAM_PATH)."""
import os
import shutil
from copy import deepcopy

import pytest

from pyakaikkr import AkaikkrJob, KKRFailedExecutionError, KKRUnknownOptionError
from option_env import TESTS_DIR, needs_specx, needs_specx_cnd, specx_path

CU_PARAM = dict(
    go="go", potentialfile="pot.dat", brvtyp="fcc", a=6.82, **{"c/a": 1.0, "b/a": 1.0},
    alpha=90, beta=90, gamma=90, edelt=1e-3, ewidth=1.0, reltyp="sra", sdftyp="mjw",
    magtyp="nmag", record="init", outtyp="update", bzqlty=4, maxitr=3, pmix=0.02,
    ntyp=1, rmt=[1.0], field=[0.0], mxl=[2], type=["Cu"], ncmp=[1], anclr=[[29]],
    conc=[[100]], natm=1, atmicx=[["0.0a", "0.0b", "0.0c", "Cu"]])


@needs_specx
def test_go_with_option(tmp_path):
    job = AkaikkrJob(str(tmp_path))
    dic = deepcopy(CU_PARAM)
    dic["option"] = {"number_emesh": 5, "tol": 1e-5}
    job.make_inputcard(dic, "inputcard_go")
    job.run(specx_path("akaikkr"), "inputcard_go", "out_go.log")
    assert job.get_option("out_go.log") == {"mse": 5, "tol": 1e-5}
    assert job.get_emesh_param("out_go.log")["mse"] == 5
    assert job.unused_option("inputcard_go", "out_go.log") == set()


@needs_specx
def test_unknown_key_is_reported(tmp_path):
    job = AkaikkrJob(str(tmp_path))
    dic = deepcopy(CU_PARAM)
    dic["option"] = {"nosuchkey": 1}
    with pytest.raises(KKRUnknownOptionError):
        job.make_inputcard(dic, "inputcard_go")
    # strict=False path: specx itself stops and run() reports the token
    from pyakaikkr.option import normalize_option, format_option_card
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        card = format_option_card(normalize_option({"nosuchkey": 1}, strict=False))
    del dic["option"]
    job.make_inputcard(dic, "inputcard_go")
    with open(os.path.join(str(tmp_path), "inputcard_go"), "a") as f:
        f.write("\n".join(card))
    with pytest.raises(KKRFailedExecutionError, match="unknown token: nosuchkey="):
        job.run(specx_path("akaikkr"), "inputcard_go", "out_go.log")


CND_CU = os.path.join(TESTS_DIR, "akaikkr_cnd", "Cu")


@needs_specx_cnd
@pytest.mark.skipif(not os.path.isfile(os.path.join(CND_CU, "pot.dat")),
                    reason="tests/akaikkr_cnd/Cu/pot.dat (converged by testrun.py) is needed")
def test_cnd_dos_with_cemesr_ref(tmp_path):
    """akaikkr_cnd: dos with cemesr_ref= 0.75 instead of the build default 0.5.

    The DOS window is [EF - ref*ewidth, EF + (1-ref)*ewidth]; with ewidth=2 the
    mesh moves from [-1, 1] to [-1.5, 0.5], and cemesr_ref is echoed by the dos run.
    """
    src = AkaikkrJob(CND_CU)
    shutil.copyfile(os.path.join(CND_CU, "pot.dat"), os.path.join(str(tmp_path), "pot.dat"))
    job = AkaikkrJob(str(tmp_path))
    base = src.read_inputcard_option("inputcard_dos")
    assert base == {}  # the reference run used no option

    # 1. default ref of the cnd build (0.5): window [-1, 1]
    text = open(os.path.join(CND_CU, "inputcard_dos")).read()
    with open(os.path.join(str(tmp_path), "inputcard_dos0"), "w") as f:
        f.write(text)
    job.run(specx_path("akaikkr_cnd"), "inputcard_dos0", "out_dos0.log")
    e0, _ = job.get_dos_as_list("out_dos0.log")
    assert job.get_option("out_dos0.log") == {}
    assert min(e0) == pytest.approx(-1.0, abs=0.02)
    assert max(e0) == pytest.approx(1.0, abs=0.02)

    # 2. cemesr_ref= 0.75: window [-1.5, 0.5]
    with open(os.path.join(str(tmp_path), "inputcard_dos"), "w") as f:
        f.write(text + "\nbegin_option\n cemesr_ref= 0.75\nend_option\n")
    assert job.read_inputcard_option("inputcard_dos") == {"cemesr_ref": 0.75}
    job.run(specx_path("akaikkr_cnd"), "inputcard_dos", "out_dos.log")
    assert job.get_option("out_dos.log") == {"cemesr_ref": 0.75}
    assert job.unused_option("inputcard_dos", "out_dos.log") == set()
    e1, _ = job.get_dos_as_list("out_dos.log")
    assert len(e1) == len(e0)
    assert min(e1) == pytest.approx(-1.5, abs=0.02)
    assert max(e1) == pytest.approx(0.5, abs=0.02)
    assert job.get_emesh_param("out_dos.log")["mse"] == 201
