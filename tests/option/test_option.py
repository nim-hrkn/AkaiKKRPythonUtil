"""pyakaikkr.option: key registry, validation, writing and reading of begin_option."""
import io
import os
import re
import warnings

import pytest

from pyakaikkr import (AkaikkrJob, OPTION_KEYS, KKRUnknownOptionError, KKROptionValueError,
                       canonical_name, list_option_keys, normalize_option, format_option_card,
                       parse_option_block, read_inputcard_option, read_inputcard_option_blocks,
                       parse_option_echo, find_option_error)
from kkr_env import DATA_DIR as SHARED_DATA_DIR, DOCS_DIR

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
DOC = os.path.join(DOCS_DIR, "akaikkr_option_keys.md")
EXAMPLE_BLOCK = """begin_option
 mse= 3
 cpaitr_show= True
 begin_klabel
  G X W K G L
 end_klabel
end_option
"""


# ---------- A. registry ----------
def test_registry_size_and_aliases():
    assert len(OPTION_KEYS) == 16
    assert canonical_name("number_emesh=") == "mse"
    assert canonical_name("thresh_go") == "tol"
    assert canonical_name("thresh_scf") == "tol"
    assert canonical_name("ndegree_cheb") == "ng"
    assert canonical_name("ndirection_cnd") == "ie"
    assert canonical_name("mse=") == "mse"
    with pytest.raises(KKRUnknownOptionError):
        canonical_name("nosuchkey")
    assert OPTION_KEYS["cemesr_ref"].default_for("akaikkr") == 0.75
    assert OPTION_KEYS["cemesr_ref"].default_for("akaikkr_cnd") == 0.5
    cnd = {k.name for k in list_option_keys("akaikkr_cnd")}
    akaikkr = {k.name for k in list_option_keys("akaikkr")}
    assert {"ie", "cpaitr_show", "cpaitr_tol", "ddos", "tempmu"} <= cnd - akaikkr


@pytest.mark.skipif(not os.path.isfile(DOC), reason="docs/akaikkr_option_keys.md not found")
def test_registry_matches_doc_table():
    names = set()
    for line in open(DOC):
        m = re.match(r"\|\s*`([a-z_0-9]+)=`", line)
        if m:
            names.add(m.group(1))
        if line.startswith("| `begin_klabel"):
            names.add("klabel")
    assert names == set(OPTION_KEYS)


# ---------- D. normalize / format ----------
def test_normalize_values_and_aliases():
    out = normalize_option({"number_emesh": 3, "cpaitr_show": True, "ddos": "false",
                            "tol": 1e-5, "klabel": ["G", "X"], "spmain_bnd2": "b"})
    assert out == {"mse": "3", "cpaitr_show": "T", "ddos": "F", "tol": "1e-05",
                   "klabel": ["G", "X"], "spmain_bnd2": "b"}


def test_normalize_errors_and_warnings():
    with pytest.raises(KKRUnknownOptionError):
        normalize_option({"nosuchkey": 1})
    with pytest.raises(KKRUnknownOptionError):
        normalize_option({"mse=5": ""})
    with pytest.raises(KKROptionValueError):
        normalize_option({"mse": "abc"})
    with pytest.raises(KKROptionValueError):
        normalize_option({"mse": 2.5})
    with pytest.raises(KKROptionValueError):
        normalize_option({"cpaitr_show": "maybe"})
    with pytest.raises(KKROptionValueError):
        normalize_option({"spmain_bnd2": "x" * 81})
    with pytest.raises(KKROptionValueError):
        normalize_option({"mse": 3, "number_emesh": 4})
    with pytest.warns(UserWarning, match="no effect in akaikkr"):
        normalize_option({"cpaitr_show": True}, code="akaikkr")
    with pytest.warns(UserWarning, match="unknown option key"):
        out = normalize_option({"nosuchkey": 1}, strict=False)
    assert out == {"nosuchkey": "1"}
    assert normalize_option(None) == {}


def test_format_option_card_is_unchanged():
    lines = format_option_card(normalize_option({"mse": 3, "cpaitr_show": True,
                                                 "klabel": ["G", "X", "W", "K", "G", "L"]}))
    assert "\n".join(lines) == "\nbegin_option\n mse= 3\n cpaitr_show= T\n begin_klabel\n G X W K G L\n end_klabel\nend_option\n"
    assert format_option_card({}) == []


def test_make_inputcard_writes_canonical_block_or_nothing(tmp_path):
    job = AkaikkrJob(str(tmp_path))
    dic = dict(job.default)
    dic["option"] = {"number_emesh": 5, "thresh_go": 1e-5}
    f = io.StringIO()
    job.make_inputcard(dic, f)
    assert f.getvalue().endswith("\n\nbegin_option\n mse= 5\n tol= 1e-05\nend_option\n")
    for opt in ({}, None):
        dic2 = dict(job.default)
        if opt is not None:
            dic2["option"] = opt
        f = io.StringIO()
        job.make_inputcard(dic2, f)
        assert "begin_option" not in f.getvalue()
    dic["option"] = {"nosuchkey": 1}
    with pytest.raises(KKRUnknownOptionError):
        job.make_inputcard(dic, io.StringIO())


# ---------- B. parse ----------
def test_parse_example_block():
    assert parse_option_block(EXAMPLE_BLOCK) == {
        "mse": 3, "cpaitr_show": True, "klabel": ["G", "X", "W", "K", "G", "L"]}
    assert parse_option_block(EXAMPLE_BLOCK, typed=False)["mse"] == "3"


def test_parse_measured_inputcard():
    path = os.path.join(DATA_DIR, "inputcard_go_option")
    with pytest.warns(UserWarning, match="unknown array label foo"):
        opt = read_inputcard_option(path, strict=False)
    assert opt == {"mse": 5, "tol": 1e-5, "critic": -2.0, "cemesr_ref": 0.6,
                   "klabel": ["G", "X", "W"], "foo": ["1", "2"]}
    with pytest.raises(KKRUnknownOptionError, match="unknown array label foo"):
        read_inputcard_option(path, strict=True)


def test_parse_errors():
    with pytest.raises(KKRUnknownOptionError, match="unknown token: mse=5"):
        parse_option_block("begin_option\n mse=5\nend_option\n")
    with pytest.raises(KKRUnknownOptionError, match="end_option"):
        parse_option_block("begin_option\n mse= 5\n")
    with pytest.raises(KKRUnknownOptionError, match="unknown token: nosuch="):
        parse_option_block("begin_option\n nosuch= 5\nend_option\n")
    with pytest.warns(UserWarning):
        assert parse_option_block("begin_option\n nosuch= 5\nend_option\n", strict=False) == {"nosuch": "5"}
    with pytest.raises(KKRUnknownOptionError, match="begin_option not found"):
        parse_option_block("go pot.dat\n")
    assert read_inputcard_option("#--- go\ngo pot.dat\n") == {}


def test_round_trip_and_two_blocks(tmp_path):
    job = AkaikkrJob(str(tmp_path))
    dic = dict(job.default)
    dic["option"] = {"number_emesh": 5, "cpaitr_show": True, "klabel": ["G", "X"]}
    job.make_inputcard(dic, "inputcard_go")
    assert job.read_inputcard_option("inputcard_go") == {"mse": 5, "cpaitr_show": True,
                                                         "klabel": ["G", "X"]}
    text = open(os.path.join(str(tmp_path), "inputcard_go")).read()
    text += "\nbegin_option\n mse= 7\n tol= 1e-4\nend_option\n"
    blocks = read_inputcard_option_blocks(text)
    assert [b["mse"] for b in blocks] == [5, 7]
    assert read_inputcard_option(text) == {"mse": 7, "cpaitr_show": True, "klabel": ["G", "X"],
                                           "tol": 1e-4}


# ---------- C. echo ----------
def test_echo_of_measured_output():
    job = AkaikkrJob(DATA_DIR)
    assert job.get_option("out_go_option.log") == {"tol": 1e-5, "mse": 5, "critic": -2.0}
    assert job.get_emesh_param("out_go_option.log") == {"meshr": 400, "mse": 5, "ng": 21, "mxl": 3}
    assert job.check_option_error("out_go_option.log") is None
    # cemesr_ref was given but is not used by a go run; klabel is never echoed
    assert job.unused_option("inputcard_go_option", "out_go_option.log") == {"cemesr_ref"}
    assert parse_option_echo(" optnwrt:cpaitr_show  T\n optnwrt:tol  1.000000000000000E-005\n"
                             " optnwrt:ref  0.750000000000000\n optnwrt:spmain_bnd2 b\n") == {
        "cpaitr_show": True, "tol": 1e-5, "cemesr_ref": 0.75, "spmain_bnd2": "b"}


def test_echo_error_output():
    job = AkaikkrJob(DATA_DIR)
    assert job.check_option_error("out_go_badoption.log") == "unknown token: mse=5"
    assert job.get_option("out_go_badoption.log") == {}
    assert find_option_error(["x", ' failed to read "end_option", but found EOF']).startswith("failed")


CU_GO = os.path.join(SHARED_DATA_DIR, "Cu", "out_go.log")


@pytest.mark.skipif(not os.path.isfile(CU_GO), reason="tests/data/Cu/out_go.log not found")
def test_output_without_option():
    job = AkaikkrJob(os.path.dirname(CU_GO))
    assert job.get_option("out_go.log") == {}
    assert job.get_emesh_param("out_go.log")["ng"] == 21


def test_echo_of_cnd_dos_with_ref():
    # akaikkr_cnd dos run with cemesr_ref= 0.75 (echoed by specx as 'optnwrt:ref')
    job = AkaikkrJob(DATA_DIR)
    assert job.read_inputcard_option("inputcard_dos_cnd_ref") == {"cemesr_ref": 0.75}
    assert job.get_option("out_dos_cnd_ref.log") == {"cemesr_ref": 0.75}
    assert job.unused_option("inputcard_dos_cnd_ref", "out_dos_cnd_ref.log") == set()
    e, _ = job.get_dos_as_list("out_dos_cnd_ref.log")
    assert min(e) == pytest.approx(-1.5, abs=0.02) and max(e) == pytest.approx(0.5, abs=0.02)
