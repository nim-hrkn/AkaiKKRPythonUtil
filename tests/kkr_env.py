# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""shared paths and markers of the pytest suites (tests/ is on sys.path through pyproject.toml pythonpath).

- AKAIKKR_PROGRAM_PATH: directory that contains akaikkr/specx, akaikkr_cnd/specx, akaikkr_cpa2021v01/specx.
  Tests marked ``specx`` (``needs_specx``, ``needs_specx_cnd``) are skipped without it (see conftest.py).
- GAES_RUN_DIR: where the GAES specx tests keep their RUN directories (default tests/gaes/RUN, kept for inspection).
"""
import os
import pytest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(TESTS_DIR)
DATA_DIR = os.path.join(TESTS_DIR, "data")                # small fixtures (Cu outputs, an out_go.log)
TESTRUN_DIR = os.path.join(TESTS_DIR, "testrun")          # regression test sets (outputs made by testrun.py)
STRUCTURE_DIR = os.path.join(TESTRUN_DIR, "structure")    # CIF files of the test sets
DOCS_DIR = os.path.join(ROOT_DIR, "docs")

PROGRAM_PATH = os.environ.get("AKAIKKR_PROGRAM_PATH")
RUN_DIR = os.environ.get("GAES_RUN_DIR", os.path.join(TESTS_DIR, "gaes", "RUN"))

# the 2019 HEA production run (legacy layout tests)
RUN0 = os.path.join(ROOT_DIR, "..", "fukushima_HEA_run_exprlattice", "production_run", "run0")
RUN0_RUN = os.path.join(RUN0, "RUN")


def specx_path(code="akaikkr"):
    """path of <AKAIKKR_PROGRAM_PATH>/<code>/specx, or None when it is not available."""
    if PROGRAM_PATH is None:
        return None
    path = os.path.join(PROGRAM_PATH, code, "specx")
    return path if os.path.isfile(path) else None


def legacy_dir(key, polytyp, iew=0, ied=0, ipm=0):
    return os.path.join(RUN0_RUN, "key_{},ew_{:03d},ed_{:03d},polytyp_{},pm_{:03d}".format(key, iew, ied, polytyp, ipm))


# markers: the skip is applied in conftest.py (pytest_collection_modifyitems)
needs_specx = pytest.mark.specx("akaikkr")
needs_specx_cnd = pytest.mark.specx("akaikkr_cnd")
needs_run0 = pytest.mark.skipif(not os.path.isdir(RUN0_RUN), reason="2019 RUN directory not available")
