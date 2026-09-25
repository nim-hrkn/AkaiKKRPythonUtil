import os
import pytest

PROGRAM_PATH = os.environ.get("AKAIKKR_PROGRAM_PATH")
TESTS_DIR = os.path.join(os.path.dirname(__file__), "..")
RUN0 = os.path.join(TESTS_DIR, "..", "..", "fukushima_HEA_run_exprlattice", "production_run", "run0")
RUN0_RUN = os.path.join(RUN0, "RUN")
# where the specx tests put their RUN directories (kept for inspection)
RUN_DIR = os.environ.get("GAES_RUN_DIR", os.path.join(os.path.dirname(__file__), "RUN"))


def specx_path(code="akaikkr"):
    if PROGRAM_PATH is None:
        return None
    path = os.path.join(PROGRAM_PATH, code, "specx")
    return path if os.path.isfile(path) else None


def legacy_dir(key, polytyp, iew=0, ied=0, ipm=0):
    return os.path.join(RUN0_RUN, "key_{},ew_{:03d},ed_{:03d},polytyp_{},pm_{:03d}".format(key, iew, ied, polytyp, ipm))


needs_specx = pytest.mark.skipif(specx_path() is None,
                                 reason="set AKAIKKR_PROGRAM_PATH (directory containing akaikkr/specx)")
needs_run0 = pytest.mark.skipif(not os.path.isdir(RUN0_RUN), reason="2019 RUN directory not available")
