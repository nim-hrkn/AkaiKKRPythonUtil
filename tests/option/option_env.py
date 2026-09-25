import os
import pytest

PROGRAM_PATH = os.environ.get("AKAIKKR_PROGRAM_PATH")
DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
TESTS_DIR = os.path.join(os.path.dirname(__file__), "..")


def specx_path(code="akaikkr"):
    if PROGRAM_PATH is None:
        return None
    path = os.path.join(PROGRAM_PATH, code, "specx")
    return path if os.path.isfile(path) else None


needs_specx = pytest.mark.skipif(
    specx_path("akaikkr") is None,
    reason="set AKAIKKR_PROGRAM_PATH (directory containing akaikkr/specx)")
needs_specx_cnd = pytest.mark.skipif(
    specx_path("akaikkr_cnd") is None,
    reason="set AKAIKKR_PROGRAM_PATH (directory containing akaikkr_cnd/specx)")
