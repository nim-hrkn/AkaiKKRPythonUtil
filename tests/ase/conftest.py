import os
import pytest

PROGRAM_PATH = os.environ.get("AKAIKKR_PROGRAM_PATH")
STRUCTURE_DIR = os.path.join(os.path.dirname(__file__), "..", "structure")
AKAIKKR_DIR = os.path.join(os.path.dirname(__file__), "..", "akaikkr")


def specx_path():
    if PROGRAM_PATH is None:
        return None
    path = os.path.join(PROGRAM_PATH, "akaikkr", "specx")
    return path if os.path.isfile(path) else None


needs_specx = pytest.mark.skipif(
    specx_path() is None,
    reason="set AKAIKKR_PROGRAM_PATH (directory containing akaikkr/specx)")
