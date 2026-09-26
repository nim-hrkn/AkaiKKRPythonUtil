# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""skip the tests marked ``specx`` (``@pytest.mark.specx("<build>")``) when that specx binary is not available.

Run only the offline tests with ``pytest -m "not specx"``.
"""
import pytest

from kkr_env import specx_path


def pytest_collection_modifyitems(config, items):
    for item in items:
        marker = item.get_closest_marker("specx")
        if marker is None:
            continue
        code = marker.args[0] if marker.args else "akaikkr"
        if specx_path(code) is None:
            item.add_marker(pytest.mark.skip(
                reason="set AKAIKKR_PROGRAM_PATH (directory containing {}/specx)".format(code)))
