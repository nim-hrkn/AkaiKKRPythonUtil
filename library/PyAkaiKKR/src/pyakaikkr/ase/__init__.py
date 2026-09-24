# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""ASE interface of AkaiKKR.

This subpackage needs the ase package. Install it with ``pip install ase``
(or ``pip install pyakaikkr[ase]``).
"""
try:
    import ase as _ase  # noqa: F401
except ImportError as e:  # pragma: no cover
    raise ImportError(
        "pyakaikkr.ase requires the ase package. Install it with 'pip install ase'.") from e

from .occupancy import set_occupancy, get_kind_occupancies, atoms_to_kkr_types, to_primitive
from .structure import (atoms_to_kkr_lattice, atoms_to_kkr_atmicx, atoms_to_kkr_param,
                        check_kkr_output_structure)
from .calculator import AkaiKKR

__all__ = ["AkaiKKR", "set_occupancy", "get_kind_occupancies", "atoms_to_kkr_types",
           "to_primitive", "atoms_to_kkr_lattice", "atoms_to_kkr_atmicx",
           "atoms_to_kkr_param", "check_kkr_output_structure"]
