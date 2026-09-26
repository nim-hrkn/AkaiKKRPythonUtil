# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""keys of the 2019 HEA run: two-digit atomic numbers concatenated ("13142122" = AlSiScTi), equal fractions."""
from ..Error import GaesError
from .composition import SiteComposition, _z_to_symbol, _symbol_to_z


def heakey_to_composition(key):
    """"13142122" or 13142122 -> SiteComposition(Al, Si, Sc, Ti), equal fractions."""
    key = str(key).strip()
    if len(key) == 0 or len(key) % 2 != 0 or not key.isdigit():
        raise GaesError("heakey must have an even number of digits: {!r}".format(key))
    z = [int(key[i:i + 2]) for i in range(0, len(key), 2)]
    return SiteComposition.from_elements([_z_to_symbol(v) for v in z])


def composition_to_heakey(comp):
    """SiteComposition (equal fractions) -> "13142122"."""
    if not comp.is_equiatomic:
        raise GaesError("heakey is defined for equal fractions only: {}".format(comp))
    return "".join("{:02d}".format(_symbol_to_z(s)) for s in comp.elements)


def load_heakeylist(path):
    """heakeylist0.csv (columns heakey, elements) -> list of heakey strings (leading zeros kept)."""
    import pandas as pd
    df = pd.read_csv(path, dtype={"heakey": str})
    return [str(k).strip() for k in df["heakey"]]
