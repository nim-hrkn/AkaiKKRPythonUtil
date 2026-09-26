# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""directory names of the scheme: key_<key>,ew_<iew>[-<ewidth>],ed_<ied>[-<edelt>],polytyp_<p>,pm_<ipm>[-<pmix>]"""
import os
from dataclasses import dataclass, replace
from typing import Optional

from ..Error import GaesError


@dataclass(frozen=True)
class RunPoint:
    key: str
    iew: int
    ewidth: Optional[float]
    ied: int
    edelt: Optional[float]
    polytyp: str
    ipm: int
    pmix: Optional[float]

    def replace(self, **kw):
        return replace(self, **kw)


class Layout:
    """version 1: the 2019 names (no parameter values); version 2: values included."""

    def __init__(self, prefix="RUN", version=2, sep=","):
        if version not in (1, 2):
            raise ValueError("version must be 1 or 2")
        if sep in "_-" or not sep:
            raise ValueError("sep must not be '_' or '-' or empty")
        self.prefix = prefix
        self.version = version
        self.sep = sep

    def dirname(self, p):
        if self.version == 1:
            parts = ["key_{}".format(p.key), "ew_{:03d}".format(p.iew), "ed_{:03d}".format(p.ied),
                     "polytyp_{}".format(p.polytyp), "pm_{:03d}".format(p.ipm)]
        else:
            parts = ["key_{}".format(p.key), "ew_{:03d}-{:.4f}".format(p.iew, p.ewidth),
                     "ed_{:03d}-{:.0e}".format(p.ied, p.edelt), "polytyp_{}".format(p.polytyp),
                     "pm_{:03d}-{:.0e}".format(p.ipm, p.pmix)]
        return self.sep.join(parts)

    def path(self, p):
        return os.path.join(self.prefix, self.dirname(p))

    def parse(self, dirname):
        name = os.path.basename(dirname.rstrip("/"))
        items = {}
        for part in name.split(self.sep):
            if "_" not in part:
                raise GaesError("cannot parse directory name {!r}".format(name))
            k, v = part.split("_", 1)
            items[k] = v
        try:
            def idx_val(s):
                if "-" in s:
                    i, v = s.split("-", 1)
                    return int(i), float(v)
                return int(s), None
            iew, ewidth = idx_val(items["ew"])
            ied, edelt = idx_val(items["ed"])
            ipm, pmix = idx_val(items["pm"])
            return RunPoint(items["key"], iew, ewidth, ied, edelt, items["polytyp"], ipm, pmix)
        except (KeyError, ValueError) as e:
            raise GaesError("cannot parse directory name {!r}: {}".format(name, e)) from e

    def find(self, **fixed):
        """RunPoints of the existing directories under prefix, filtered by field values."""
        if not os.path.isdir(self.prefix):
            return []
        out = []
        for entry in os.scandir(self.prefix):
            if not entry.is_dir() or not entry.name.startswith("key_"):
                continue
            try:
                p = self.parse(entry.name)
            except GaesError:
                continue
            if all(getattr(p, k) == v for k, v in fixed.items()):
                out.append(p)
        return sorted(out, key=lambda p: (p.key, p.polytyp, p.iew, p.ied, p.ipm))
