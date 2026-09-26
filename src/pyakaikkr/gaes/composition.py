# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""single-site composition with arbitrary elements and ratios, and the single-site CPA input."""
import re
from copy import deepcopy
from dataclasses import dataclass
from typing import Tuple

from pymatgen.core.periodic_table import Element

from ..Error import GaesError
from ..option import normalize_option

TYPE_NAME_MAX_LEN = 40   # len_type of AkaiKKR 2022.0721 (m_param.f); longer names are silently truncated
VACANCY = "Vc"           # AkaiKKR vacancy, Z = 0

_TOKEN = re.compile(r"([A-Z][a-z]?)([0-9]*\.?[0-9]*)")


def _symbol_to_z(symbol):
    if symbol == VACANCY:
        return 0
    try:
        return Element(symbol).Z
    except ValueError as e:
        raise GaesError("unknown element symbol {!r}".format(symbol)) from e


def _z_to_symbol(z):
    return VACANCY if z == 0 else str(Element.from_Z(z))


def _fmt_fraction(x, ndigits):
    s = "{:.{}f}".format(x, ndigits).rstrip("0").rstrip(".")
    return s if s else "0"


@dataclass(frozen=True)
class SiteComposition:
    """elements and fractions (sum 1) of one CPA site, in the given order."""
    elements: Tuple[str, ...]
    fractions: Tuple[float, ...]

    def __post_init__(self):
        if len(self.elements) == 0 or len(self.elements) != len(self.fractions):
            raise GaesError("elements and fractions must be non-empty and of equal length")
        for s in self.elements:
            _symbol_to_z(s)
        if any(f <= 0 for f in self.fractions):
            raise GaesError("fractions must be positive: {}".format(self.fractions))
        if abs(sum(self.fractions) - 1.0) > 1e-6:
            raise GaesError("fractions must sum to 1: {} (sum {})".format(self.fractions, sum(self.fractions)))

    # ---- constructors ----
    @classmethod
    def from_elements(cls, elements):
        """equal fractions. elements: ["Al", "Si"] or "AlSi"."""
        if isinstance(elements, str):
            elements = [m.group(1) for m in _TOKEN.finditer(elements)]
        elements = tuple(elements)
        n = len(elements)
        if n == 0:
            raise GaesError("no element given")
        return cls(elements, tuple([1.0 / n] * n))

    @classmethod
    def from_dict(cls, comp):
        """{"Rh": 0.5, "Pt": 0.5} (fractions are normalized to sum 1)."""
        total = float(sum(comp.values()))
        if total <= 0:
            raise GaesError("fractions must be positive")
        return cls(tuple(comp.keys()), tuple(float(v) / total for v in comp.values()))

    @classmethod
    def from_type_name(cls, name):
        """"Rh0.5Pt0.5", "Rh0.5Pt0.5_1", "Fe", "B0.975Vc0.025". A trailing "_..." is dropped.
        Missing fractions mean 1 (single element) or equal shares when none is given."""
        body = name.split("_", 1)[0].strip()
        if not body:
            raise GaesError("empty type name")
        pos = 0
        elements, fractions = [], []
        for m in _TOKEN.finditer(body):
            if m.start() != pos:
                raise GaesError("cannot parse type name {!r} at {!r}".format(name, body[pos:]))
            pos = m.end()
            elements.append(m.group(1))
            fractions.append(float(m.group(2)) if m.group(2) not in ("", ".") else None)
        if pos != len(body):
            raise GaesError("cannot parse type name {!r} at {!r}".format(name, body[pos:]))
        if all(f is None for f in fractions):
            return cls.from_elements(elements)
        if any(f is None for f in fractions):
            raise GaesError("type name {!r}: every element needs a fraction, or none".format(name))
        return cls.from_dict(dict(zip(elements, fractions)))

    # ---- views ----
    @property
    def z(self):
        return tuple(_symbol_to_z(s) for s in self.elements)

    @property
    def conc(self):
        """concentrations in percent as AkaiKKR writes them; integers when exact."""
        out = []
        for f in self.fractions:
            c = 100.0 * f
            out.append(int(round(c)) if abs(c - round(c)) < 1e-9 else round(c, 6))
        return tuple(out)

    @property
    def is_equiatomic(self):
        return all(abs(f - self.fractions[0]) < 1e-9 for f in self.fractions)

    def type_name(self, suffix="", ndigits=3, max_len=TYPE_NAME_MAX_LEN):
        """AkaiKKR type name: element + fraction (omitted for a single element), plus suffix.

        Raises:
            GaesError: longer than max_len, or containing a space or comma.
        """
        if len(self.elements) == 1:
            name = self.elements[0]
        else:
            name = "".join(s + _fmt_fraction(f, ndigits) for s, f in zip(self.elements, self.fractions))
        name += suffix
        check_type_name(name, max_len=max_len)
        return name

    def key(self, ndigits=3):
        """file-system friendly label ("." -> "p")."""
        return self.type_name(ndigits=ndigits, max_len=10 ** 6).replace(".", "p")

    def __str__(self):
        return self.type_name(max_len=10 ** 6)


def check_type_name(name, max_len=TYPE_NAME_MAX_LEN):
    """AkaiKKR accepts at most 40 characters (truncated silently beyond) and no space or comma."""
    if len(name) > max_len:
        raise GaesError("type name {!r} has {} characters; AkaiKKR keeps only {}".format(
            name, len(name), max_len))
    if any(c.isspace() for c in name) or "," in name:
        raise GaesError("type name {!r} must not contain spaces or commas".format(name))
    return name


def make_single_site_param(comp, brvtyp, lattice="expr", go="go", ewidth=1.2, edelt=1e-4,
                           pmix=0.005, maxitr=500, bzqlty=10, sdftyp="pbe", reltyp="sra",
                           magtyp="mag", record="2nd", outtyp="update", rmt=1.0, mxl=2,
                           field=0.0, type_name=None, potentialfile="pot.dat", option=None):
    """AkaiKKR input dict of a single-site CPA (one type, one atom) for make_inputcard.

    Defaults are those of the 2019 HEA run (fukushima_HEA_run_exprlattice run0).

    Args:
        comp (SiteComposition): site composition.
        brvtyp (str): fcc, bcc, hcp, ...
        lattice (str or float, optional): "expr" writes a=1000000 (specx takes the
            concentration average of the experimental atomic volumes), "mjw" writes
            a=0 (MJW table), a float is a in bohr. Defaults to "expr".
        type_name (str, optional): AkaiKKR type name. Defaults to comp.type_name().
        option (dict, optional): begin_option keys (validated by normalize_option).

    Returns:
        dict: parameter dict.
    """
    from ..AkaiKkr import AkaikkrJob
    if type_name is None:
        type_name = comp.type_name()
    else:
        check_type_name(type_name)
    if lattice == "expr":
        a = 1000000
    elif lattice == "mjw":
        a = 0
    else:
        a = float(lattice)
        if a <= 0:
            raise GaesError("lattice must be 'expr', 'mjw' or a positive length in bohr")
    dic = deepcopy(AkaikkrJob("dummy").default)
    dic.update({
        "go": go, "potentialfile": potentialfile,
        "brvtyp": brvtyp, "a": a, "c/a": 1.0, "b/a": 1.0, "alpha": 90, "beta": 90, "gamma": 90,
        "edelt": edelt, "ewidth": ewidth, "reltyp": reltyp, "sdftyp": sdftyp, "magtyp": magtyp,
        "record": record, "outtyp": outtyp, "bzqlty": bzqlty, "maxitr": maxitr, "pmix": pmix,
        "ntyp": 1, "type": [type_name], "ncmp": [len(comp.elements)],
        "rmt": [rmt], "field": [field], "mxl": [mxl],
        "anclr": [list(comp.z)], "conc": [list(comp.conc)],
        "natm": 1, "atmicx": [["0.0a", "0.0b", "0.0c", type_name]],
    })
    if option:
        dic["option"] = normalize_option(option, strict=True)
    return dic
