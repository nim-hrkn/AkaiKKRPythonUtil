# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""per-orbital valence / core rules and the ewidth range they imply (docs/ewidth_tuning_scheme.md section 15).

A rule "Rb4p=valence" (alias "occupied") asks for the Rb 4p state inside the energy contour:
ewidth >= |E_4p - E_F| + ediff.  "Rb4p=core" (alias "unoccupied") asks for it below the contour:
ewidth <= |E_4p - E_F| - ediff.  The levels come from the component blocks of a go output while
GAES runs, or from the tables in ``data/`` before the first go.
"""
import csv
import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

from ..AkaiKkr import AkaikkrJob
from ..Error import GaesError

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
TABLE_CONVERGED = os.path.join(DATA_DIR, "converged_core_levels_2019.csv")   # E - E_F, 2019 RUN medians
TABLE_ATOMIC = os.path.join(DATA_DIR, "atomic_core_levels_dsp.csv")         # absolute, go=dsp atomic start
EF_ASSUMED = 0.6    # Ry, E_F assumed when only the atomic table is available

ROLES = {"valence": "valence", "occupied": "valence", "core": "core", "unoccupied": "core"}
_RULE = re.compile(r"^\s*([A-Z][a-z]?)\s*([1-7][spdf])\s*(?:[=:]|\s)\s*([A-Za-z]+)\s*$")


@dataclass(frozen=True)
class OrbitalRule:
    element: str        # "Rb"
    orbital: str        # "4p"
    role: str           # "valence" | "core"

    @property
    def key(self):
        return self.element + self.orbital

    def __str__(self):
        return "{}{}={}".format(self.element, self.orbital, self.role)


@dataclass
class Level:
    """one core-configured state of one element, E - E_F in Ry (spin / component spread kept)."""
    element: str
    orbital: str
    e: float                    # representative value (mean)
    e_min: float                # deepest value seen (spin down/up, several components)
    e_max: float                # shallowest value seen
    star: Optional[bool]        # None when unknown (tables)
    source: str                 # "go" | "table_converged" | "table_atomic"

    @property
    def key(self):
        return self.element + self.orbital

    def as_list(self):
        return [self.e, self.star]


@dataclass
class Bounds:
    min_ewidth: Optional[float]
    max_ewidth: Optional[float]
    details: List[dict] = field(default_factory=list)   # {"rule", "level", "bound", "side"}

    def as_list(self):
        return [self.min_ewidth, self.max_ewidth]


def parse_orbital_rule(text):
    """"Rb4p=valence", "Se 4s:core", "Bi6s occupied" -> OrbitalRule (role names are case-insensitive)."""
    if isinstance(text, OrbitalRule):
        return text
    m = _RULE.match(str(text))
    if m is None:
        raise GaesError("cannot parse orbital rule {!r} (expected e.g. Rb4p=valence, Se4s=core, "
                        "Bi6s=occupied, In4d=unoccupied)".format(text))
    element, orbital, role = m.group(1), m.group(2), m.group(3).lower()
    if role not in ROLES:
        raise GaesError("unknown role {!r} in orbital rule {!r}: use valence|occupied|core|unoccupied".format(role, text))
    return OrbitalRule(element, orbital, ROLES[role])


def parse_orbital_rules(items):
    """list of rule strings (or OrbitalRule) -> list of OrbitalRule; a contradiction (the same orbital
    asked to be valence and core) is a GaesError."""
    rules = []
    for it in items or []:
        r = parse_orbital_rule(it)
        for other in rules:
            if other.key == r.key and other.role != r.role:
                raise GaesError("orbital {} is asked to be both {} and {}".format(r.key, other.role, r.role))
        if r not in rules:
            rules.append(r)
    return rules


def _merge(levels, element, orbital, values, star, source):
    """accumulate values (E - E_F) of one (element, orbital) into levels."""
    key = element + orbital
    if key in levels:
        lv = levels[key]
        lv.e_min = min(lv.e_min, min(values))
        lv.e_max = max(lv.e_max, max(values))
        lv.e = 0.5 * (lv.e_min + lv.e_max)
        if star is not None:
            lv.star = star if lv.star is None else (lv.star or star)
    else:
        levels[key] = Level(element, orbital, sum(values) / len(values), min(values), max(values), star, source)


def levels_from_records(records, levels: Optional[Dict[str, Level]] = None) -> Dict[str, Level]:
    """accumulate AkaikkrJob.get_core_levels_by_component records (dicts with element, orbital,
    e_minus_ef_Ry, star) into a {key: Level} dict."""
    levels = {} if levels is None else levels
    for rec in records:
        _merge(levels, rec["element"], rec["orbital"], [float(rec["e_minus_ef_Ry"])], rec.get("star"), "go")
    return levels


def levels_from_go(outfiles) -> Dict[str, Level]:
    """core levels (E - E_F, both spins, every component) of one or several go outputs.

    Args:
        outfiles: path of an out_go.log, or a sequence of paths / (AkaikkrJob, filename) pairs.
    """
    if isinstance(outfiles, (str, tuple)):
        outfiles = [outfiles]
    levels: Dict[str, Level] = {}
    for item in outfiles:
        if isinstance(item, tuple):
            job, fname = item
        else:
            d, fname = os.path.split(item)
            job = AkaikkrJob(d or ".")
        levels_from_records(job.get_core_levels_by_component(fname), levels)
    return levels


def levels_from_tables(elements: Sequence[str], ef_assumed=EF_ASSUMED) -> Dict[str, Level]:
    """E - E_F of the core-configured states of the given elements before any go has run:
    the converged 2019 table (median / min / max) where available, else the atomic go=dsp table
    minus ef_assumed (0.2-0.5 Ry too deep; use only to choose the first ewidth)."""
    wanted = set(elements)
    levels: Dict[str, Level] = {}
    with open(TABLE_CONVERGED) as f:
        for row in csv.DictReader(f):
            if row["element"] in wanted:
                levels[row["element"] + row["orbital"]] = Level(
                    row["element"], row["orbital"], float(row["E_minus_EF_median_Ry"]),
                    float(row["E_minus_EF_min_Ry"]), float(row["E_minus_EF_max_Ry"]),
                    float(row["valence_star_fraction"]) > 0.5, "table_converged")
    with open(TABLE_ATOMIC) as f:
        for row in csv.DictReader(f):
            key = row["element"] + row["orbital"]
            if row["element"] in wanted and row["orbital"] and row["energy_Ry"] and key not in levels:
                e = float(row["energy_Ry"]) - ef_assumed
                levels[key] = Level(row["element"], row["orbital"], e, e, e, None, "table_atomic")
    return levels


def levels_for_step0(elements: Sequence[str], ef_assumed=EF_ASSUMED) -> Dict[str, Level]:
    """levels used to choose the first ewidth: only the converged 2019 table.  The atomic table is
    0.2-0.5 Ry too deep and a wrong first range costs a go (a valence rule would start far too deep,
    a core rule too shallow), whereas the first judgement corrects any start from the real levels."""
    return {k: v for k, v in levels_from_tables(elements, ef_assumed).items() if v.source == "table_converged"}


def bounds_from_rules(rules: Sequence[OrbitalRule], levels: Dict[str, Level], ediff=0.2,
                      user_min=None, user_max=None, default_min=None, default_max=None, strict=True) -> Bounds:
    """ewidth range implied by the rules.

    valence: min_ewidth = |e_min| + ediff (deepest value seen); core: max_ewidth = |e_max| - ediff
    (shallowest).  Bounds from several rules are combined (max of mins, min of maxes).  Explicit
    user bounds are intersected; default bounds apply only when no rule contributes to that side
    and no user bound is given... except that as soon as any rule contributes, the defaults are
    dropped on both sides (the rule replaces the default range, section 15.2).
    A rule whose orbital is not in `levels` gives no bound for valence; for core it is a GaesError
    when strict (levels read from a go: the state is not core-configured, so ewidth cannot make it
    core) and no bound otherwise (step 0: the element is simply absent from the tables).
    """
    mins, maxs, details = [], [], []
    for r in rules:
        lv = levels.get(r.key)
        if lv is None:
            if r.role == "core" and strict:
                raise GaesError("{}: no core level for {} {} (not in the core configuration); "
                                "it cannot be made core by ewidth".format(r, r.element, r.orbital))
            details.append({"rule": str(r), "level": None, "bound": None, "side": "min" if r.role == "valence" else "max"})
            continue
        if r.role == "valence":
            b = abs(lv.e_min) + ediff
            mins.append(b)
            details.append({"rule": str(r), "level": lv.e_min, "bound": b, "side": "min"})
        else:
            b = abs(lv.e_max) - ediff
            maxs.append(b)
            details.append({"rule": str(r), "level": lv.e_max, "bound": b, "side": "max"})
    if mins or maxs:
        lo = max(mins + ([user_min] if user_min is not None else [])) if (mins or user_min is not None) else None
        hi = min(maxs + ([user_max] if user_max is not None else [])) if (maxs or user_max is not None) else None
    else:
        lo = user_min if user_min is not None else default_min
        hi = user_max if user_max is not None else default_max
    if lo is not None and hi is not None and lo > hi + 1e-9:
        raise GaesError("orbital rules contradict: min_ewidth {:.4f} > max_ewidth {:.4f} ({})".format(
            lo, hi, "; ".join("{}: {} {:.4f}".format(d["rule"], d["side"], d["bound"]) for d in details if d["bound"] is not None)))
    return Bounds(lo, hi, details)


def check_rules(rules: Sequence[OrbitalRule], levels: Dict[str, Level]) -> List[dict]:
    """rules that the go output does not satisfy: valence needs the '*' mark, core needs no mark.
    Returns [{"rule", "star", "level"}] (empty when everything matches)."""
    bad = []
    for r in rules:
        lv = levels.get(r.key)
        if lv is None or lv.star is None:
            if r.role == "core":
                bad.append({"rule": str(r), "star": None, "level": None})
            continue
        if (r.role == "valence") != bool(lv.star):
            bad.append({"rule": str(r), "star": lv.star, "level": lv.e})
    return bad


def initial_ewidth(ewidth_init, bounds: Bounds, ediff=0.2):
    """first ewidth: ewidth_init raised to min_ewidth + ediff (or the middle of a closed range) when
    a valence rule puts the lower bound above it.  It is never lowered for a core rule: a first go
    with a contour cutting into the valence band (Tl 5d core: 0.39 Ry) diverges or gives no gap,
    while a go at ewidth_init yields the DOS from which the next candidate inside the range is taken."""
    lo, hi = bounds.min_ewidth, bounds.max_ewidth
    if lo is None or ewidth_init >= lo - 1e-9:
        return ewidth_init
    if hi is not None and hi < lo + ediff:
        return round(0.5 * (lo + hi), 4)
    return round(lo + ediff, 4)


def levels_as_dict(levels: Dict[str, Level]) -> Dict[str, list]:
    return {k: [round(v.e, 4), v.star] for k, v in sorted(levels.items())}
