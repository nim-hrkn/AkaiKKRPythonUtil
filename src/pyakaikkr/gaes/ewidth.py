# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""is E_F - ewidth_go anchored in a gap region, and if not, which ewidth to try next."""
from typing import Optional

import numpy as np

from .gap import GapRegion


ETH = 0.30     # Ry, minimum gap width accepted (2019 script)
EDIFF = 0.20   # Ry, required distance from the upper edge of the gap (2019 script)
MARGIN = 0.01  # Ry, extra margin of a new ewidth (2019 script)


def check_ewidth(regions, ewidth, eth=ETH, ediff=EDIFF):
    """the gap region (highest first) that anchors E_F - ewidth, or None.

    A region qualifies when its width exceeds eth, e1 < -ewidth < e2 and
    -ewidth < e2 - ediff.
    """
    for g in sorted(regions, key=lambda g: g.e2, reverse=True):
        if g.width > eth and g.e1 < -ewidth < g.e2 and -ewidth < g.e2 - ediff:
            return g
    return None


def within_bounds(ewidth, min_ewidth=None, max_ewidth=None, tol=1e-9):
    """True when min_ewidth <= ewidth <= max_ewidth (None = no bound)."""
    return ((min_ewidth is None or ewidth >= min_ewidth - tol)
            and (max_ewidth is None or ewidth <= max_ewidth + tol))


def _clip_candidate(ew, lo, hi, min_ewidth, max_ewidth):
    """candidate ew (Ry) with -ew allowed in the open interval (lo, hi): if ew is outside
    [min_ewidth, max_ewidth] it is moved to the nearest bound provided -bound still lies
    inside (lo, hi); otherwise None."""
    if min_ewidth is not None and ew < min_ewidth:
        ew = min_ewidth
    if max_ewidth is not None and ew > max_ewidth:
        ew = max_ewidth
    return ew if lo < -ew < hi else None


def ewidth_candidates(regions, eth=ETH, ediff=EDIFF, margin=MARGIN, min_ewidth=None, max_ewidth=None):
    """new ewidth candidates, one per region wider than eth, highest region first:
    -(e2 - ediff - margin), moved to min_ewidth / max_ewidth when outside the limits as long as
    E_F - ewidth stays inside the region (below e2 - ediff)."""
    cands = []
    for g in sorted(regions, key=lambda g: g.e2, reverse=True):
        if g.width > eth:
            ew = _clip_candidate(-(g.e2 - ediff - margin), g.e1, g.e2 - ediff, min_ewidth, max_ewidth)
            if ew is not None:
                cands.append(ew)
    return cands


def choose_ewidth(regions, ewidth, eth=ETH, ediff=EDIFF, margin=MARGIN, min_ewidth=None, max_ewidth=None):
    """decision of the 2019 scheme (calc_new_ewidth).

    Args:
        regions (Sequence[GapRegion]): gap regions.
        ewidth (float or None): current ewidth_go. None means "no current value":
            the decision is never "old".

    Returns:
        (str, float or None, list[float]): ("old", ewidth, candidates) when the current
        ewidth is anchored, ("new", candidates[0], candidates) when another ewidth is
        proposed, ("fail", None, []) when no region is wide enough.
    """
    cands = ewidth_candidates(regions, eth=eth, ediff=ediff, margin=margin, min_ewidth=min_ewidth, max_ewidth=max_ewidth)
    # "old" only inside [min_ewidth, max_ewidth]: an ewidth that is anchored in a gap but violates
    # a bound (e.g. one derived from an orbital rule, section 15) must move
    if (ewidth is not None and within_bounds(ewidth, min_ewidth, max_ewidth)
            and check_ewidth(regions, ewidth, eth=eth, ediff=ediff) is not None):
        return "old", ewidth, cands
    if not cands:
        return "fail", None, []
    return "new", cands[0], cands


# ---------------------------------------------------------------- Method 2
from dataclasses import dataclass, field
from typing import List

from .gap import gap_regions

DOSTH2 = 1e-3          # fine threshold of Method 2
DOSTH2_RELAX = 2.0     # factor by which dosth2 is relaxed when no fine sub-region exists


@dataclass
class Decision:
    """result of choose_ewidth2 (also usable to report Method 1)."""
    flag: str                                  # old | new | fail
    ewidth: Optional[float]                    # current (old) or proposed (new) ewidth
    candidates: List[float] = field(default_factory=list)
    coarse: List[GapRegion] = field(default_factory=list)   # regions below dosth, wider than eth
    fine: List[GapRegion] = field(default_factory=list)     # sub-regions below dosth2 inside the coarse ones
    relaxed: bool = False                      # dosth2 was relaxed to find a fine sub-region
    gap_used: Optional[GapRegion] = None       # fine sub-region anchoring -ewidth when flag == old
    dosth2_used: Optional[float] = None
    window_limited: bool = False               # a coarse region touches the bottom of the dos window (ignored)
    reasons: List[str] = field(default_factory=list)   # why a region gave no candidate (fail diagnosis)


def _sorted_coarse(regions):
    """highest first; regions touching the bottom of the mesh (i1 == 0) last."""
    return sorted(regions, key=lambda g: (g.i1 == 0, -g.e2))


def choose_ewidth2(energy, curves, ewidth, dosth=2e-2, dosth2=DOSTH2, eth=ETH, ediff=EDIFF,
                   margin=MARGIN, dosth2_relax=DOSTH2_RELAX, min_ewidth=None, max_ewidth=None):
    """Method 2: two-threshold decision (docs/ewidth_tuning_scheme.md section 4.3).

    1. coarse regions: DOS < dosth, wider than eth (regions touching the mesh bottom are ranked last);
    2. fine sub-regions: DOS < dosth2 inside each coarse region (dosth2 relaxed by dosth2_relax
       when there is none); in a coarse region touching the bottom of the mesh only positions
       with at least eth of verified low DOS below them (-ewidth >= e_min + eth) are used, and
       Decision.window_limited is set when that leaves nothing (the caller may widen the window);
    3. the upper margin is measured from the coarse region's upper edge: top = e2 - ediff;
    4. "old" if -ewidth lies in a fine sub-region and below top, and ewidth is inside
       [min_ewidth, max_ewidth] (a bound violated by the current ewidth forces a move);
    5. candidates: -(min(f2, top) - margin) for each fine sub-region [f1, f2] when it lies above f1;
       a candidate outside [min_ewidth, max_ewidth] is moved to the bound if E_F - bound is still
       inside the sub-region (below top), otherwise dropped.

    Args:
        energy, curves: as for gap_regions.
        ewidth (float or None): current ewidth_go (None: never "old").

    Returns:
        Decision
    """
    coarse_all = gap_regions(energy, curves, dosth=dosth)
    coarse = _sorted_coarse([g for g in coarse_all if g.width > eth])
    dec = Decision(flag="fail", ewidth=None, coarse=coarse, dosth2_used=dosth2)
    if not coarse_all:
        dec.reasons.append("no region with DOS < {:g} in the window [{:.3f}, {:.3f}]".format(dosth, float(np.min(energy)), float(np.max(energy))))
    for g in coarse_all:
        if g.width <= eth:
            dec.reasons.append("region [{:.3f}, {:.3f}] (DOS < {:g}) is {:.3f} Ry wide < eth {:.2f}".format(g.e1, g.e2, dosth, g.width, eth))
    cands = []
    e_min = float(np.min(energy))
    may_keep = ewidth is not None and within_bounds(ewidth, min_ewidth, max_ewidth)
    for g in coarse:
        # a region touching the bottom of the dos window has a real upper edge but an unknown
        # lower edge: it may anchor or propose an ewidth only if at least eth of low DOS is
        # verified below E_F - ewidth (-ewidth >= e_min + eth); otherwise the window must be widened
        floor = e_min + eth if g.i1 == 0 else None
        top = g.e2 - ediff
        th2 = dosth2
        fine = [f for f in gap_regions(energy, curves, dosth=th2) if f.e1 >= g.e1 - 1e-12 and f.e2 <= g.e2 + 1e-12]
        if not fine and dosth2_relax and dosth2_relax > 1.0:
            th2 = dosth2 * dosth2_relax
            fine = [f for f in gap_regions(energy, curves, dosth=th2) if f.e1 >= g.e1 - 1e-12 and f.e2 <= g.e2 + 1e-12]
            if fine:
                dec.relaxed = True
                dec.dosth2_used = th2
        if not fine:
            sel = (energy >= g.e1) & (energy <= g.e2)
            dmin = float(min(np.min(c[sel]) for c in curves)) if np.any(sel) else float("nan")
            dec.reasons.append("coarse region [{:.3f}, {:.3f}]: no sub-region with DOS < {:g} (minimum DOS {:.2e})".format(g.e1, g.e2, th2, dmin))
        fine = sorted(fine, key=lambda f: -f.e2)
        dec.fine.extend(fine)
        for f in fine:
            lo = f.e1 if floor is None else max(f.e1, floor)
            if may_keep and dec.gap_used is None and lo < -ewidth < f.e2 and -ewidth < top:
                dec.gap_used = f
            hi = min(f.e2, top)
            if hi - margin > lo:
                # shallowest position inside the fine sub-region below top; moved to the
                # min/max bound when outside the limits as long as it stays inside (lo, hi)
                raw = -(hi - margin)
                c = _clip_candidate(raw, lo, hi + 1e-12, min_ewidth, max_ewidth)
                if c is not None:
                    cands.append(c)
                else:
                    dec.reasons.append("sub-region [{:.3f}, {:.3f}]: candidate {:.4f} is outside [min_ewidth, max_ewidth] = [{}, {}] "
                                       "and the bound is not inside the sub-region".format(f.e1, f.e2, raw, min_ewidth, max_ewidth))
            elif floor is not None:
                dec.window_limited = True   # usable part of a window-bottom region is too small
                dec.reasons.append("sub-region [{:.3f}, {:.3f}] touches the window bottom: less than eth {:.2f} of low DOS verified "
                                   "below E_F - ewidth (usable part above {:.3f} is empty)".format(f.e1, f.e2, eth, floor))
            else:
                # the sub-region lies within ediff of the coarse upper edge: ediff excludes it.
                # ediff < e2 - f1 - margin would leave room (candidate about -f1)
                ediff_max = g.e2 - lo - margin
                dec.reasons.append("sub-region [{:.3f}, {:.3f}] excluded by ediff {:.2f} (E_F - ewidth must be below {:.3f} = coarse upper edge "
                                   "{:.3f} - ediff); ediff < {:.3f} would give a candidate near {:.4f}".format(
                                       f.e1, f.e2, ediff, top, g.e2, ediff_max, -lo - margin if ediff_max > 0 else float("nan")))
    dec.candidates = cands
    if dec.gap_used is not None:
        dec.flag, dec.ewidth = "old", ewidth
    elif cands:
        dec.flag, dec.ewidth = "new", cands[0]
    return dec


def decide(method, energy, curves, ewidth, dosth=1e-3, dosth2=DOSTH2, eth=ETH, ediff=EDIFF, margin=MARGIN,
           dosth2_relax=DOSTH2_RELAX, min_ewidth=None, max_ewidth=None):
    """Method 1 or 2 with a common Decision result."""
    if method == 2:
        return choose_ewidth2(energy, curves, ewidth, dosth=dosth, dosth2=dosth2, eth=eth, ediff=ediff,
                              margin=margin, dosth2_relax=dosth2_relax, min_ewidth=min_ewidth, max_ewidth=max_ewidth)
    regions = gap_regions(energy, curves, dosth=dosth)
    flag, ew, cands = choose_ewidth(regions, ewidth, eth=eth, ediff=ediff, margin=margin, min_ewidth=min_ewidth,
                                    max_ewidth=max_ewidth)
    used = check_ewidth(regions, ewidth, eth=eth, ediff=ediff) if flag == "old" else None
    return Decision(flag=flag, ewidth=ew, candidates=list(cands), coarse=list(regions), fine=[], relaxed=False,
                    gap_used=used, dosth2_used=None, window_limited=any(g.i1 == 0 and g.width > eth for g in regions))
