# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""is E_F - ewidth_go anchored in a gap region, and if not, which ewidth to try next."""
from typing import Optional

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


def ewidth_candidates(regions, eth=ETH, ediff=EDIFF, margin=MARGIN, min_ewidth=None, max_ewidth=None):
    """new ewidth candidates, one per region wider than eth, highest region first:
    -(e2 - ediff - margin). Candidates below min_ewidth are dropped."""
    cands = []
    for g in sorted(regions, key=lambda g: g.e2, reverse=True):
        if g.width > eth:
            ew = -(g.e2 - ediff - margin)
            if (min_ewidth is None or ew >= min_ewidth) and (max_ewidth is None or ew <= max_ewidth):
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
    if ewidth is not None and check_ewidth(regions, ewidth, eth=eth, ediff=ediff) is not None:
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


def _sorted_coarse(regions):
    """highest first; regions touching the bottom of the mesh (i1 == 0) last."""
    return sorted(regions, key=lambda g: (g.i1 == 0, -g.e2))


def choose_ewidth2(energy, curves, ewidth, dosth=2e-2, dosth2=DOSTH2, eth=ETH, ediff=EDIFF,
                   margin=MARGIN, dosth2_relax=DOSTH2_RELAX, min_ewidth=None, max_ewidth=None):
    """Method 2: two-threshold decision (docs/ewidth_tuning_scheme.md section 4.3).

    1. coarse regions: DOS < dosth, wider than eth (regions touching the mesh bottom are ranked last);
    2. fine sub-regions: DOS < dosth2 inside each coarse region (dosth2 relaxed by dosth2_relax
       when there is none); a coarse region touching the bottom of the mesh is ignored
       (Decision.window_limited) because its lower edge is the window, not a band edge;
    3. the upper margin is measured from the coarse region's upper edge: top = e2 - ediff;
    4. "old" if -ewidth lies in a fine sub-region and below top;
    5. candidates: -(min(f2, top) - margin) for each fine sub-region [f1, f2] when it lies above f1.

    Args:
        energy, curves: as for gap_regions.
        ewidth (float or None): current ewidth_go (None: never "old").

    Returns:
        Decision
    """
    coarse_all = gap_regions(energy, curves, dosth=dosth)
    coarse = _sorted_coarse([g for g in coarse_all if g.width > eth])
    dec = Decision(flag="fail", ewidth=None, coarse=coarse, dosth2_used=dosth2)
    cands = []
    for g in coarse:
        if g.i1 == 0:
            # touches the bottom of the dos window: its lower edge is the window, not a band edge,
            # so it is neither an anchor nor a source of candidates (the caller may widen the window)
            dec.window_limited = True
            continue
        top = g.e2 - ediff
        th2 = dosth2
        fine = [f for f in gap_regions(energy, curves, dosth=th2) if f.e1 >= g.e1 - 1e-12 and f.e2 <= g.e2 + 1e-12]
        if not fine and dosth2_relax and dosth2_relax > 1.0:
            th2 = dosth2 * dosth2_relax
            fine = [f for f in gap_regions(energy, curves, dosth=th2) if f.e1 >= g.e1 - 1e-12 and f.e2 <= g.e2 + 1e-12]
            if fine:
                dec.relaxed = True
                dec.dosth2_used = th2
        fine = sorted(fine, key=lambda f: -f.e2)
        dec.fine.extend(fine)
        for f in fine:
            if ewidth is not None and dec.gap_used is None and f.e1 < -ewidth < f.e2 and -ewidth < top:
                dec.gap_used = f
            hi = min(f.e2, top)
            if hi - margin > f.e1:
                c = -(hi - margin)
                if (min_ewidth is None or c >= min_ewidth) and (max_ewidth is None or c <= max_ewidth):
                    cands.append(c)
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
