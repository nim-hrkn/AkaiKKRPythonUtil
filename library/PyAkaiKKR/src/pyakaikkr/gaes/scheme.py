# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""GAES (Gap-Anchored Ewidth Search): the sequential ewidth tuning scheme.

STEP1 runs go -> dos once per polytyp with the current ewidth_go and judges the DOS:
"old" (E_F - ewidth_go anchored in a gap) and all converged -> finished; "new" ->
retry STEP1 with the proposed ewidth (in a new directory); "fail" -> ewidth_fail.
STEP2 tightens the SCF (edelt from large to small, pmix steps under each edelt,
continuing from the previous potential while the run still converges) and judges
the DOS after each edelt.
"""
import json
import logging
import os
from copy import deepcopy
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional

from ..Error import GaesError, KKRFailedExecutionError
from .convergence import is_converging
from .ewidth import decide, ETH, EDIFF, MARGIN, DOSTH2, DOSTH2_RELAX
from .orbital import (parse_orbital_rules, levels_from_go, levels_for_step0, bounds_from_rules, check_rules,
                      initial_ewidth, levels_as_dict, EF_ASSUMED)
from .gap import dos_curves_from_outputs
from .layout import Layout, RunPoint
from .runner import KkrRunner

STATUS = ("finished", "ewidth_fail", "ewidth_exhausted", "not_converged", "error")


@dataclass
class Judgement:
    step: str                   # "step1" or "step2"
    iew: int
    ied: int
    ewidth: float
    converged: Dict[str, bool]
    regions: List[list]         # [[e1, e2], ...] of the AND of all polytyps
    flag: str                   # old | new | fail
    next_ewidth: Optional[float]
    candidates: List[float]
    gap_used: Optional[list]    # [e1, e2] anchoring -ewidth when flag == old
    directories: Dict[str, str]
    ewidth_dos: Optional[float] = None
    method: int = 1
    fine_regions: List[list] = field(default_factory=list)   # Method 2: sub-regions below dosth2
    relaxed: bool = False                                    # Method 2: dosth2 was relaxed
    dosth2_used: Optional[float] = None
    window_limited: bool = False
    orbital_levels: Dict[str, list] = field(default_factory=dict)   # 'Rb4p': [E - E_F, star] from the go outputs
    orbital_bounds: Optional[list] = None                            # [min_ewidth, max_ewidth] used for this judgement
    orbital_mismatch: List[dict] = field(default_factory=list)      # rules the go output does not satisfy


@dataclass
class KeyResult:
    key: str
    status: str = "error"
    polytyps: List[str] = field(default_factory=list)
    ewidth_final: Optional[float] = None
    ewidth_tried: List[float] = field(default_factory=list)
    gap_used: Optional[list] = None
    final: Dict[str, str] = field(default_factory=dict)      # polytyp -> directory
    converged: Dict[str, bool] = field(default_factory=dict)
    judgements: List[Judgement] = field(default_factory=list)
    results: Dict[str, dict] = field(default_factory=dict)   # polytyp -> KkrRunner.result()
    message: str = ""
    parameters: dict = field(default_factory=dict)

    def to_dict(self):
        return asdict(self)

    def save(self, path):
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=1, default=str)


DEFAULT_MIN_EWIDTH = 1.0   # Ry, user-set defaults (2026-09-25)
DEFAULT_MAX_EWIDTH = 2.0


class Gaes:
    """sequential Gap-Anchored Ewidth Search (port of run_scheme2 of the 2019 HEA run)."""

    def __init__(self, akaikkr_exe, layout=None, *, ewidth_init=1.2, ewidth_dos=3.0,
                 ewidth_dos_auto=True, ewidth_dos_max=4.5, ref=0.75, method=2, dosth=2e-2, dosth2=DOSTH2, dosth2_relax=DOSTH2_RELAX,
                 eth=ETH, ediff=EDIFF, margin=MARGIN,
                 min_ewidth=None, max_ewidth=None, orbitals=None, ef_assumed=EF_ASSUMED, edelt_init=1e-4, edelt_dos=1e-4, edelt_steps=(1e-2, 1e-3, 1e-4),
                 pmix_steps=(1e-2, 5e-3, 1e-3, 5e-4, 1e-4), pmix_init=0.005,
                 maxitr_init=500, maxitr_2nd=200, maxitr_pm=300, max_pm_iter=20, max_ew=10,
                 fresh_retry=True, bzqlty_steps=("+4",),
                 with_j=True, compat=False, tighten_before_fail=False, logger=None):
        self.akaikkr_exe = akaikkr_exe
        self.layout = layout if layout is not None else Layout("RUN", version=1 if compat else 2)
        self.ewidth_init = ewidth_init
        self.ewidth_dos = ewidth_dos
        self.ewidth_dos_auto = ewidth_dos_auto
        self.ewidth_dos_max = ewidth_dos_max   # widening limit when a gap touches the window bottom
        self.ref = ref
        self.method = 1 if compat else method
        self.dosth = 2e-2 if compat else dosth
        self.dosth2, self.dosth2_relax = dosth2, dosth2_relax
        self.eth, self.ediff, self.margin = eth, ediff, margin
        # ewidth is chosen inside [min_ewidth, max_ewidth] (candidates outside are moved to the bound,
        # "old" requires the current ewidth inside); gap regions are judged independently of the bounds.
        # None = the defaults 1.0 / 2.0 Ry, which per-orbital rules (section 15) replace.
        self._user_min, self._user_max = min_ewidth, max_ewidth
        self.min_ewidth = DEFAULT_MIN_EWIDTH if min_ewidth is None else min_ewidth
        self.max_ewidth = DEFAULT_MAX_EWIDTH if max_ewidth is None else max_ewidth
        self.orbitals = parse_orbital_rules(orbitals)   # e.g. ["Rb4p=valence", "Bi6s=core"]
        self.ef_assumed = ef_assumed                    # E_F assumed with the atomic table (step 0)
        self.edelt_init = edelt_init          # STEP1 (2019: 1e-4)
        # edelt of every dos used for the gap judgement, independent of the go edelt of STEP2:
        # a large edelt broadens the DOS and shrinks the gap regions, which would change the
        # criterion between stages. None = the go edelt (2019 behaviour, compat).
        self.edelt_dos = None if compat else edelt_dos
        self.edelt_steps = tuple(edelt_steps)  # STEP2, large to small
        self.pmix_steps = tuple(pmix_steps)
        self.pmix_init = pmix_init
        self.maxitr_init, self.maxitr_2nd, self.maxitr_pm = maxitr_init, maxitr_2nd, maxitr_pm
        self.max_pm_iter, self.max_ew = max_pm_iter, max_ew
        # potential handling (user rules, 2026-09-25): a new ewidth starts from a fresh potential;
        # changed convergence parameters continue from the previous potential; if that does not
        # converge, the same setting is retried from a fresh potential; then bzqlty is raised.
        self.fresh_retry = fresh_retry and not compat
        self.bzqlty_steps = tuple(bzqlty_steps) if not compat else ()
        self.with_j = with_j
        self.compat = compat
        # The DOS is judged whether or not the SCF converged (2019 behaviour, kept by design):
        # an ewidth outside the gap itself causes a density/potential mismatch that prevents
        # convergence, and the gap moves with ewidth. tighten_before_fail=True instead runs
        # STEP2 before accepting "fail" on an unconverged run.
        self.tighten_before_fail = tighten_before_fail and not compat
        self.logger = logger or logging.getLogger("pyakaikkr.gaes")

    # ---------------------------------------------------------------- helpers
    def parameters(self):
        return {k: getattr(self, k) for k in ("ewidth_init", "ewidth_dos", "ewidth_dos_auto", "ewidth_dos_max", "ref", "method", "dosth", "dosth2", "dosth2_relax",
                                               "eth", "ediff", "margin", "min_ewidth", "max_ewidth", "edelt_init", "edelt_dos", "edelt_steps", "pmix_steps",
                                               "pmix_init", "maxitr_init", "maxitr_2nd", "maxitr_pm",
                                               "max_pm_iter", "max_ew", "fresh_retry", "bzqlty_steps", "with_j", "compat",
                                               "tighten_before_fail")} | {"orbitals": [str(r) for r in self.orbitals],
                                                                          "ef_assumed": self.ef_assumed}

    def _bounds(self, levels, strict=True):
        """[min_ewidth, max_ewidth] for the decision: the orbital rules applied to `levels`, else the
        (user or default) bounds. strict=False at step 0 (a level missing from the tables is no error)."""
        return bounds_from_rules(self.orbitals, levels, self.ediff, self._user_min, self._user_max,
                                 self.min_ewidth, self.max_ewidth, strict=strict)

    def _ewidth_dos_for(self, ewidth_go):
        """dos window wide enough to see a gap of width eth below E_F - ewidth_go."""
        need = (ewidth_go + self.eth + self.ediff) / self.ref
        if self.ewidth_dos_auto and need > self.ewidth_dos:
            return round(need + 1e-9, 4)
        return self.ewidth_dos

    def _min_ewidth_dos(self, ewidth_go):
        """narrowest dos window still covering E_F - ewidth_go with margin ediff (used when
        specx cannot run a wide window)."""
        return (ewidth_go + self.ediff) / self.ref

    @staticmethod
    def _bzqlty(current, step):
        """bzqlty of a bzqlty_steps entry: "+4" adds to the current value, a number replaces it."""
        cur = int(current) if not isinstance(current, str) else int(str(current).rstrip("abcdefghijklmnopqrstuvwxyz") or 6)
        if isinstance(step, str) and step.startswith("+"):
            return cur + int(step[1:])
        return int(step)

    def _runner(self, point, param_go):
        dic = deepcopy(param_go)
        dic.update({"ewidth": point.ewidth, "edelt": point.edelt, "pmix": point.pmix})
        return KkrRunner(self.akaikkr_exe, self.layout.path(point), dic,
                         ewidth_dos=self._ewidth_dos_for(point.ewidth), edelt_dos=self.edelt_dos, compat=self.compat,
                         logger=self.logger)

    def _judge(self, step, iew, ied, ewidth, runners, converged):
        """gap regions of the AND of all polytyps' DOS, and the ewidth decision.

        If Method 2 finds no candidate only because the low-DOS region touches the bottom of the
        dos window, the dos is recomputed with a wider window (x1.5, up to ewidth_dos_max) and
        judged again."""
        widened = 0
        levels = {}
        mismatch = []
        if self.orbitals:
            # the levels of this go decide the mismatch; the bounds use every level seen so far in this
            # key (a state treated as core sits 0.7-0.9 Ry deeper than the same state treated as
            # valence, so the shallowest value seen bounds a core rule, the deepest a valence rule)
            levels = levels_from_go([(r.job, r.files["out_go"]) for r in runners.values() if r.has("go")])
            mismatch = check_rules(self.orbitals, levels)
            for k, v in levels.items():
                seen = self._levels_seen.get(k)
                if seen is None:
                    self._levels_seen[k] = v
                else:
                    seen.e_min, seen.e_max = min(seen.e_min, v.e_min), max(seen.e_max, v.e_max)
                    seen.e, seen.star = 0.5 * (seen.e_min + seen.e_max), v.star
        bounds = self._bounds(self._levels_seen if self.orbitals else {})
        if bounds.min_ewidth is not None and bounds.min_ewidth > ewidth:
            # a valence rule asks for a deeper contour than this go used: the dos window must reach
            # E_F - min_ewidth - eth - ediff to see the gap there (Rb 4s valence: min 2.3 Ry)
            need = round(min((bounds.min_ewidth + self.eth + self.ediff) / self.ref + 1e-9, self.ewidth_dos_max), 4)
            for r in runners.values():
                if (r.ewidth_dos_used or r.ewidth_dos) < need - 1e-6:
                    self.logger.info("orbital bounds %s: widening ewidth_dos %.3f -> %.3f", bounds.as_list(), r.ewidth_dos_used or r.ewidth_dos, need)
                    try:
                        r.run_dos(ewidth_dos=need, force=True)
                    except KKRFailedExecutionError as e:
                        self.logger.warning("wider dos window %.3f failed (%s)", need, e)
        while True:
            pairs = []
            for p, r in runners.items():
                if not r.has("dos"):
                    raise GaesError("no dos output for polytyp {} in {}".format(p, r.directory))
                pairs.append((r.job, r.files["out_dos"]))
            energy, curves, _ = dos_curves_from_outputs(pairs)
            dec = decide(self.method, energy, curves, ewidth, dosth=self.dosth, dosth2=self.dosth2, eth=self.eth,
                         ediff=self.ediff, margin=self.margin, dosth2_relax=self.dosth2_relax, min_ewidth=bounds.min_ewidth,
                         max_ewidth=bounds.max_ewidth)
            if not (dec.flag == "fail" and dec.window_limited and self.method == 2):
                break
            current = max(r.ewidth_dos_used or r.ewidth_dos for r in runners.values())
            wider = round(min(current * 1.5, self.ewidth_dos_max), 4)
            if wider <= current + 1e-6 or widened >= 3:
                self.logger.info("window-limited gap at ewidth_dos=%.3f; not widening further", current)
                break
            widened += 1
            self.logger.info("low-DOS region touches the window bottom; widening ewidth_dos %.3f -> %.3f", current, wider)
            try:
                for r in runners.values():
                    r.run_dos(ewidth_dos=wider, force=True)
            except KKRFailedExecutionError as e:
                # specx cannot run the wider window (e.g. reconf): restore the previous dos and judge with it
                self.logger.warning("wider dos window %.3f failed (%s); keeping %.3f", wider, e, current)
                for r in runners.values():
                    r.run_dos(ewidth_dos=current, force=True)
                break
        if mismatch and dec.flag == "old":
            # the go output does not treat the orbitals as asked (a level moved across E_F - ewidth):
            # the bounds were re-derived from the new levels above, so move to the next candidate
            self.logger.info("orbital rules not satisfied at ewidth %.4f: %s; bounds now %s", ewidth, mismatch, bounds.as_list())
            dec.flag = "new" if dec.candidates else "fail"
            dec.ewidth = dec.candidates[0] if dec.candidates else None
            dec.gap_used = None
        flag, next_ew, cands = dec.flag, dec.ewidth if dec.flag == "new" else (ewidth if dec.flag == "old" else None), dec.candidates
        j = Judgement(step=step, iew=iew, ied=ied, ewidth=ewidth, converged=dict(converged),
                      regions=[list(g.as_tuple()) for g in dec.coarse], flag=flag, next_ewidth=next_ew,
                      candidates=list(cands), gap_used=list(dec.gap_used.as_tuple()) if dec.gap_used else None,
                      directories={p: r.directory for p, r in runners.items()},
                      ewidth_dos=next(iter(runners.values())).ewidth_dos_used if runners else None,
                      method=self.method, fine_regions=[list(f.as_tuple()) for f in dec.fine], relaxed=dec.relaxed,
                      dosth2_used=dec.dosth2_used, window_limited=dec.window_limited,
                      orbital_levels=levels_as_dict(levels), orbital_bounds=bounds.as_list(), orbital_mismatch=mismatch)
        for p, r in runners.items():
            if r.ref_effective is not None and abs(r.ref_effective - self.ref) > 0.05:
                self.logger.warning("polytyp %s: dos window gives ref=%.3f but Gaes(ref=%.3f); the window "
                                    "[%.3f, %.3f] may not reach E_F - ewidth_go - eth - ediff = %.3f",
                                    p, r.ref_effective, self.ref, energy.min(), energy.max(),
                                    -ewidth - self.eth - self.ediff)
        self.logger.info("judge %s iew=%d ied=%d ewidth=%.4f converged=%s regions=%s bounds=%s -> %s %s",
                         step, iew, ied, ewidth, converged, j.regions, bounds.as_list(), flag, next_ew)
        return j

    # ---------------------------------------------------------------- steps
    def _step1(self, key, params, iew, ewidth):
        """go -> dos (-> j) once per polytyp with maxitr_init and pmix_init."""
        runners, converged = {}, {}
        for polytyp, param_go in params.items():
            point = RunPoint(key, iew, ewidth, 0, self.edelt_init, polytyp, 0, self.pmix_init)
            dic = deepcopy(param_go)
            dic["maxitr"] = self.maxitr_init
            r = self._runner(point, dic)
            converged[polytyp] = r.run_all(with_j=self.with_j, dos_always=True,
                                           min_ewidth_dos=self._min_ewidth_dos(ewidth))
            runners[polytyp] = r
        return runners, converged

    def _chain(self, key, param_go, iew, ewidth, ied, edelt, polytyp, pmix, start_dir, ipm):
        """up to max_pm_iter go runs with one pmix, each continued from the previous run while
        it still converges. start_dir=None starts from a fresh potential.
        Returns (converged, last runner, next ipm, last directory that was converging or None)."""
        prev_dir = start_dir
        last = None
        for _it in range(self.max_pm_iter):
            point = RunPoint(key, iew, ewidth, ied, edelt, polytyp, ipm, pmix)
            dic = deepcopy(param_go)
            dic["maxitr"] = self.maxitr_pm
            r = self._runner(point, dic)
            copy_from = prev_dir if (prev_dir is not None and r.directory != prev_dir) else None
            if prev_dir is None and r.ran.get("go") is None:
                # fresh start: make sure no old potential is picked up by record=2nd
                pot = os.path.join(r.directory, r.files["potential"])
                if os.path.isfile(pot) and not (r._same_inputcard(r.files["inputcard_go"], r._card_text(r._param("go")))
                                                and r.has("go")):
                    os.remove(pot)
            conv = r.run_go(copy_potential_from=copy_from)
            last = r
            ipm += 1
            if conv:
                return True, r, ipm, r.directory
            job = r.job
            if is_converging(job.get_rms_error(r.files["out_go"]), job.get_moment_history(r.files["out_go"])):
                prev_dir = r.directory
                continue
            return False, r, ipm, None
        return False, last, ipm, prev_dir

    def _step2_polytyp(self, key, param_go, iew, ewidth, ied, edelt, polytyp, start_dir):
        """edelt fixed; for each pmix: continue from the previous potential; if that does not
        converge, retry the same pmix from a fresh potential (fresh_retry). Returns (converged, last runner)."""
        ipm = 0
        prev_dir = start_dir
        last = None
        for pmix in self.pmix_steps:
            sources = [("previous", prev_dir)]
            if self.fresh_retry:
                sources.append(("fresh", None))
            for source, start in sources:
                self.logger.info("key %s %s: edelt=%g pmix=%g from %s potential", key, polytyp, edelt, pmix, source)
                conv, r, ipm, good_dir = self._chain(key, param_go, iew, ewidth, ied, edelt, polytyp, pmix, start, ipm)
                last = r
                if conv:
                    r.run_dos(min_ewidth_dos=self._min_ewidth_dos(ewidth))
                    if self.with_j:
                        r.run_j()
                    return True, r
                if good_dir is not None:
                    prev_dir = good_dir     # the chain was still converging: continue from it next
                    break                   # no fresh retry needed for a chain that only ran out of iterations
        # unconverged: still need a DOS to judge the ewidth
        last.run_dos(min_ewidth_dos=self._min_ewidth_dos(ewidth))
        if self.compat and self.with_j:
            last.run_j()
        return False, last

    # ---------------------------------------------------------------- driver
    def run(self, key, params, save=True):
        """tune ewidth for one key. params: polytyp -> go parameter dict (make_inputcard format).

        Returns:
            KeyResult (also written to <prefix>/key_<key>.json when save=True).
        """
        res = KeyResult(key=key, polytyps=list(params), parameters=self.parameters())
        try:
            self._run(key, params, res)
        except (GaesError, KKRFailedExecutionError, OSError, ValueError) as e:
            res.status = "error"
            res.message = "{}: {}".format(type(e).__name__, e)
            self.logger.error("key %s: %s", key, res.message)
        if save:
            os.makedirs(self.layout.prefix, exist_ok=True)
            res.save(os.path.join(self.layout.prefix, "key_{}.json".format(key)))
        return res

    def _finish(self, res, status, ewidth, runners, converged, judgement):
        res.status = status
        if status == "ewidth_fail" and judgement is not None and judgement.orbital_bounds and not res.message:
            res.message = "no gap candidate inside [min_ewidth, max_ewidth] = {}{}".format(
                judgement.orbital_bounds, " (orbital rules {})".format([str(r) for r in self.orbitals]) if self.orbitals else "")
        res.ewidth_final = ewidth
        res.converged = dict(converged)
        res.final = {p: r.directory for p, r in runners.items()}
        res.gap_used = judgement.gap_used if judgement is not None else None
        res.results = {p: r.result(dosth=self.dosth) for p, r in runners.items()}

    def _next_ewidth(self, candidates, tried):
        for c in candidates:
            if round(c, 6) not in tried:
                return c
        return None

    def _run(self, key, params, res):
        ewidth = self.ewidth_init
        self._levels_seen = {}
        if self.orbitals:
            # step 0: the converged table gives a first range; the go outputs replace it at every judgement
            levels0 = levels_for_step0({r.element for r in self.orbitals}, self.ef_assumed)
            bounds0 = self._bounds(levels0, strict=False)
            ewidth = initial_ewidth(self.ewidth_init, bounds0, self.ediff)
            res.parameters["orbital_bounds_step0"] = bounds0.as_list()
            res.parameters["orbital_levels_step0"] = {k: [round(v.e, 4), v.source] for k, v in levels0.items()}
            self.logger.info("key %s: orbital rules %s -> step-0 bounds %s from %s; first ewidth %.4f", key,
                             [str(r) for r in self.orbitals], bounds0.as_list(),
                             {k: v.source for k, v in levels0.items()}, ewidth)
        tried = set()
        iew = 0
        runners, converged, judgement = {}, {}, None
        while iew < self.max_ew:
            # ---- STEP1: rough ewidth
            runners, converged = self._step1(key, params, iew, ewidth)
            tried.add(round(ewidth, 6))
            res.ewidth_tried.append(ewidth)
            judgement = self._judge("step1", iew, 0, ewidth, runners, converged)
            res.judgements.append(judgement)
            if judgement.flag == "fail" and (all(converged.values()) or not self.tighten_before_fail):
                return self._finish(res, "ewidth_fail", ewidth, runners, converged, judgement)
            if judgement.flag == "new":
                nxt = self._next_ewidth(judgement.candidates, tried)
                if nxt is None:
                    return self._finish(res, "not_converged", ewidth, runners, converged, judgement)
                if self.compat:
                    # 2019 behaviour: iew is not advanced, so the same directory (and result) is reused
                    ewidth = nxt
                    continue
                ewidth = nxt
                iew += 1
                continue
            if all(converged.values()):
                return self._finish(res, "finished", ewidth, runners, converged, judgement)
            # ---- STEP2: tighten the SCF (edelt rounds; then the same with a larger bzqlty)
            start_dirs = {p: r.directory for p, r in runners.items()}
            go_new = False
            n_ed = len(self.edelt_steps)
            for bz_round, bz in enumerate((None,) + self.bzqlty_steps):
                if bz is not None and all(converged.values()):
                    break
                for i_ed, edelt in enumerate(self.edelt_steps):
                    ied = bz_round * n_ed + i_ed
                    for polytyp, param_go in params.items():
                        if converged.get(polytyp):
                            continue   # already converged at a previous stage; keep it
                        dic = deepcopy(param_go)
                        if bz is not None:
                            dic["bzqlty"] = self._bzqlty(dic.get("bzqlty", 6), bz)
                            self.logger.info("key %s %s: bzqlty %s -> %s", key, polytyp, param_go.get("bzqlty"), dic["bzqlty"])
                        conv, r = self._step2_polytyp(key, dic, iew, ewidth, ied, edelt, polytyp, start_dirs[polytyp])
                        converged[polytyp] = conv
                        runners[polytyp] = r
                        start_dirs[polytyp] = r.directory
                    judgement = self._judge("step2", iew, ied, ewidth, runners, converged)
                    res.judgements.append(judgement)
                    if judgement.flag == "fail" and (all(converged.values()) or not self.tighten_before_fail):
                        return self._finish(res, "ewidth_fail", ewidth, runners, converged, judgement)
                    if judgement.flag == "new":
                        go_new = True
                        break
                    if judgement.flag == "old" and all(converged.values()):
                        return self._finish(res, "finished", ewidth, runners, converged, judgement)
                if go_new:
                    break
            if judgement.flag == "fail":
                return self._finish(res, "ewidth_fail", ewidth, runners, converged, judgement)
            nxt = self._next_ewidth(judgement.candidates, tried)
            if nxt is None:
                return self._finish(res, "not_converged", ewidth, runners, converged, judgement)
            ewidth = nxt
            iew += 1
            if not go_new:
                self.logger.info("key %s: SCF not converged with ewidth %.4f; trying alternative %.4f",
                                 key, res.ewidth_tried[-1], ewidth)
        self._finish(res, "ewidth_exhausted", ewidth, runners, converged, judgement)


def collect(prefix):
    """all key_<key>.json under prefix as a DataFrame (one row per key and polytyp)."""
    import glob
    import pandas as pd
    rows = []
    for path in sorted(glob.glob(os.path.join(prefix, "key_*.json"))):
        with open(path) as f:
            d = json.load(f)
        for polytyp in d.get("polytyps", []):
            row = {"key": d["key"], "polytyp": polytyp, "status": d["status"], "ewidth_final": d.get("ewidth_final"),
                   "n_ewidth_tried": len(d.get("ewidth_tried", [])), "gap_used": d.get("gap_used"),
                   "directory": d.get("final", {}).get(polytyp), "converged": d.get("converged", {}).get(polytyp)}
            r = d.get("results", {}).get(polytyp, {})
            for k in ("total_energy_Ry", "total_moment", "a_bohr", "volume_bohr3", "n_iter", "last_err",
                      "Tc_K", "ewidth_dos", "low_dos_regions"):
                row[k] = r.get(k)
            rows.append(row)
    return pd.DataFrame(rows)


def collect_legacy(prefix, akaikkr_exe="", dosth=2e-2):
    """rows of the 2019 RUN directories (Layout version 1); ewidth etc. are read from the outputs."""
    import pandas as pd
    layout = Layout(prefix, version=1)
    rows = []
    for p in layout.find():
        r = KkrRunner(akaikkr_exe, layout.path(p), {"go": "go"}, compat=True)
        row = r.result(dosth=dosth, save_csv=False)
        row.update({"key": p.key, "polytyp": p.polytyp, "iew_dir": p.iew, "ied_dir": p.ied, "ipm_dir": p.ipm})
        rows.append(row)
    return pd.DataFrame(rows)
