# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""kkr-gaes: command line of the ewidth tuning scheme."""
import argparse
import json
import logging
import os
import sys

from ..AkaiKkr import AkaikkrJob
from .composition import SiteComposition, make_single_site_param, check_type_name
from .ewidth import decide, ETH, EDIFF, MARGIN
from .gap import dos_curves_from_outputs
from .layout import Layout
from .legacy_hea import heakey_to_composition, composition_to_heakey, load_heakeylist
from .scheme import Gaes, collect, collect_legacy


def _add_scheme_options(p):
    p.add_argument("--exe", required=True, help="path of specx")
    p.add_argument("--prefix", default="RUN")
    p.add_argument("--method", type=int, default=2, choices=(1, 2), help="1: single threshold (2019), 2: two thresholds")
    p.add_argument("--dosth", type=float, default=2e-2, help="threshold (Method 1) / coarse threshold (Method 2)")
    p.add_argument("--dosth2", type=float, default=1e-3, help="fine threshold of Method 2")
    p.add_argument("--min-ewidth", type=float, default=1.0, help="drop candidates shallower than this (Ry)")
    p.add_argument("--max-ewidth", type=float, default=2.0, help="drop candidates deeper than this (Ry)")
    p.add_argument("--ewidth-init", type=float, default=1.2)
    p.add_argument("--ewidth-dos", type=float, default=3.0)
    p.add_argument("--ref", type=float, default=0.75, help="cemesr ref of the build (akaikkr 0.75, akaikkr_cnd 0.5)")
    p.add_argument("--max-ew", type=int, default=10)
    p.add_argument("--max-pm-iter", type=int, default=20)
    p.add_argument("--maxitr-init", type=int, default=500)
    p.add_argument("--edelt-init", type=float, default=1e-4, help="edelt of STEP1")
    p.add_argument("--edelt-dos", type=float, default=1e-4, help="edelt of every dos used for the gap judgement")
    p.add_argument("--edelt-steps", type=float, nargs="+", default=[1e-2, 1e-3, 1e-4], help="edelt of STEP2, large to small")
    p.add_argument("--no-j", action="store_true")
    p.add_argument("--compat", action="store_true", help="reproduce the 2019 behaviour (layout v1, dosth 2e-2)")
    p.add_argument("--threads", type=int, default=None, help="OMP_NUM_THREADS of specx")
    p.add_argument("--polytyp", nargs="+", default=["bcc", "fcc"])
    p.add_argument("--lattice", default="expr", help="expr | mjw | <a in bohr>")
    p.add_argument("--type-name", default=None, help="AkaiKKR type name (default: from the composition; 2019 used HEA)")


def _gaes(args):
    if args.threads:
        os.environ["OMP_NUM_THREADS"] = str(args.threads)
    os.makedirs(args.prefix, exist_ok=True)
    logging.basicConfig(filename=os.path.join(args.prefix, "kkr-gaes.log"), level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    layout = Layout(args.prefix, version=1 if args.compat else 2)
    return Gaes(args.exe, layout, ewidth_init=args.ewidth_init, ewidth_dos=args.ewidth_dos, ref=args.ref,
                method=args.method, dosth=args.dosth, dosth2=args.dosth2, min_ewidth=args.min_ewidth, max_ewidth=args.max_ewidth,
                max_ew=args.max_ew, max_pm_iter=args.max_pm_iter, maxitr_init=args.maxitr_init,
                edelt_init=args.edelt_init, edelt_dos=args.edelt_dos, edelt_steps=tuple(args.edelt_steps),
                with_j=not args.no_j, compat=args.compat)


def _lattice(s):
    return s if s in ("expr", "mjw") else float(s)


def _site_params(comp, args):
    return {pt: make_single_site_param(comp, pt, lattice=_lattice(args.lattice), type_name=args.type_name)
            for pt in args.polytyp}


def cmd_run(args):
    with open(args.input) as f:
        items = json.load(f)
    g = _gaes(args)
    for key, params in items.items():
        r = g.run(key, params)
        print(key, r.status, r.ewidth_final, r.gap_used)


def cmd_site(args):
    g = _gaes(args)
    for spec in args.comp:
        comp = SiteComposition.from_type_name(spec)
        key = comp.key()
        r = g.run(key, _site_params(comp, args))
        print(key, r.status, r.ewidth_final, r.gap_used)


def cmd_hea(args):
    g = _gaes(args)
    keys = load_heakeylist(args.keys)[args.start:args.stop]
    for key in keys:
        comp = heakey_to_composition(key)
        r = g.run(key, _site_params(comp, args))
        print(key, r.status, r.ewidth_final, r.gap_used)


def cmd_check(args):
    pairs = []
    for path in args.dos:
        d, f = os.path.split(path)
        pairs.append((AkaikkrJob(d or "."), f))
    energy, curves, labels = dos_curves_from_outputs(pairs)
    dec = decide(args.method, energy, curves, args.ewidth, dosth=args.dosth, dosth2=args.dosth2, eth=args.eth,
                 ediff=args.ediff, margin=MARGIN)
    print("window: [{:.4f}, {:.4f}] Ry, {} points; method {}".format(energy.min(), energy.max(), len(energy), args.method))
    print("gap regions (DOS < {:g}{}):".format(args.dosth, ", wider than eth" if args.method == 2 else ""))
    for g in dec.coarse:
        print("  [{:.4f}, {:.4f}]  width {:.4f}".format(g.e1, g.e2, g.width))
    if args.method == 2:
        print("fine sub-regions (DOS < {:g}{}):".format(dec.dosth2_used or args.dosth2, ", relaxed" if dec.relaxed else ""))
        for f in dec.fine:
            print("  [{:.4f}, {:.4f}]".format(f.e1, f.e2))
    print("ewidth {:.4f}: {}".format(args.ewidth, dec.flag),
          "" if dec.flag == "fail" else "next {:.4f}".format(dec.ewidth) if dec.flag == "new" else "keep")
    print("candidates:", ", ".join("{:.4f}".format(c) for c in dec.candidates))


def cmd_collect(args):
    df = collect_legacy(args.prefix) if args.legacy else collect(args.prefix)
    df.to_csv(args.output, index=False)
    print("wrote", args.output, len(df), "rows")


def cmd_comp(args):
    for s in args.spec:
        if s.isdigit():
            comp = heakey_to_composition(s)
        else:
            comp = SiteComposition.from_type_name(s)
        line = {"elements": comp.elements, "fractions": comp.fractions, "z": comp.z, "conc": comp.conc,
                "type_name": comp.type_name(max_len=10 ** 6), "key": comp.key()}
        try:
            check_type_name(line["type_name"])
            line["type_name_ok"] = True
        except Exception as e:  # noqa: BLE001
            line["type_name_ok"] = str(e)
        if comp.is_equiatomic:
            line["heakey"] = composition_to_heakey(comp)
        print(json.dumps(line, ensure_ascii=False))


def main(argv=None):
    p = argparse.ArgumentParser(prog="kkr-gaes", description="GAES: gap-anchored ewidth search for AkaiKKR")
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("run", help="tune ewidth for {key: {polytyp: param_go}} in a JSON file")
    _add_scheme_options(s); s.add_argument("--input", required=True); s.set_defaults(func=cmd_run)
    s = sub.add_parser("site", help="single-site CPA from compositions (Rh0.5Pt0.5, AlSiScTi, ...)")
    _add_scheme_options(s); s.add_argument("--comp", nargs="+", required=True); s.set_defaults(func=cmd_site)
    s = sub.add_parser("hea", help="2019 heakey list (heakeylist0.csv)")
    _add_scheme_options(s); s.add_argument("--keys", required=True)
    s.add_argument("--start", type=int, default=0); s.add_argument("--stop", type=int, default=None)
    s.set_defaults(func=cmd_hea)
    s = sub.add_parser("check", help="judge an ewidth against existing dos outputs (no specx run)")
    s.add_argument("--dos", nargs="+", required=True); s.add_argument("--ewidth", type=float, required=True)
    s.add_argument("--method", type=int, default=2, choices=(1, 2))
    s.add_argument("--dosth", type=float, default=2e-2); s.add_argument("--dosth2", type=float, default=1e-3)
    s.add_argument("--eth", type=float, default=ETH)
    s.add_argument("--ediff", type=float, default=EDIFF); s.set_defaults(func=cmd_check)
    s = sub.add_parser("collect", help="collect key_*.json (or the 2019 RUN with --legacy) into a CSV")
    s.add_argument("--prefix", default="RUN"); s.add_argument("--legacy", action="store_true")
    s.add_argument("-o", "--output", default="result.csv"); s.set_defaults(func=cmd_collect)
    s = sub.add_parser("comp", help="composition / type name / heakey conversion and the 40-character check")
    s.add_argument("spec", nargs="+"); s.set_defaults(func=cmd_comp)
    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
