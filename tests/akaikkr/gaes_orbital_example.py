# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""GAES example 2: ask for one orbital to be valence (inside the SCF contour) or core (below it).

usage: python gaes_orbital_example.py <program_path> [--comp SeMnFeCo] [--orbital Se4s] [--threads 20]

Runs GAES twice on the same system: "<orbital>=valence" (alias occupied) and "<orbital>=core"
(alias unoccupied). The rule turns the core level of that orbital, read from every go output, into
the ewidth range: valence -> min_ewidth = |E - E_F| + ediff, core -> max_ewidth = |E - E_F| - ediff
(docs/ewidth_tuning_scheme.md section 15). Runs go to RUN_gaes_orbital/<comp>_<polytyp>_<rule>/,
the figure to gaes_orbital_<comp>_<polytyp>_<orbital>.png (left valence, right core).
"""
import argparse
import os

from pyakaikkr.gaes import Gaes, Layout, SiteComposition, make_single_site_param
from gaes_example_plot import plot_key_results


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("program_path", help="directory containing akaikkr/specx")
    p.add_argument("--comp", default="SeMnFeCo")
    p.add_argument("--polytyp", default="fcc")
    p.add_argument("--orbital", default="Se4s", help="element + orbital, e.g. Se4s, Rb4p, Bi6s")
    p.add_argument("--ewidth-init", type=float, default=1.2)
    p.add_argument("--eth", type=float, default=0.3, help="minimum gap width (Ry); core rules on a shallow level may need 0.2")
    p.add_argument("--build", default="akaikkr")
    p.add_argument("--threads", type=int, default=None)
    args = p.parse_args()
    if args.threads:
        os.environ["OMP_NUM_THREADS"] = str(args.threads)
    exe = os.path.join(args.program_path, args.build, "specx")
    ref = 0.75 if args.build == "akaikkr" else 0.5

    comp = SiteComposition.from_type_name(args.comp)
    key = "{}_{}".format(args.comp, args.polytyp)
    params = {args.polytyp: make_single_site_param(comp, args.polytyp, type_name="HEA")}
    results = {}
    for role in ("valence", "core"):
        rule = "{}={}".format(args.orbital, role)
        g = Gaes(exe, Layout(os.path.join("RUN_gaes_orbital", "{}_{}".format(key, rule.replace("=", "-"))), version=2),
                 ewidth_init=args.ewidth_init, ref=ref, method=2, dosth=2e-2, dosth2=1e-3, eth=args.eth,
                 orbitals=[rule],                     # the range [min_ewidth, max_ewidth] comes from the rule
                 edelt_init=1e-4, edelt_dos=1e-4, with_j=False)
        res = g.run(key, params)
        results[rule] = res
        print("==", rule)
        print("status  :", res.status, " ewidth:", res.ewidth_final, " tried:", [round(t, 4) for t in res.ewidth_tried])
        print("step-0 bounds from the level table:", res.parameters.get("orbital_bounds_step0"))
        for j in res.judgements:
            lv = j.orbital_levels.get(args.orbital)
            print("  ewidth %.4f %-4s converged=%s bounds=%s %s=%s mismatch=%s fine=%s" % (
                j.ewidth, j.flag, j.converged, [None if b is None else round(b, 3) for b in j.orbital_bounds],
                args.orbital, None if lv is None else (round(lv[0], 3), "valence*" if lv[1] else "core"),
                [m["rule"] for m in j.orbital_mismatch], [[round(a, 3), round(b, 3)] for a, b in j.fine_regions]))
            for r in j.reasons:
                print("     note:", r)
        if res.message:
            print("message :", res.message)
    png = plot_key_results(results, "gaes_orbital_{}_{}.png".format(key, args.orbital),
                           title="GAES with an orbital rule: {} {}".format(key, args.orbital), highlight=[args.orbital])
    print("figure  :", png)


if __name__ == "__main__":
    main()
