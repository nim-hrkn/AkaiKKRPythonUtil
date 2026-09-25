# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""GAES example 1: decide the ewidth of go by the gap search only (no orbital rule).

usage: python gaes_ewidth_example.py <program_path> [--comp SeMnFeCo] [--polytyp fcc] [--ewidth-init 1.2]
                                     [--min-ewidth 1.0] [--max-ewidth 2.0] [--threads 20]

<program_path> is the directory that contains akaikkr/specx (as for testrun.py). The run goes to
RUN_gaes_ewidth/<comp>_<polytyp>/ (key_<comp>.json + one directory per go), the figure to
gaes_ewidth_<comp>_<polytyp>.png. See docs/gaes_usage.md and docs/ewidth_tuning_scheme.md.
"""
import argparse
import os

from pyakaikkr.gaes import Gaes, Layout, SiteComposition, make_single_site_param
from gaes_example_plot import plot_key_result


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("program_path", help="directory containing akaikkr/specx")
    p.add_argument("--comp", default="SeMnFeCo", help="single-site CPA composition (equiatomic elements, or Rh0.5Pt0.5)")
    p.add_argument("--polytyp", default="fcc")
    p.add_argument("--ewidth-init", type=float, default=1.2)
    p.add_argument("--min-ewidth", type=float, default=1.0)
    p.add_argument("--max-ewidth", type=float, default=2.0)
    p.add_argument("--build", default="akaikkr", help="akaikkr | akaikkr_cpa2021v01 (ref 0.5) | akaikkr_cnd (ref 0.5)")
    p.add_argument("--threads", type=int, default=None)
    args = p.parse_args()
    if args.threads:
        os.environ["OMP_NUM_THREADS"] = str(args.threads)
    exe = os.path.join(args.program_path, args.build, "specx")
    ref = 0.75 if args.build == "akaikkr" else 0.5          # dos window of the build (cemesr ref)

    comp = SiteComposition.from_type_name(args.comp)
    key = "{}_{}".format(args.comp, args.polytyp)
    params = {args.polytyp: make_single_site_param(comp, args.polytyp, type_name="HEA")}
    g = Gaes(exe, Layout(os.path.join("RUN_gaes_ewidth", key), version=2),
             ewidth_init=args.ewidth_init, ref=ref, method=2, dosth=2e-2, dosth2=1e-3,
             min_ewidth=args.min_ewidth, max_ewidth=args.max_ewidth,
             edelt_init=1e-4, edelt_dos=1e-4,        # go and judgement dos both at 1e-4 (2019 setting)
             with_j=False)
    res = g.run(key, params)

    print("status  :", res.status)
    print("ewidth  :", res.ewidth_final, " tried:", [round(t, 4) for t in res.ewidth_tried])
    print("gap used:", res.gap_used)
    for j in res.judgements:
        print("  ewidth %.4f %-4s converged=%s coarse=%s fine=%s candidates=%s" % (
            j.ewidth, j.flag, j.converged, [[round(a, 3), round(b, 3)] for a, b in j.regions],
            [[round(a, 3), round(b, 3)] for a, b in j.fine_regions], [round(c, 4) for c in j.candidates]))
        for r in j.reasons:
            print("     note:", r)
    if res.message:
        print("message :", res.message)
    png = plot_key_result(res, "gaes_ewidth_{}.png".format(key), title="GAES (ewidth only): " + key)
    print("figure  :", png)


if __name__ == "__main__":
    main()
