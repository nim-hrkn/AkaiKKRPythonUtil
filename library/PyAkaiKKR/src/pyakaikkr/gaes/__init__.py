# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""GAES (Gap-Anchored Ewidth Search): automatic ewidth tuning for AkaiKKR.

The method was used for the HEA screening of T. Fukushima, H. Akai, T. Chikyow, H. Kino,
Phys. Rev. Materials 6, 023802 (2022), doi:10.1103/PhysRevMaterials.6.023802, but is itself
unpublished; docs/ewidth_tuning_scheme.md is its first description. Not imported by ``import pyakaikkr``; use
``from pyakaikkr.gaes import ...``.
"""
from .gap import GapRegion, gap_regions, dos_curves_from_outputs, pdos_curves_from_output
from .ewidth import check_ewidth, choose_ewidth, ewidth_candidates, choose_ewidth2, decide, Decision
from .convergence import is_converging
from .layout import Layout, RunPoint
from .runner import KkrRunner, output_finished
from .composition import SiteComposition, make_single_site_param, check_type_name, TYPE_NAME_MAX_LEN
from .legacy_hea import heakey_to_composition, composition_to_heakey, load_heakeylist
from .scheme import Gaes, KeyResult, Judgement, collect, collect_legacy

__all__ = ["GapRegion", "gap_regions", "dos_curves_from_outputs", "pdos_curves_from_output",
           "check_ewidth", "choose_ewidth", "ewidth_candidates", "choose_ewidth2", "decide", "Decision", "is_converging",
           "Layout", "RunPoint", "KkrRunner", "output_finished",
           "SiteComposition", "make_single_site_param", "check_type_name", "TYPE_NAME_MAX_LEN",
           "heakey_to_composition", "composition_to_heakey", "load_heakeylist",
           "Gaes", "KeyResult", "Judgement", "collect", "collect_legacy"]
