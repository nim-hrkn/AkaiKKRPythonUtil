# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""GAES figure helpers: re-exported from pyakaikkr.plot (the array-based drawing module shared with
aiida-akaikkr). Kept so that ``from pyakaikkr.gaes.plot import draw_gaes_dos`` keeps working."""
from ..plot import (DEFAULT_STYLE, shade_regions, shade_ewidth_bounds, draw_ewidth, draw_levels, draw_thresholds,
                    plot_gaes_dos, gaes_legend_text)

draw_gaes_dos = plot_gaes_dos
legend_text = gaes_legend_text

__all__ = ["DEFAULT_STYLE", "draw_gaes_dos", "plot_gaes_dos", "shade_regions", "shade_ewidth_bounds", "draw_ewidth",
           "draw_levels", "draw_thresholds", "legend_text", "gaes_legend_text"]
