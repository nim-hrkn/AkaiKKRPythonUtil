# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""is the SCF still heading to convergence? (linear-regression test of the 2019 script)"""
import numpy as np


def _fit(y):
    """(r2 of the standardized linear fit, MAE of the raw fit) of y against its index."""
    y = np.asarray(y, dtype=float)
    x = np.arange(len(y), dtype=float)
    sy = y.std()
    if len(y) < 2 or sy == 0.0 or x.std() == 0.0:
        return 0.0, 0.0
    xs = (x - x.mean()) / x.std()
    ys = (y - y.mean()) / sy
    a, b = np.polyfit(xs, ys, 1)
    yp = a * xs + b
    ss_res = float(((ys - yp) ** 2).sum())
    ss_tot = float(((ys - ys.mean()) ** 2).sum())
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    mae = float(np.abs((ys - yp) * sy).mean())
    return r2, mae


def is_converging(err_history, moment_history, last=100, r2_th=0.80, mae_err_th=5e-4,
                  mae_moment_th=5e-5, mae_err_loose_th=6e-3):
    """True if the last `last` points of err (log10 rms) and moment still look like a
    converging run: r2_err > r2_th, or mae_err < mae_err_th, or
    (mae_moment < mae_moment_th and mae_err < mae_err_loose_th).

    Fewer than 3 points give False.
    """
    err = list(err_history)[-last:]
    mom = list(moment_history)[-last:]
    if len(err) < 3 or len(mom) < 3:
        return False
    r2_err, mae_err = _fit(err)
    _r2_mom, mae_mom = _fit(mom)
    if r2_err > r2_th:
        return True
    if mae_err < mae_err_th:
        return True
    if mae_mom < mae_moment_th and mae_err < mae_err_loose_th:
        return True
    return False
