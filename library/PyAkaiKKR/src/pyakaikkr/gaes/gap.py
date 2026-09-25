# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""band-gap regions of a DOS: contiguous energy-mesh regions with DOS(E) < threshold.

The criterion is the DOS value only. DOS derivatives are not used because the
AkaiKKR DOS is sometimes jagged (docs/ewidth_tuning_scheme.md section 0.2).
"""
from dataclasses import dataclass

import numpy as np

from ..Error import GaesError


@dataclass(frozen=True)
class GapRegion:
    """contiguous region [e1, e2] (Ry, E - EF) where every curve is below the threshold.

    i1, i2 are the mesh indices of e1 and e2. e2 is the first point after the run of
    low-DOS points (the same convention as the 2019 script), or the last mesh point
    when the run reaches the end of the mesh.
    """
    e1: float
    e2: float
    i1: int
    i2: int

    @property
    def width(self):
        return self.e2 - self.e1

    def contains(self, e):
        return self.e1 < e < self.e2

    def as_tuple(self):
        return (self.e1, self.e2)


def gap_regions(energy, curves, dosth=1e-3, labels=None):
    """contiguous mesh regions where all curves are below dosth.

    Args:
        energy (array-like): energy mesh, E - EF in Ry (increasing).
        curves (Sequence[array-like]): DOS curves on the same mesh (total DOS of each
            polytyp, PDOS of each component, ...). The AND of all curves is taken.
        dosth (float, optional): threshold. Defaults to 1e-3.
        labels (Sequence[str], optional): names of the curves for error messages.

    Raises:
        GaesError: no curve, or a curve whose length differs from the mesh.

    Returns:
        list[GapRegion]: regions in increasing energy order.
    """
    e = np.asarray(energy, dtype=float)
    if len(curves) == 0:
        raise GaesError("gap_regions: no DOS curve")
    mask = np.ones(len(e), dtype=bool)
    for k, c in enumerate(curves):
        c = np.asarray(c, dtype=float)
        if c.shape != e.shape:
            name = labels[k] if labels and k < len(labels) else str(k)
            raise GaesError("gap_regions: curve {} has {} points, mesh has {}".format(
                name, c.shape, e.shape))
        mask &= c < dosth
    regions = []
    start = None
    for i, low in enumerate(mask):
        if low and start is None:
            start = i
        elif not low and start is not None:
            regions.append(GapRegion(float(e[start]), float(e[i]), start, i))
            start = None
    if start is not None:
        i = len(e) - 1
        regions.append(GapRegion(float(e[start]), float(e[i]), start, i))
    return regions


def same_mesh(meshes, tol=1e-6):
    """True if all energy meshes are equal within tol."""
    m0 = np.asarray(meshes[0], dtype=float)
    for m in meshes[1:]:
        m = np.asarray(m, dtype=float)
        if m.shape != m0.shape or np.max(np.abs(m - m0)) > tol:
            return False
    return True


def dos_curves_from_outputs(jobs_and_files, spin_sum=True):
    """total DOS curves (spin sum) of several dos outputs on a common mesh.

    Args:
        jobs_and_files (Sequence[tuple[AkaikkrJob, str]]): (job, dos output filename).
        spin_sum (bool, optional): add up and down. Defaults to True.

    Raises:
        GaesError: the meshes differ.

    Returns:
        (np.ndarray, list[np.ndarray], list[str]): energy, curves, labels (path_dir/outfile).
    """
    energies, curves, labels = [], [], []
    for job, outfile in jobs_and_files:
        e, block = job.get_dos_as_list(outfile)
        d = np.asarray(block[0], dtype=float)
        if spin_sum and len(block) > 1:
            d = d + np.asarray(block[1], dtype=float)
        energies.append(np.asarray(e, dtype=float))
        curves.append(d)
        labels.append("{}/{}".format(job.path_dir, outfile))
    if not same_mesh(energies):
        raise GaesError("dos energy meshes differ: {}".format(labels))
    return energies[0], curves, labels


def pdos_curves_from_output(job, outfile, l_sum=True, spin_sum=True):
    """PDOS curves of every component of one dos output (entry point of the PDOS extension).

    Returns:
        (np.ndarray, list[np.ndarray], list[str]): energy, one curve per component
        (summed over l and spin by default), component short names.
    """
    typeofsite = job.get_type_of_site(outfile)
    names = [s for site in typeofsite for s in site["comp_shortname"]]
    energy, pdos_block = job.get_pdos_as_list(outfile, output_format="spin_separation")
    energy = np.asarray(energy, dtype=float)
    curves = []
    for icmp in range(len(pdos_block[0])):
        c = np.asarray(pdos_block[0][icmp], dtype=float)
        if spin_sum and len(pdos_block) > 1:
            c = c + np.asarray(pdos_block[1][icmp], dtype=float)
        curves.append(c.sum(axis=1) if l_sum else c)
    return energy, curves, names[:len(curves)]
