# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""ase.Atoms -> AkaiKKR lattice and atomic positions."""
from __future__ import annotations

import numpy as np
from ase import Atoms
from ase.calculators.calculator import InputError

from ..Unit import Unit
from ..Error import KKRStructureMismatchError
from .occupancy import atoms_to_kkr_types

_BOHR = Unit().length_au2ang  # the same conversion factor as Cif2Kkr


def atoms_to_kkr_lattice(atoms: Atoms) -> dict:
    """make the lattice section of AkaiKKR with brvtyp=aux.

    a = |cell[0]| in bohr, r1,r2,r3 = cell / |cell[0]| (in units of a).

    Args:
        atoms (Atoms): periodic structure

    Raises:
        InputError: non-periodic, left-handed or degenerate cell

    Returns:
        dict: brvtyp, a, c/a, b/a, alpha, beta, gamma, r1, r2, r3
    """
    if len(atoms) == 0:
        raise InputError("no atoms.")
    if not np.all(atoms.pbc):
        raise InputError("AkaiKKR needs a 3D periodic structure (atoms.pbc must be all True).")
    cell = np.array(atoms.cell, dtype=float)
    if np.linalg.det(cell) <= 1e-12:
        raise InputError("the cell must be right-handed and non-degenerate.")
    a_ang = float(np.linalg.norm(cell[0]))
    if a_ang < 1e-8:
        raise InputError("cell[0] is zero.")
    r = cell / a_ang
    return {"brvtyp": "aux", "a": a_ang / _BOHR,
            "c/a": 1.0, "b/a": 1.0, "alpha": 90.0, "beta": 90.0, "gamma": 90.0,
            "r1": r[0].tolist(), "r2": r[1].tolist(), "r3": r[2].tolist()}


def atoms_to_kkr_atmicx(atoms: Atoms, type_of_atom: list, type_names: list | None = None) -> dict:
    """make the atom section of AkaiKKR.

    Fractional coordinates with respect to r1, r2, r3 (a/b/c suffix), as Cif2Kkr.
    If type_names is given, the atoms are grouped by type in that order (stable),
    as Cif2Kkr does. The automatic muffin-tin radii (rmt=0) depend on this order.

    Args:
        atoms (Atoms): structure
        type_of_atom (list): type name of each atom
        type_names (list, optional): type names in the type order. Defaults to None (atom order).

    Returns:
        dict: natm, atmicx, atom_order (index of the atom of each row)
    """
    frac = atoms.get_scaled_positions(wrap=True)
    order = list(range(len(atoms)))
    if type_names is not None:
        order.sort(key=lambda i: type_names.index(type_of_atom[i]))
    atmicx = []
    for i in order:
        x = frac[i]
        atmicx.append(["%10.8fa" % x[0], "%10.8fb" % x[1], "%10.8fc" % x[2], type_of_atom[i]])
    return {"natm": len(atoms), "atmicx": atmicx, "atom_order": order}


def atoms_to_kkr_param(atoms: Atoms, Vc: str = "Og", type_mode: str = "symmetry",
                       symprec: float = 1e-5, type_params: dict | None = None,
                       type_order: str = "electronegativity"):
    """make the structure part of the AkaiKKR input dict from ase.Atoms.

    Returns:
        dict, list: parameters (the same keys as Cif2Kkr.ak_cif2kkrparam), type name of each atom
    """
    param = atoms_to_kkr_lattice(atoms)
    types, type_of_atom = atoms_to_kkr_types(atoms, Vc=Vc, type_mode=type_mode,
                                             symprec=symprec, type_params=type_params,
                                             type_order=type_order)
    param.update(types)
    param.update(atoms_to_kkr_atmicx(atoms, type_of_atom, types["type"]))
    return param, type_of_atom


def check_kkr_output_structure(job, outfile, param: dict, tol: float = 1e-4):
    """check that AkaiKKR interpreted the structure as given.

    The primitive translation vectors and the atomic positions written in the
    output are compared with r1..r3 and atmicx of the input dict.

    Args:
        job (AkaikkrJob): job of the run directory
        outfile (str): output filename
        param (dict): input dict with r1, r2, r3 and atmicx
        tol (float, optional): tolerance in units of a. Defaults to 1e-4.

    Raises:
        KKRStructureMismatchError: mismatch
    """
    r_in = np.array([param["r1"], param["r2"], param["r3"]], dtype=float)
    r_out = np.array(job.get_prim_vec(outfile, unitof="relative"), dtype=float)
    if r_out.shape != (3, 3) or np.max(np.abs(r_out - r_in)) > tol:
        raise KKRStructureMismatchError(
            "primitive vectors differ.\ninput={}\noutput={}".format(r_in.tolist(), r_out.tolist()))

    coords, names = job.get_atom_coord(outfile)
    if len(coords) != len(param["atmicx"]):
        raise KKRStructureMismatchError(
            "number of atoms differ. input={} output={}".format(len(param["atmicx"]), len(coords)))
    inv_r = np.linalg.inv(r_in)
    for i, (c_out, name_out, atm) in enumerate(zip(coords, names, param["atmicx"])):
        frac = np.array([float(s[:-1]) for s in atm[:3]])
        c_in = np.dot(frac, r_in)
        d = np.dot(np.array(c_out) - c_in, inv_r)
        if np.max(np.abs(d - np.round(d))) > tol:
            raise KKRStructureMismatchError(
                "position of atom {} differs. input={} output={}".format(i, c_in.tolist(), c_out))
        if name_out != atm[3]:
            raise KKRStructureMismatchError(
                "type of atom {} differs. input={} output={}".format(i, atm[3], name_out))
