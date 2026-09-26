# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""partial occupancy (CPA) support for ase.Atoms.

ASE convention for partial occupancies (ase.io.cif, ase.spacegroup.crystal):

- ``atoms.info['occupancy']``: ``{"<kind>": {"<symbol>": occupancy, ...}, ...}``.
  The keys are kind numbers as str.
- ``atoms.arrays['spacegroup_kinds']``: kind number of each atom.
  If it is absent, kind = atom index.
- ``atoms.symbols[i]`` is one of the elements of its kind.

They are mapped to the AkaiKKR type/ncmp/anclr/conc parameters here.
"""
from __future__ import annotations

import numpy as np
import spglib
from ase import Atoms
from ase.calculators.calculator import InputError

from ..ElementKkr import ElementKKR

_OCC_TOL = 1e-6
_VC_LABEL = "Vc"


def _dataset_attr(dataset, name):
    """spglib>=2 returns a SpglibDataset object, older versions a dict."""
    if dataset is None:
        raise InputError("spglib failed to find the symmetry of the structure.")
    if isinstance(dataset, dict):
        return dataset[name]
    return getattr(dataset, name)


def get_kind_occupancies(atoms: Atoms):
    """get the occupancy of each atom and classify atoms into kinds.

    Atoms with the same occupancy dict belong to the same kind.

    Args:
        atoms (Atoms): structure

    Raises:
        InputError: inconsistent occupancy information

    Returns:
        np.ndarray, list: kind id (0..K-1, in the order of first appearance) of each atom,
            occupancy dict {symbol: occupancy} of each kind
    """
    n = len(atoms)
    symbols = list(atoms.get_chemical_symbols())
    occ_info = atoms.info.get("occupancy", None)
    kinds_arr = atoms.arrays.get("spacegroup_kinds", None)

    occ_per_atom = []
    for i in range(n):
        if occ_info is None:
            occ = {symbols[i]: 1.0}
        else:
            if kinds_arr is not None:
                key = str(int(kinds_arr[i]))
            else:
                key = str(i)
            if key not in occ_info:
                raise InputError(
                    f"atoms.info['occupancy'] has no entry for kind '{key}' (atom {i}).")
            occ = {str(s): float(o) for s, o in occ_info[key].items()}
            if symbols[i] not in occ:
                raise InputError(
                    f"symbol {symbols[i]} of atom {i} is not in its occupancy {occ}.")
        if len(occ) == 0:
            raise InputError(f"empty occupancy for atom {i}.")
        total = 0.0
        for s, o in occ.items():
            if not (0.0 < o <= 1.0 + _OCC_TOL):
                raise InputError(f"occupancy of {s} at atom {i} must be in (0,1], got {o}.")
            total += o
        if total > 1.0 + _OCC_TOL:
            raise InputError(f"sum of occupancies at atom {i} is {total} > 1.")
        occ_per_atom.append(occ)

    keys = [tuple(sorted((s, round(o, 8)) for s, o in occ.items())) for occ in occ_per_atom]
    uniq = []
    kind_of_atom = []
    for k in keys:
        if k not in uniq:
            uniq.append(k)
        kind_of_atom.append(uniq.index(k))
    kind_occ = [dict(k) for k in uniq]
    return np.array(kind_of_atom, dtype=int), kind_occ


def _attach_kinds(atoms: Atoms, kind_of_atom, kind_occ) -> Atoms:
    """write kind information back to atoms in the ASE convention."""
    atoms.info["occupancy"] = {str(k): dict(occ) for k, occ in enumerate(kind_occ)}
    if atoms.has("spacegroup_kinds"):
        atoms.set_array("spacegroup_kinds", np.asarray(kind_of_atom, dtype=int))
    else:
        atoms.new_array("spacegroup_kinds", np.asarray(kind_of_atom, dtype=int))
    return atoms


def set_occupancy(atoms: Atoms, occupancy: dict) -> Atoms:
    """attach partial occupancies to a copy of atoms in the ASE convention.

    Example: set_occupancy(bulk("Ni", "fcc", a=3.571), {0: {"Fe": 0.1, "Ni": 0.9}})

    Args:
        atoms (Atoms): structure
        occupancy (dict): {atom index: {symbol: occupancy}}. Atoms not listed keep occupancy 1.

    Returns:
        Atoms: a copy with atoms.info['occupancy'] and atoms.arrays['spacegroup_kinds'].
            The symbol of each listed atom is set to the element with the largest occupancy
            if the current symbol is not in its occupancy dict.
    """
    atoms = atoms.copy()
    n = len(atoms)
    symbols = list(atoms.get_chemical_symbols())
    occ_per_atom = [{symbols[i]: 1.0} for i in range(n)]
    for i, occ in occupancy.items():
        i = int(i)
        if not (0 <= i < n):
            raise InputError(f"atom index {i} is out of range.")
        occ = {str(s): float(o) for s, o in occ.items()}
        occ_per_atom[i] = occ
        if symbols[i] not in occ:
            symbols[i] = max(occ.items(), key=lambda so: so[1])[0]
    atoms.set_chemical_symbols(symbols)
    atoms.info.pop("occupancy", None)
    if atoms.has("spacegroup_kinds"):
        del atoms.arrays["spacegroup_kinds"]
    # classify with a temporary info
    tmp = atoms.copy()
    tmp.info["occupancy"] = {str(i): occ for i, occ in enumerate(occ_per_atom)}
    kind_of_atom, kind_occ = get_kind_occupancies(tmp)
    return _attach_kinds(atoms, kind_of_atom, kind_occ)


def to_primitive(atoms: Atoms, symprec: float = 1e-5) -> Atoms:
    """reduce atoms to a primitive cell keeping the partial occupancies.

    spglib.find_primitive is called with the kind ids as atomic numbers so that
    sites with different occupancies are never merged. The orientation of the
    input cell is kept. The atoms of the primitive cell are ordered by the
    index of the corresponding atom in the input cell.

    Args:
        atoms (Atoms): structure
        symprec (float, optional): symmetry tolerance in Angstrom. Defaults to 1e-5.

    Returns:
        Atoms: primitive cell with occupancy information
    """
    kind_of_atom, kind_occ = get_kind_occupancies(atoms)
    cell = (np.array(atoms.cell), atoms.get_scaled_positions(wrap=True),
            (kind_of_atom + 1).tolist())
    res = spglib.find_primitive(cell, symprec=symprec)
    if res is None:
        raise InputError("spglib.find_primitive failed.")
    lattice, positions, numbers = res
    kinds_prim = np.array(numbers, dtype=int) - 1

    # map each primitive atom to the smallest index of the equivalent input atom
    inv_lat = np.linalg.inv(lattice)
    cart_in = atoms.get_positions(wrap=True)
    cart_prim = np.dot(positions, lattice)
    origin = []
    for j in range(len(kinds_prim)):
        found = None
        for i in range(len(atoms)):
            if kind_of_atom[i] != kinds_prim[j]:
                continue
            d = np.dot(cart_in[i] - cart_prim[j], inv_lat)
            if np.all(np.abs(d - np.round(d)) < 1e-4):
                found = i
                break
        if found is None:
            raise InputError("failed to map a primitive atom to the input atoms.")
        origin.append(found)
    order = np.argsort(origin, kind="stable")

    symbols = []
    for j in order:
        occ = kind_occ[kinds_prim[j]]
        cand = [s for s in occ if s == atoms[origin[j]].symbol]
        symbols.append(cand[0] if cand else max(occ.items(), key=lambda so: so[1])[0])
    prim = Atoms(symbols=symbols, scaled_positions=positions[order], cell=lattice, pbc=True)
    prim.info.update({k: v for k, v in atoms.info.items() if k != "occupancy"})
    return _attach_kinds(prim, kinds_prim[order], kind_occ)


def _format_occ(occ: float) -> str:
    return "%.6g" % occ


def _average_electronegativity(occ: dict) -> float:
    """occupancy-weighted Pauling electronegativity (as pymatgen Composition.average_electroneg).
    NaN (e.g. Og used as vacancy) is mapped to +inf so that such types come last."""
    from pymatgen.core.periodic_table import Element
    num = 0.0
    den = 0.0
    for s, o in occ.items():
        x = Element(s).X
        if x is None or x != x:
            return float("inf")
        num += x * o
        den += o
    return num / den if den > 0 else float("inf")


def atoms_to_kkr_types(atoms: Atoms, Vc: str = "Og", type_mode: str = "symmetry",
                       symprec: float = 1e-5, type_params: dict | None = None,
                       maxlen: int = 40, type_order: str = "electronegativity"):
    """make the type section of AkaiKKR (ntyp, type, ncmp, rmt, field, mxl, anclr, conc).

    A type is a set of atoms which share one potential. Atoms belong to the same type
    when they have the same occupancy (kind) and are symmetrically equivalent
    (type_mode="symmetry"), or when they have the same kind (type_mode="kind"),
    or never (type_mode="atom").

    Components in a type are ordered by the atomic number. If the sum of occupancies
    is less than 1, a vacancy (element Vc, Z=0) is added as the last component.
    rmt is the average, field and mxl are the maxima over the elements (as Cif2Kkr).

    Args:
        atoms (Atoms): structure with optional occupancy information
        Vc (str, optional): element treated as Z=0 (vacancy). Defaults to "Og".
        type_mode (str, optional): "symmetry", "kind" or "atom". Defaults to "symmetry".
        symprec (float, optional): symmetry tolerance in Angstrom. Defaults to 1e-5.
        type_params (dict, optional): {type name or element symbol: {"rmt":, "field":, "mxl":}}
            to override the default values. Element keys are applied first.
        maxlen (int, optional): maximum length of a type name. Defaults to 40.
        type_order (str, optional): order of the types. "electronegativity": ascending
            occupancy-weighted Pauling electronegativity, ties in the order of appearance
            (the same order as pymatgen sorts sites, so the same as Cif2Kkr).
            "appearance": order of the first atom of each type. Defaults to "electronegativity".
            The order matters when rmt=0 (muffin-tin radii are determined type by type).

    Returns:
        dict, list: type parameters, type name of each atom
    """
    n = len(atoms)
    if n == 0:
        raise InputError("no atoms.")
    kind_of_atom, kind_occ = get_kind_occupancies(atoms)
    elementkkr = ElementKKR(Vc=Vc)
    for occ in kind_occ:
        for s in occ:
            if s not in elementkkr.dict:
                raise InputError(f"unknown element {s}.")

    if type_mode == "symmetry":
        cell = (np.array(atoms.cell), atoms.get_scaled_positions(wrap=True),
                (kind_of_atom + 1).tolist())
        dataset = spglib.get_symmetry_dataset(cell, symprec=symprec)
        labels = np.array(_dataset_attr(dataset, "equivalent_atoms"), dtype=int)
    elif type_mode == "kind":
        labels = kind_of_atom
    elif type_mode == "atom":
        labels = np.arange(n)
    else:
        raise ValueError(f"unknown type_mode {type_mode}")

    uniq = []
    type_of_atom = []
    for lab in labels:
        if lab not in uniq:
            uniq.append(lab)
        type_of_atom.append(uniq.index(lab))
    ntyp = len(uniq)
    rep_atoms = [type_of_atom.index(t) for t in range(ntyp)]
    for i in range(n):
        if kind_of_atom[i] != kind_of_atom[rep_atoms[type_of_atom[i]]]:
            raise InputError("atoms of different kinds fell into one type.")

    if type_order == "electronegativity":
        xs = [_average_electronegativity(kind_occ[kind_of_atom[rep_atoms[it]]]) for it in range(ntyp)]
        order = sorted(range(ntyp), key=lambda it: xs[it])   # stable: ties keep appearance order
        newindex = {old: new for new, old in enumerate(order)}
        type_of_atom = [newindex[t] for t in type_of_atom]
        rep_atoms = [rep_atoms[old] for old in order]
    elif type_order != "appearance":
        raise ValueError(f"unknown type_order {type_order}")

    type_params = type_params or {}
    param = {"ntyp": ntyp, "type": [], "ncmp": [], "rmt": [], "field": [], "mxl": [],
             "anclr": [], "conc": []}
    names = []
    for it in range(ntyp):
        occ = kind_occ[kind_of_atom[rep_atoms[it]]]
        comps = sorted(occ.items(), key=lambda so: (elementkkr.getAtomicNumber(so[0]), so[0]))
        total = sum(o for _, o in comps)
        anclr = [elementkkr.getAtomicNumber(s) for s, _ in comps]
        conc = [round(o * 100.0, 8) for _, o in comps]
        formula = ""
        for s, o in comps:
            formula += s
            if not (len(comps) == 1 and abs(o - 1.0) < _OCC_TOL):
                formula += _format_occ(o)
        if total < 1.0 - _OCC_TOL:
            anclr.append(0)
            conc.append(round((1.0 - total) * 100.0, 8))
            formula += _VC_LABEL + _format_occ(1.0 - total)
        diff = 100.0 - sum(conc)
        if abs(diff) < 1e-6:
            conc[-1] = round(conc[-1] + diff, 8)
        name = "{}_{}".format(formula, it)
        if len(name) > maxlen:
            name = "T{}".format(it)
        names.append(name)

        elements = [s for s, _ in comps]
        rmt = float(np.mean([elementkkr.getAtomicRMT(s) for s in elements]))
        field = float(max(abs(elementkkr.getAtomicField(s)) for s in elements))
        mxl = int(max(elementkkr.getAtomicLMax(s) for s in elements))
        tp = {"rmt": rmt, "field": field, "mxl": mxl}
        for s in elements:
            if s in type_params:
                tp.update(type_params[s])
        if name in type_params:
            tp.update(type_params[name])

        param["type"].append(name)
        param["ncmp"].append(len(anclr))
        param["rmt"].append(tp["rmt"])
        param["field"].append(tp["field"])
        param["mxl"].append(tp["mxl"])
        param["anclr"].append(anclr)
        param["conc"].append(conc)
    return param, [names[t] for t in type_of_atom]


def same_occupancy(atoms1: Atoms, atoms2: Atoms) -> bool:
    """compare the occupancy information of two Atoms (used by the calculator cache)."""
    if atoms1 is None or atoms2 is None or len(atoms1) != len(atoms2):
        return False
    try:
        k1, o1 = get_kind_occupancies(atoms1)
        k2, o2 = get_kind_occupancies(atoms2)
    except InputError:
        return False
    if len(o1) != len(o2) or not np.array_equal(k1, k2):
        return False
    for a, b in zip(o1, o2):
        if set(a) != set(b):
            return False
        for s in a:
            if abs(a[s] - b[s]) > 1e-8:
                return False
    return True
