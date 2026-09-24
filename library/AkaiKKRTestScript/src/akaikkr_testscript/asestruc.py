# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""structure parameters of AkaiKKR from ase.Atoms (the ASE backend of the test script)."""
import os
import warnings

from pyakaikkr import AkaikkrJob
from pyakaikkr.ase import atoms_to_kkr_param, to_primitive, check_kkr_output_structure

_GEOM_INPUTCARD = "inputcard_geom"
_GEOM_OUTPUTCARD = "out_geom.log"


def get_kkr_struc_from_ase(structure, akaikkr_exe: str, displc: bool,
                           use_bravais=False, remove_temperaryfiles=True,
                           Vc: str = "Og", directory: str = "temporary",
                           fmt=None, use_primitive=True,
                           type_mode="symmetry", symprec=1e-5,
                           check_geom=True) -> dict:
    """get kkr structure parameters from ase.Atoms or a structure file.

    The returned dict has the same keys as get_kkr_struc_from_cif().
    The lattice is always given as brvtyp=aux with r1, r2, r3 (use_bravais is ignored).
    If check_geom is True, specx is run with go=geom in directory and the structure
    in the output is compared with the input (as CompareCifKkr does for cif files).

    Args:
        structure (ase.Atoms or str): structure or filename readable by ase.io.read
        akaikkr_exe (str): specx path
        displc (bool): not used (displc is added by GoGo.execute). kept for compatibility.
        use_bravais (bool, optional): ignored. Defaults to False.
        remove_temperaryfiles (bool, optional): delete the geom files. Defaults to True.
        Vc (str, optional): element treated as vacancy (Z=0). Defaults to "Og".
        directory (str, optional): directory of the geom run. Defaults to "temporary".
        fmt (str, optional): format for ase.io.read. Defaults to None (guess).
        use_primitive (bool, optional): reduce to the primitive cell. Defaults to True.
        type_mode (str, optional): see pyakaikkr.ase.atoms_to_kkr_types. Defaults to "symmetry".
        symprec (float, optional): symmetry tolerance. Defaults to 1e-5.
        check_geom (bool, optional): run specx in geom mode to check. Defaults to True.

    Returns:
        dict: kkr structure parameters
    """
    import ase.io
    from ase import Atoms

    if use_bravais:
        warnings.warn("use_bravais is ignored by the ASE backend. brvtyp=aux is used.")
    if isinstance(structure, Atoms):
        atoms = structure
    else:
        atoms = ase.io.read(structure, format=fmt)
    if use_primitive:
        atoms = to_primitive(atoms, symprec=symprec)

    param, _ = atoms_to_kkr_param(atoms, Vc=Vc, type_mode=type_mode, symprec=symprec)

    if check_geom:
        os.makedirs(directory, exist_ok=True)
        job = AkaikkrJob(directory)
        dic = dict(job.default)
        dic.update(param)
        dic["go"] = "geom"
        job.make_inputcard(dic, _GEOM_INPUTCARD)
        job.run(akaikkr_exe, _GEOM_INPUTCARD, _GEOM_OUTPUTCARD)
        check_kkr_output_structure(job, _GEOM_OUTPUTCARD, dic)
        if remove_temperaryfiles:
            for name in [_GEOM_INPUTCARD, _GEOM_OUTPUTCARD]:
                path = os.path.join(directory, name)
                if os.path.isfile(path):
                    os.remove(path)
            try:
                os.rmdir(directory)
            except OSError:
                pass
    return param
