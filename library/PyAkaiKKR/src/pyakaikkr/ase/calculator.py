# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""ASE calculator for AkaiKKR (energy only, CPA supported)."""
from __future__ import annotations

import os
from copy import deepcopy

import numpy as np
from ase import units
from ase.calculators.calculator import (FileIOCalculator, InputError, CalculationFailed,
                                        SCFError, ReadError, all_changes)

from ..AkaiKkr import AkaikkrJob
from .occupancy import same_occupancy
from .structure import atoms_to_kkr_param, check_kkr_output_structure

_INPUT_KEYS = ["go", "potentialfile", "edelt", "ewidth", "reltyp", "sdftyp", "magtyp",
               "outtyp", "bzqlty", "maxitr", "pmix"]
_TYPE_KEYS = ["ntyp", "type", "ncmp", "anclr", "conc", "magtyp"]


class AkaiKKR(FileIOCalculator):
    """ASE calculator running specx of AkaiKKR.

    Properties: energy (total energy in eV per cell), magmom (total moment per cell),
    magmoms (per atom; concentration-weighted average of the component moments
    of its type). Forces and stress are not available.

    The structure is passed with brvtyp=aux, a=|cell[0]| (bohr) and r1..r3=cell/|cell[0]|.
    rmt (if given through type_params) is in units of this a. Types are ordered by
    electronegativity (type_order); with rmt=0 the automatic muffin-tin radii depend on it.
    Partial occupancies follow the ASE convention (atoms.info['occupancy'] and
    atoms.arrays['spacegroup_kinds']); see pyakaikkr.ase.set_occupancy.

    The command is given by ``command=``, the environment variable ASE_AKAIKKR_COMMAND
    or the [akaikkr] section of the ASE config. PREFIX is replaced by the label.
    Example: command="/path/to/specx < PREFIX.in > PREFIX.out".
    """
    implemented_properties = ["energy", "magmom", "magmoms"]
    _legacy_default_command = "specx < PREFIX.in > PREFIX.out"
    default_parameters = dict(
        go="go", potentialfile="pot.dat",
        edelt=1e-3, ewidth=1.0, reltyp="sra", sdftyp="mjw", magtyp="mag",
        record=None, outtyp="update", bzqlty=6, maxitr=200, pmix=0.02,
        option=None,
        type_mode="symmetry", symprec=1e-5, type_params=None, Vc="Og",
        type_order="electronegativity",
        reuse_potential=True, allow_unconverged=False,
    )

    def _get_name(self) -> str:
        return "akaikkr"

    def __init__(self, restart=None, label="akaikkr", atoms=None, command=None,
                 profile=None, directory=".", **kwargs):
        super().__init__(restart=restart, label=label, atoms=atoms, command=command,
                         profile=profile, directory=directory, **kwargs)
        if self.prefix is None:
            self.prefix = "akaikkr"
        self._last_dic = None
        self._type_of_atom = None
        self.job = None

    # ----- cache -----
    def check_state(self, atoms, tol=1e-15):
        changes = super().check_state(atoms, tol)
        if self.atoms is not None and not same_occupancy(self.atoms, atoms):
            return all_changes[:]
        return changes

    # ----- files -----
    @property
    def infile(self):
        return self.prefix + ".in"

    @property
    def outfile(self):
        return self.prefix + ".out"

    @staticmethod
    def _same_types(dic1, dic2):
        if dic1 is None or dic2 is None:
            return False
        return all(dic1.get(k) == dic2.get(k) for k in _TYPE_KEYS)

    def _decide_record(self, dic):
        p = self.parameters
        if p["record"] is not None:
            return p["record"]
        potfile = os.path.join(self.directory, dic["potentialfile"])
        if p["reuse_potential"] and os.path.isfile(potfile) and self._same_types(self._last_dic, dic):
            return "2nd"
        return "init"

    def make_input_dict(self, atoms):
        """make the AkaiKKR input dict for atoms (without deciding record)."""
        p = self.parameters
        if str(p["go"])[:3] == "spc":
            raise InputError("go=spc* (needs a k-path) is not supported by the calculator.")
        param, type_of_atom = atoms_to_kkr_param(
            atoms, Vc=p["Vc"], type_mode=p["type_mode"], symprec=p["symprec"],
            type_params=p["type_params"], type_order=p["type_order"])
        job = AkaikkrJob(self.directory)
        dic = deepcopy(job.default)
        dic.update(param)
        for key in _INPUT_KEYS:
            dic[key] = p[key]
        if p["option"]:
            dic["option"] = p["option"]
        if p["go"] == "fsm":
            if "fspin" not in p:
                raise InputError("go=fsm needs the parameter fspin.")
            dic["fspin"] = p["fspin"]
        return dic, type_of_atom

    def write_input(self, atoms, properties=None, system_changes=None):
        super().write_input(atoms, properties, system_changes)
        dic, type_of_atom = self.make_input_dict(atoms)
        dic["record"] = self._decide_record(dic)
        job = AkaikkrJob(self.directory)
        with open(os.path.join(self.directory, self.infile), "w") as f:
            job.make_inputcard(dic, f)
        self._last_dic = deepcopy(dic)
        self._type_of_atom = list(type_of_atom)

    def read_results(self):
        p = self.parameters
        outpath = os.path.join(self.directory, self.outfile)
        if not os.path.isfile(outpath):
            raise ReadError(f"{outpath} not found.")
        job = AkaikkrJob(self.directory)
        if job.check_stopped_by_errtrp(self.outfile):
            raise CalculationFailed("specx stopped by errtrp: " + job._read(self.outfile)[-1])
        converged = job.get_convergence(self.outfile)
        if not converged and not p["allow_unconverged"]:
            raise SCFError("AkaiKKR did not converge (set allow_unconverged=True to accept).")
        if self._last_dic is not None:
            check_kkr_output_structure(job, self.outfile, self._last_dic)

        te = job.get_total_energy(self.outfile)
        self.results["energy"] = te * units.Rydberg
        self.results["magmom"] = job.get_total_moment(self.outfile)

        typeofsite = job.get_type_of_site(self.outfile)
        comps = job.get_component_moment(self.outfile)
        component_moments = {}
        weighted = {}
        k = 0
        for t in typeofsite:
            name = t["type"]
            component_moments[name] = {}
            wsum = 0.0
            for c in t["component"]:
                if k >= len(comps) or comps[k]["type"] != name:
                    raise ReadError("component moments do not match the type of site.")
                e = comps[k]
                component_moments[name][e["element"]] = {
                    "spin": e["spin"], "orbital": e["orbital"], "conc": c["conc"]}
                wsum += c["conc"] * e["spin"]
                k += 1
            weighted[name] = wsum
        if self._type_of_atom is not None:
            self.results["magmoms"] = np.array([weighted[n] for n in self._type_of_atom])
        self.results["component_moments"] = component_moments
        self.results["converged"] = converged
        self.results["type_of_site"] = typeofsite
        self.job = job
