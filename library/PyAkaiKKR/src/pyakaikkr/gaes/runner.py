# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""go / dos / j in one directory with reuse of finished, identical runs."""
import io
import os
import shutil
from copy import deepcopy

from ..AkaiKkr import AkaikkrJob
from ..Error import KKRFailedExecutionError, KKRValueAquisitionError, GaesError
from .gap import gap_regions

DEFAULT_FILES = {"inputcard_go": "inputcard_go", "out_go": "out_go.log",
                 "inputcard_dos": "inputcard_dos", "out_dos": "out_dos.log",
                 "inputcard_j": "inputcard_j", "out_j": "out_j.log",
                 "potential": "pot.dat"}


def output_finished(job, outfile):
    """True if the output exists, ends normally (sbtime report) and has no ***err / option error."""
    path = os.path.join(job.path_dir, outfile)
    if not os.path.isfile(path):
        return False
    try:
        data = job._read(outfile)
    except OSError:
        return False
    if not any("sbtime report" in line for line in data):
        return False
    if any(line.startswith(" ***err in") for line in data):
        return False
    if job.check_option_error(outfile) is not None:
        return False
    return True


class KkrRunner:
    """run go -> dos -> j of one parameter set in one directory.

    A step is skipped when its inputcard on disk is identical to the one to be written
    and its output finished normally; otherwise the old output is kept as
    <out>.bak-<n> and specx runs again.
    """

    def __init__(self, akaikkr_exe, directory, param_go, ewidth_dos=3.0, edelt_dos=None, files=None,
                 compat=False, logger=None):
        self.akaikkr_exe = akaikkr_exe
        self.directory = directory
        self.param_go = deepcopy(param_go)
        self.param_go["go"] = "go"
        self.ewidth_dos = ewidth_dos
        self.edelt_dos = edelt_dos     # None: the go edelt (2019 behaviour)
        self.files = dict(DEFAULT_FILES)
        if files:
            self.files.update(files)
        self.compat = compat
        self.logger = logger
        self.job = AkaikkrJob(directory)
        self.ran = {}          # step -> bool (specx executed this time)
        self.ewidth_dos_used = None
        self.ref_effective = None   # -emin / ewidth_dos of the dos actually run

    # ---- helpers ----
    def _log(self, msg):
        if self.logger is not None:
            self.logger.info("%s: %s", self.directory, msg)

    def _param(self, go, ewidth=None):
        dic = deepcopy(self.param_go)
        dic["go"] = go
        if ewidth is not None:
            dic["ewidth"] = ewidth
        if go == "dos" and self.edelt_dos is not None:
            dic["edelt"] = self.edelt_dos
        if go != "go":
            dic["record"] = "2nd"
        return dic

    def _card_text(self, dic):
        buf = io.StringIO()
        self.job.make_inputcard(dic, buf)
        return buf.getvalue()

    def _same_inputcard(self, name, text):
        path = os.path.join(self.directory, name)
        if not os.path.isfile(path):
            return False
        with open(path) as f:
            return f.read() == text

    def _backup(self, name):
        path = os.path.join(self.directory, name)
        if not os.path.isfile(path):
            return
        n = 1
        while os.path.exists("{}.bak-{}".format(path, n)):
            n += 1
        shutil.move(path, "{}.bak-{}".format(path, n))

    def _run_step(self, step, dic, force=False, copy_potential_from=None):
        """write inputcard and run specx unless the identical run already finished."""
        os.makedirs(self.directory, exist_ok=True)
        inputcard = self.files["inputcard_" + step]
        outfile = self.files["out_" + step]
        text = self._card_text(dic)
        if not force and self._same_inputcard(inputcard, text) and output_finished(self.job, outfile):
            self._log("{}: reuse finished run".format(step))
            self.ran[step] = False
            return
        self._backup(outfile)
        if copy_potential_from is not None:
            src = os.path.join(copy_potential_from, self.files["potential"])
            if not os.path.isfile(src):
                raise GaesError("potential to copy not found: {}".format(src))
            shutil.copyfile(src, os.path.join(self.directory, self.files["potential"]))
        with open(os.path.join(self.directory, inputcard), "w") as f:
            f.write(text)
        self._log("{}: run specx".format(step))
        self.job.data = None
        self.job.outfile = None
        self.job.run(self.akaikkr_exe, inputcard, outfile)
        self.job.data = None
        self.job.outfile = None
        self.ran[step] = True

    # ---- steps ----
    def run_go(self, copy_potential_from=None, force=False):
        """go; returns True if converged."""
        self._run_step("go", self._param("go"), force=force, copy_potential_from=copy_potential_from)
        return self.converged()

    def converged(self):
        return self.job.get_convergence(self.files["out_go"])

    def run_dos(self, force=False, ewidth_dos=None, min_ewidth_dos=None, step=0.25):
        """dos with ewidth_dos (default self.ewidth_dos).

        If specx stops (e.g. '***err in reconf' of AkaiKKR >= 2020 when the window
        crosses a core level) the window is narrowed by `step` down to min_ewidth_dos
        and retried. The value actually used is in self.ewidth_dos_used.
        """
        ew = self.ewidth_dos if ewidth_dos is None else ewidth_dos
        while True:
            try:
                self._run_step("dos", self._param("dos", ewidth=ew), force=force)
                self.ewidth_dos_used = ew
                try:
                    e, _ = self.dos_curve()
                    # the mesh starts half a step inside the window: emin = -(ref*ew) + step/2
                    step = float(e[1] - e[0]) if len(e) > 1 else 0.0
                    self.ref_effective = (-float(e.min()) + 0.5 * step) / ew
                except Exception:  # noqa: BLE001
                    self.ref_effective = None
                return ew
            except KKRFailedExecutionError as e:
                self._log("dos with ewidth_dos={} failed: {}".format(ew, e))
                if min_ewidth_dos is None or ew - step < min_ewidth_dos - 1e-9:
                    raise
                ew = round(ew - step, 6)
                force = True

    def run_j(self, force=False):
        self._run_step("j", self._param("j"), force=force)

    def run_all(self, copy_potential_from=None, with_j=True, dos_always=False, min_ewidth_dos=None):
        """go, then dos and j.

        dos and j run when go converged; with compat=True or dos_always=True the dos
        runs regardless (the scheme needs a DOS to judge an unconverged run too).
        Returns True if go converged.
        """
        converged = self.run_go(copy_potential_from=copy_potential_from)
        if converged or self.compat or dos_always:
            self.run_dos(min_ewidth_dos=min_ewidth_dos)
        if with_j and (converged or self.compat):
            self.run_j()
        return converged

    # ---- results ----
    def has(self, step):
        return output_finished(self.job, self.files["out_" + step])

    def dos_curve(self, spin_sum=True):
        """(energy, total DOS) of this directory's dos output."""
        from .gap import dos_curves_from_outputs
        e, curves, _ = dos_curves_from_outputs([(self.job, self.files["out_dos"])], spin_sum=spin_sum)
        return e, curves[0]

    def result(self, dosth=1e-3, save_csv=True):
        """one result row (see docs/ewidth_tuning_scheme.md section 9)."""
        job = self.job
        og = self.files["out_go"]
        row = {"directory": self.directory}
        if not os.path.isfile(os.path.join(self.directory, og)):
            row["go_output"] = False
            return row
        row["go_output"] = True
        try:
            row["converged"] = job.get_convergence(og)
        except KKRValueAquisitionError:
            row["converged"] = None
        for name, fn in (("ewidth", job.get_ewidth), ("edelt", job.get_edelt),
                         ("a_bohr", job.get_lattice_constant), ("volume_bohr3", job.get_unitcell_volume),
                         ("total_energy_Ry", job.get_total_energy), ("total_moment", job.get_total_moment),
                         ("emesh", job.get_emesh_param)):
            try:
                row[name] = fn(og)
            except (KKRValueAquisitionError, IndexError, ValueError):
                row[name] = None
        try:
            err = job.get_rms_error(og)
            row["n_iter"] = len(err)
            row["last_err"] = err[-1]
        except KKRValueAquisitionError:
            row["n_iter"] = row["last_err"] = None
        try:
            row["component_moment"] = job.get_component_moment(og)
        except KKRValueAquisitionError:
            row["component_moment"] = None
        row["option"] = job.get_option(og)
        if self.has("j"):
            oj = self.files["out_j"]
            try:
                row["Tc_K"] = job.get_curie_temperature(oj)
            except KKRValueAquisitionError:
                row["Tc_K"] = None
            if save_csv:
                try:
                    df = job.get_jij_as_dataframe(oj)
                    p = os.path.join(self.directory, "jij.csv")
                    df.to_csv(p, index=False)
                    row["jij_csv"] = p
                except (KKRValueAquisitionError, StopIteration, ValueError):
                    row["jij_csv"] = None
        if self.has("dos"):
            od = self.files["out_dos"]
            try:
                e, d = self.dos_curve()
                regions = gap_regions(e, [d], dosth=dosth)
                row["ewidth_dos"] = job.get_ewidth(od)
                row["dos_window"] = (float(e.min()), float(e.max()))
                row["low_dos_regions"] = [g.as_tuple() for g in regions]
                if save_csv:
                    p = os.path.join(self.directory, "dos.csv")
                    job.get_dos(od, "dataframe").to_csv(p, index=False)
                    row["dos_csv"] = p
            except (KKRValueAquisitionError, GaesError):
                row["low_dos_regions"] = None
        return row
