# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.


from .BasePlotter import BasePlotter
import os
import warnings
import matplotlib.pyplot as plt
import numpy as np

from .BasePlotter import BaseEXPlotter
from .AkaiKkr import AkaikkrJob
from .Error import KKRValueAquisitionError

from .plot import mark_ewidth_go, plot_dos as _draw_dos, plot_pdos as _draw_pdos, EWIDTH_LINE_LABEL

_EWIDTH_LINE_LABEL = EWIDTH_LINE_LABEL   # label of the ewidth line (kept for callers and tests)

_DEFAULT_GO_OUTFILE = "out_go.log"


def _mark_ewidth_go(ax, energy, ewidth_go):
    """dash-dot vertical line at E - EF = -|ewidth_go| (see pyakaikkr.plot.mark_ewidth_go).

    AkaiKKR integrates the charge over [EF - ewidth, EF] in the go run, so the line shows the
    bottom of the SCF energy contour; the ewidth of the go run, not of the dos run, must be given.
    Returns True if the line is drawn.
    """
    return mark_ewidth_go(ax, energy, ewidth_go)


def resolve_ewidth_go(directory, ewidth_go=None, go_outfile=_DEFAULT_GO_OUTFILE,
                      read_go_outfile=True):
    """decide the ewidth of the go run used for the vertical line of DOS plots.

    Priority: ewidth_go argument, then the ewidth read from go_outfile in directory.
    A warning is issued and None is returned if neither is available. The ewidth
    of the dos run is never used because -ewidth_dos always lies outside the DOS mesh.

    Args:
        directory (str): directory containing go_outfile.
        ewidth_go (float, optional): explicit ewidth of the go run. Defaults to None.
        go_outfile (str, optional): output filename of the go run. Defaults to "out_go.log".
        read_go_outfile (bool, optional): read go_outfile when ewidth_go is None. Defaults to True.

    Returns:
        float or None: ewidth of the go run, or None.
    """
    if ewidth_go is not None:
        return float(ewidth_go)
    if not read_go_outfile or go_outfile is None:
        return None
    filepath = os.path.join(directory, go_outfile)
    if not os.path.isfile(filepath):
        warnings.warn(
            "{} not found. No ewidth line is drawn on the DOS plot. "
            "Give ewidth_go or go_outfile.".format(filepath))
        return None
    try:
        return AkaikkrJob(directory).get_ewidth(go_outfile)
    except KKRValueAquisitionError:
        warnings.warn(
            "ewidth is not found in {}. No ewidth line is drawn on the DOS plot.".format(filepath))
        return None


def _plot_dos(energy, dos_block, output_direcotry=None,
              yscale="log", figsize=(5, 3), ewidth_go=None):
    """total DOS figure (dos.png): one panel per spin, drawn by pyakaikkr.plot.plot_dos."""
    nspin = len(dos_block)
    fig, axes = plt.subplots(1, nspin, figsize=(figsize[0] * nspin, figsize[1]), squeeze=False)
    for ispin, ax in enumerate(axes[0]):
        _draw_dos(ax, energy, [dos_block[ispin]], yscale=yscale, efermi=False)
        if nspin > 1:
            ax.set_title("up" if ispin == 0 else "dn")
        ax.set_xlabel("E(Ry)-EF")
        ax.set_ylabel("DOS")
        if _mark_ewidth_go(ax, energy, ewidth_go):
            ax.legend(frameon=False, fontsize=8)
    if output_direcotry is None:
        output_direcotry = "."
    fig.tight_layout()
    imgfilepath = os.path.join(output_direcotry, "dos.png")
    fig.savefig(imgfilepath)
    print("  saved to", imgfilepath)
    fig.clf()
    plt.close(fig)
    return imgfilepath


def _EX_plot_dos(directory, outfile, output_direcotry=None,
                 yscale="log", figsize=(5, 3), ewidth_go=None):
    job = AkaikkrJob(directory)
    energy, dos_block = job.get_dos_as_list(outfile)
    if output_direcotry is None:
        output_direcotry = directory
    return _plot_dos(energy, dos_block, output_direcotry=output_direcotry,
                     yscale=yscale, figsize=figsize, ewidth_go=ewidth_go)


def _plot_pdos_all(energy, pdos_block, typeofsite, output_direcotry=None,
                   yscale="log", figsize=(5, 3), ewidth_go=None):
    """PDOS figures (pdos_<icmp>.png): one figure per component, one panel per spin, one curve per l,
    drawn by pyakaikkr.plot.plot_pdos."""
    serial_site = [shortname for site in typeofsite for shortname in site["comp_shortname"]]
    energy = np.array(energy)
    nspin = len(pdos_block)
    if output_direcotry is None:
        output_direcotry = "."
    imgfilepaths = []
    for icmp, title in enumerate(serial_site):
        fig, axes = plt.subplots(1, nspin, figsize=(figsize[0] * nspin, figsize[1]), squeeze=False)
        for ispin, ax in enumerate(axes[0]):
            _draw_pdos(ax, energy, np.array(pdos_block[ispin][icmp]), yscale=yscale, efermi=False)
            _mark_ewidth_go(ax, energy, ewidth_go)
            ax.legend()
            if nspin > 1:
                ax.set_title("up" if ispin == 0 else "dn")
            ax.set_xlabel("E(Ry)-EF")
            ax.set_ylabel("PDOS")
        if nspin > 1:
            fig.suptitle(title)
        else:
            axes[0][0].set_title(title)
        fig.tight_layout()
        imgfilepath = os.path.join(output_direcotry, "pdos_{}.png".format(icmp))
        fig.savefig(imgfilepath)
        print("  saved to", imgfilepath)
        imgfilepaths.append(imgfilepath)
        fig.clf()
        plt.close(fig)
    print()
    return imgfilepaths


def _EX_plot_pdos_all(directory, outfile, output_direcotry=None,
                      yscale="log", figsize=(5, 3), ewidth_go=None):
    job = AkaikkrJob(directory)
    typeofsite = job.get_type_of_site(outfile)
    energy, pdos_block = job.get_pdos_as_list(
        outfile, output_format="spin_separation")
    if output_direcotry is None:
        output_direcotry = directory
    return _plot_pdos_all(energy, pdos_block, typeofsite,
                          output_direcotry=output_direcotry,
                          yscale=yscale, figsize=figsize, ewidth_go=ewidth_go)


class PDosPlotter(BasePlotter):
    def __init__(self, output_directory):
        super().__init__(output_directory)

    def make(self, energy, pdos_block, typeofsites,
             yscale="log", figsize=(5, 3), ewidth_go=None):
        """make PDOS from output of akaikkr.get_pdos_as_list().

        energy, pdos_block is output of  akaikkr.get_pdos_as_list().
        typeofsites is output of akaikkr.get_type_of_site().

        Args:
            energy ([float]): energy.
            pdos_block ([[float]]): pdos_block.
            typeofsites (list): type of sites.
            yscale (str, optional): y scale. Defaults to "log".
            figsize (tuple, optional): figure size. Defaults to (5, 3).
            ewidth_go (float, optional): ewidth of the go run. A vertical line is
                drawn at E-EF=-|ewidth_go| if given. Defaults to None.
        """
        return _plot_pdos_all(energy, pdos_block, typeofsites,
                              output_direcotry=self.output_directory,
                              yscale=yscale, figsize=figsize, ewidth_go=ewidth_go)


class PDosEXPlotter(BaseEXPlotter):
    def __init__(self, directory, outfile, output_directory=None,
                 go_outfile=_DEFAULT_GO_OUTFILE, read_go_outfile=True):
        """
        Args:
            directory (str): directory containing outfile (and go_outfile).
            outfile (str): output filename of the dos run.
            output_directory (str, optional): output directory. Defaults to directory.
            go_outfile (str, optional): output filename of the go run in directory,
                from which ewidth of the go run is read. Defaults to "out_go.log".
            read_go_outfile (bool, optional): read go_outfile automatically. Defaults to True.
        """
        super().__init__(directory, outfile, output_directory)
        self.go_outfile = go_outfile
        self.read_go_outfile = read_go_outfile

    def make(self, output_directory=None, yscale="log", figsize=(5, 3), ewidth_go=None):
        """make PDOS from output file.

        A vertical line is drawn at E-EF=-|ewidth_go|. ewidth_go is taken from
        the argument, then from go_outfile. No line is drawn (with a warning) if
        neither is available.

        Args:
            output_directory (str, optional): output directory of the image file. Defaults to None.
            yscale (str, optional): y scale. Defaults to "log".
            figsize (tuple, optional): figure size. Defaults to (5, 3).
            ewidth_go (float, optional): ewidth of the go run. Defaults to None.
        """
        if output_directory is None:
            output_directory = self.output_directory
        ewidth_go = resolve_ewidth_go(self.directory, ewidth_go,
                                      self.go_outfile, self.read_go_outfile)
        return _EX_plot_pdos_all(self.directory, self.outfile,
                                 output_direcotry=output_directory, yscale=yscale,
                                 figsize=figsize, ewidth_go=ewidth_go)


class DosPlotter(BasePlotter):
    def __init__(self, output_directory):
        super().__init__(output_directory)

    def make(self, energy, dos_block, yscale="log",
             figsize=(5, 3), ewidth_go=None):
        """make DOS from akaikkr.get_dos_as_list()

        energy and dos_block is output of akaikkr.get_dos_as_list()
        Args:
            energy ([float]): energy mesh.
            dos_block ([[float]]): dos block.
            yscale (str, optional): y scale. Defaults to "log".
            figsize (tuple, optional): figure size. Defaults to (5, 3).
            ewidth_go (float, optional): ewidth of the go run. A vertical line is
                drawn at E-EF=-|ewidth_go| if given. Defaults to None.
        """
        return _plot_dos(energy, dos_block,
                         output_direcotry=self.output_directory, yscale=yscale,
                         figsize=figsize, ewidth_go=ewidth_go)


class DosEXPlotter(BaseEXPlotter):
    def __init__(self, directory, outfile, output_directory=None,
                 go_outfile=_DEFAULT_GO_OUTFILE, read_go_outfile=True):
        """
        Args:
            directory (str): directory containing outfile (and go_outfile).
            outfile (str): output filename of the dos run.
            output_directory (str, optional): output directory. Defaults to directory.
            go_outfile (str, optional): output filename of the go run in directory,
                from which ewidth of the go run is read. Defaults to "out_go.log".
            read_go_outfile (bool, optional): read go_outfile automatically. Defaults to True.
        """
        super().__init__(directory, outfile, output_directory)
        self.go_outfile = go_outfile
        self.read_go_outfile = read_go_outfile

    def make(self, output_directory=None, yscale="log", figsize=(5, 3), ewidth_go=None):
        """make DOS from outputfile

        A vertical line is drawn at E-EF=-|ewidth_go|. ewidth_go is taken from
        the argument, then from go_outfile. No line is drawn (with a warning) if
        neither is available.

        Args:
            output_directory (str, optional): output directory of the image file. Defaults to None.
            yscale (str, optional): yscale. Defaults to "log".
            figsize (tuple, optional): figure size. Defaults to (5, 3).
            ewidth_go (float, optional): ewidth of the go run. Defaults to None.
        """
        if output_directory is None:
            output_directory = self.output_directory
        ewidth_go = resolve_ewidth_go(self.directory, ewidth_go,
                                      self.go_outfile, self.read_go_outfile)
        return _EX_plot_dos(self.directory, self.outfile,
                            output_direcotry=output_directory, yscale=yscale,
                            figsize=figsize, ewidth_go=ewidth_go)
