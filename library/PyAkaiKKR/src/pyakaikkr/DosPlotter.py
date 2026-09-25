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

# style of the vertical line at E - EF = -|ewidth_go|
_EWIDTH_LINE_STYLE = {"color": "tab:red", "lw": 1.0, "ls": "-."}
_EWIDTH_LINE_LABEL = "$-$ewidth (go)"
_DEFAULT_GO_OUTFILE = "out_go.log"


def _mark_ewidth_go(ax, energy, ewidth_go):
    """draw a dash-dot vertical line at E - EF = -|ewidth_go|.

    AkaiKKR integrates the charge over [EF - ewidth, EF] in the go run,
    so the line shows the bottom of the SCF energy contour. The dos run
    itself uses another ewidth to define its energy mesh
    [EF - ref*ewidth, EF + (1 - ref)*ewidth], so the ewidth of the go run,
    not of the dos run, must be given.

    If the line lies below the energy mesh, the x range is widened to show it.

    Args:
        ax (matplotlib.axes.Axes): axes to draw on.
        energy ([float]): energy mesh (E - EF) of the DOS.
        ewidth_go (float, None): ewidth of the go run. Nothing is drawn if None.

    Returns:
        bool: True if the line is drawn.
    """
    if ewidth_go is None:
        return False
    ebtm = -abs(float(ewidth_go))
    ax.axvline(ebtm, label=_EWIDTH_LINE_LABEL, **_EWIDTH_LINE_STYLE)
    energy = np.asarray(energy, dtype=float)
    lo, hi = float(energy.min()), float(energy.max())
    if ebtm < lo:
        margin = 0.02 * (hi - ebtm)
        ax.set_xlim(ebtm - margin, hi + margin)
    return True


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
    if len(dos_block) > 1:  # mag
        dos_up, dos_dn = dos_block[0], dos_block[1]
        _figsize = (figsize[0]*2, figsize[1])
        fig, axes = plt.subplots(1, 2, figsize=_figsize)
        ax = axes[0]
        ax.plot(energy, dos_up)
        ax.set_title("up")
        ax.set_yscale(yscale)
        ax.set_xlabel("E(Ry)-EF")
        ax.set_ylabel("DOS")
        if _mark_ewidth_go(ax, energy, ewidth_go):
            ax.legend(frameon=False, fontsize=8)
        ax = axes[1]
        ax.plot(energy, dos_dn)
        ax.set_title("dn")
        ax.set_yscale(yscale)
        ax.set_xlabel("E(Ry)-EF")
        ax.set_ylabel("DOS")
        if _mark_ewidth_go(ax, energy, ewidth_go):
            ax.legend(frameon=False, fontsize=8)

    else:  # monmag
        dos_up = dos_block[0]
        _figsize = figsize
        fig, ax = plt.subplots(figsize=_figsize)
        ax.plot(energy, dos_up)
        ax.set_yscale(yscale)
        ax.set_xlabel("E(Ry)-EF")
        ax.set_ylabel("DOS")
        if _mark_ewidth_go(ax, energy, ewidth_go):
            ax.legend(frameon=False, fontsize=8)

    if output_direcotry is None:
        output_direcotry = "."
    fig.tight_layout()
    imgfile = "dos.png"
    imgfilepath = os.path.join(output_direcotry, imgfile)
    fig.savefig(imgfilepath)
    print("  saved to", imgfilepath)
    # fig.show()
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
    l_label = ["s", "p", "d", "f", "g", "h", "i", "j", "k", "l", "m"]

    serial_site = []
    for site in typeofsite:
        for shortname in site["comp_shortname"]:
            serial_site.append(shortname)

    energy = np.array(energy)

    if output_direcotry is None:
        output_direcotry = "."

    imgfilepaths = []
    if len(pdos_block) > 1:  # up and down
        pdos_up, pdos_dn = pdos_block[0], pdos_block[1]

        _figsize = (figsize[0]*2, figsize[1])
        for icmp, (up_atom, dn_atom, title) in enumerate(zip(pdos_up, pdos_dn, serial_site)):
            fig, axes = plt.subplots(1, 2, figsize=_figsize)

            up_atom = np.array(up_atom)
            dn_atom = np.array(dn_atom)
            ax = axes[0]
            for l in range(up_atom.shape[1]):
                ax.plot(energy, up_atom[:, l], label=l_label[l])
            _mark_ewidth_go(ax, energy, ewidth_go)
            ax.legend()
            ax.set_title("up")
            ax.set_yscale(yscale)
            ax.set_xlabel("E(Ry)-EF")
            ax.set_ylabel("PDOS")
            ax = axes[1]
            for l in range(dn_atom.shape[1]):
                ax.plot(energy, dn_atom[:, l], label=l_label[l])
            _mark_ewidth_go(ax, energy, ewidth_go)
            ax.legend()
            ax.set_title("dn")
            ax.set_yscale(yscale)
            ax.set_xlabel("E(Ry)-EF")
            ax.set_ylabel("PDOS")
            fig.suptitle(title)
            fig.tight_layout()
            imgfile = "pdos_{}.png".format(icmp)
            imgfilepath = os.path.join(output_direcotry, imgfile)
            fig.savefig(imgfilepath)
            print("  saved to", imgfilepath)
            imgfilepaths.append(imgfilepath)
            # fig.show()
            fig.clf()
            plt.close(fig)
        print()

    else:  # up only
        pdos_up = pdos_block[0]

        _figsize = (figsize[0], figsize[1])
        for icmp, (up_atom, title) in enumerate(zip(pdos_up, serial_site)):
            fig, ax = plt.subplots(1, 1, figsize=_figsize)

            up_atom = np.array(up_atom)
            for l in range(up_atom.shape[1]):
                ax.plot(energy, up_atom[:, l], label=l_label[l])
            _mark_ewidth_go(ax, energy, ewidth_go)
            ax.legend()
            ax.set_yscale(yscale)
            ax.set_xlabel("E(Ry)-EF")
            ax.set_ylabel("PDOS")
            ax.set_title(title)
            fig.tight_layout()
            imgfile = "pdos_{}.png".format(icmp)
            imgfilepath = os.path.join(output_direcotry, imgfile)
            fig.savefig(imgfilepath)
            print("  saved to", imgfilepath)
            imgfilepaths.append(imgfilepath)
            # fig.show()
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
