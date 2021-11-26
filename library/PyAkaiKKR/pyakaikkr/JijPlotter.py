# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

import pandas as pd
import matplotlib.pyplot as plt
import os


class JijPlotter:
    """load jij.csv and plot Jij
    """

    def __init__(self, directory, filename="jij.csv"):
        """initialization routine

        Args:
            directory (str): directory to run
            filename (str): csv filename to read
        """
        self.directory = directory
        filepath = os.path.join(directory, filename)
        df = pd.read_csv(filepath)
        self.df = df

    def plot_typepair(self, a=1.0, marker="o", output_directory=None,
                      figsize=(5, 3)):
        """plot Jij of all type pairs
        save all combinations of type pairs

        Args:
            a (float, optional): a length. Defaults to 1.0.
            marker (str, optional): matplotlib.plot() maker. Defaults to "o".
            output_directory (str, optional): output directory.
            figsize (tuple, optional): figure size. Defaults to (5, 3).
        """
        directory = self.directory
        if output_directory is None:
            output_directory = directory
        df = self.df.copy()

        xlabel = "distance"
        ylabel = "J_ij(meV)"
        xlabel_fig = "$R$"
        ylabel_fig = "$J_{ij}$ (meV)"

        # plot range
        values = df[xlabel].astype(float).values*a
        xlim = (values.min(), values.max())
        dx = (xlim[1]-xlim[0])*0.05
        xlim = (xlim[0]-dx, xlim[1]+dx)
        values = df[ylabel].values
        ylim = (values.min(), values.max())
        dy = (ylim[1]-ylim[0])*0.05
        ylim = (ylim[0]-dy, ylim[1]+dy)

        # make type pair
        type1 = df["type1"]
        type2 = df["type2"]
        type_pair_list = []
        for t1, t2 in zip(type1, type2):
            s = "-".join([t1, t2])
            type_pair_list.append(s)
        df["typepair"] = type_pair_list
        uniq_type_pair = list(set(type_pair_list))

        # make figures
        for pair in uniq_type_pair:
            _df = df.query("typepair=='{}'".format(pair))

            distance = _df[xlabel]*a
            Jij = _df[ylabel]

            fig, ax = plt.subplots(figsize=figsize)
            ax.plot(distance, Jij, linestyle="-", marker=marker)
            ax.axhline(y=0, linestyle="--", linewidth=1)
            ax.set_xlabel(xlabel_fig)
            ax.set_ylabel(ylabel_fig)
            ax.set_title(pair)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            if True:
                imgfile = "Jij_{}.png".format(pair)
                imgpath = os.path.join(output_directory,
                                       imgfile)
                fig.tight_layout()
                fig.savefig(imgpath)
                print("  saved to", imgpath)
                fig.clf()
                plt.close(fig)

    def plot_comppair(self, type1, type2, typeofsite, a=1.0,
                      marker="o",  output_directory=None,
                      figsize=(5, 3)):
        """plot Jij of specified type1 and type2
        save png images of all the (comp1,comp2) combinations

        Args:
            type1 (str): type1 name
            type2 (str): type2 name
            typeofsite (dict): typeofsite of AkaikkrJob
            a (float, optional): a scale of length. Defaults to 1.0.
            marker (str, optional): matplotlib.plot() marker. Defaults to "o".
            output_directory (str, optional): output directory.
            figsize (tuple, optional): figure size. Defaults to (5, 3).
        """
        directory = self.directory
        if output_directory is None:
            output_directory = directory
        df = self.df.query("type1=='{}' and type2=='{}'".format(
            type1, type2)).reset_index(drop=True)
        xlabel = "distance"
        ylabel = "J_ij(meV)"
        xlabel_fig = "$R$"
        ylabel_fig = "$J_{ij}$ (meV)"

        # plot range
        values = df[xlabel].astype(float).values*a
        xlim = (values.min(), values.max())
        dx = (xlim[1]-xlim[0])*0.05
        xlim = (xlim[0]-dx, xlim[1]+dx)
        values = df[ylabel].values
        ylim = (values.min(), values.max())
        dy = (ylim[1]-ylim[0])*0.05
        ylim = (ylim[0]-dy, ylim[1]+dy)

        # make comp pair
        comp1 = df["comp1"]
        comp2 = df["comp2"]
        type_pair_list = []
        for t1, t2 in zip(comp1, comp2):
            s = "-".join([str(t1), str(t2)])
            type_pair_list.append(s)
        df_comppair = pd.DataFrame({"comppair": type_pair_list})
        df = pd.concat([df, df_comppair], axis=1)
#        df["comppair"] = type_pair_list # warning occurs
        uniq_type_pair = list(set(type_pair_list))

        shortname_dic = {}
        for component in typeofsite:
            type = component["type"]
            shortname_dic[type] = component["comp_shortname"]

        # make figures

        for pairname in uniq_type_pair:
            s = pairname.split("-")
            comp1, comp2 = int(s[0])-1, int(s[1])-1
            comp1name = shortname_dic[type1][comp1]
            comp2name = shortname_dic[type2][comp2]
            label = "{}-{}".format(comp1name, comp2name)
            _df = df.query("comppair=='{}'".format(pairname))

            distance = _df[xlabel]*a
            Jij = _df[ylabel]

            fig, ax = plt.subplots(figsize=figsize)
            ax.plot(distance, Jij, linestyle="-", marker=marker, label=label)
            ax.axhline(y=0, linewidth=1, linestyle="--")
            ax.set_xlabel(xlabel_fig)
            ax.set_ylabel(ylabel_fig)
            ax.set_title(label)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            if True:
                imgfile = "Jij_{}_{}_{}.png".format(type1, type2, pairname)
                imgpath = os.path.join(output_directory, imgfile)
                fig.tight_layout()
                fig.savefig(imgpath)
                fig.clf()
                plt.close(fig)
                print("  saved to", imgpath)

# todo:
# (type, comp) pair function must be made.
