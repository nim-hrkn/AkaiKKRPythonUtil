# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""run the tests with the cif (pymatgen) backend.

usage: python testrun.py <program_path> [--set <name>] [--create_ref] [--compiler ifort]
"""
from akaikkr_testscript import all_go
from testsets import make_exe

if __name__ == "__main__":
    akaikkr_exe = "specx"
    fmg_exe = "fmg"
    exe_dic = make_exe()
    all_go(akaikkr_exe, fmg_exe=fmg_exe, exe_dic=exe_dic, displc=False)
