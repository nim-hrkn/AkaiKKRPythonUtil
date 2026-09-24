# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""run the tests with the ase backend.

The structures are read with ase.io.read, reduced to the primitive cell and passed
to AkaiKKR as brvtyp=aux. The reference is reference/<compiler>_ase.json and the
result is saved to result_ase.json. Make the reference first with --create_ref.

usage: python testrun_ase.py <program_path> [--set <name>] [--create_ref] [--compiler ifort]
"""
from akaikkr_testscript import all_go
from testsets import make_exe

if __name__ == "__main__":
    akaikkr_exe = "specx"
    fmg_exe = "fmg"
    exe_dic = make_exe()
    all_go(akaikkr_exe, fmg_exe=fmg_exe, exe_dic=exe_dic, displc=False, backend="ase")
