# AkaiKKRPythonUtil
Python utilities for AkaiKKR

- `library/PyAkaiKKR` (`pyakaikkr`): input-card generation, output parsing, `go/dos/spc/j/tc/fsm/cnd` runner classes, A(w,k)/DOS plotters, an [ASE calculator](#ase-interface).
- `library/AkaiKKRTestScript` (`akaikkr_testscript`): the per-material test definitions and result comparison used by `tests/`.
- `tests/akaikkr`, `tests/akaikkr_cnd`, `tests/akaikkr_cpa2021v01`: regression tests against `reference/*.json`.
- `docs/`: documentation index in [docs/README.md](docs/README.md) (usage of the ASE calculator and the test scripts, the `begin_option` keys, design specs).
- AiiDA: the separate plugin [aiida-akaikkr](../aiida-akaikkr) runs the same test set as AiiDA CalcJobs.

# License

 Copyright (c) 2021-2023 AkaiKKRteam.
 Distributed under the terms of the Apache License, Version 2.0.



# installation
{PREFIX} is the directory you installed AkaiKKRPythonUtil.

install PyAkaiKKR by
```
cd {PREFIX}/library/PyAkaiKKR
pip install .
```

If you use {PREFIX}/tests, also install AkaiKKRTestScript by
```
cd {PREFIX}/library/AkaiKKRTestScript
pip install .
```

## for the CPA2021V01 user

{AKAIKKRCPA2021V01} is the directory where you installed CPA2021V01.
{AKAIKKRCPA2021V01} can be an absolute or a relative path.

1. Please apply patch to akaikkr_cpa2021v01.
```
$ cd {AKAIKKRCPA2021V01}
$ patch -p1 < {PREFIX}/tests/akaikkr_cpa2021v01/akaikkr_cpa2021v01.patch
```

The patched cpa2021v01 has 
- 'geom' mode
- Site labels exteded to 40 characters.
- additional information on the "dispersion" files.

2. make specx and fmg in CPA2021V01.
```
$ make specx fmg
```

# test script

- It tests AkaiKKR package.
- It generates AkaiKKR input file from cif files.

 
For example, if you use CPA2021V01, (Note that the directory name specified by {PREFIX} must be .../akaikkr_cpa2021v01.)
```
$ cd {PREFIX}/tests/akaikkr_cpa2021v01
$ python testrun.py <program_path> [--set <name>] [--create_ref] [--compiler ifort]
```
`<program_path>` is the directory that contains `akaikkr/specx`, `akaikkr_cnd/specx`, ...
See [docs/testscript_usage.md](docs/testscript_usage.md) for the options, the reference files and the known last-digit differences.

The following display appears at the end of the execution.
```
SHORT SUMMARY
             dos fsm go gofmg spc31
AlMnFeCo_bcc   O   O  O     O     -
Co             O   O  O           -
Co2MnSi        O   O  O           -
Cu             O      O           -
Fe             O   O  O           -
FeB195         O      O           -
FeRh05Pt05     O   O  O           -
Fe_lmd         O      O           -
GaAs           O      O           -
Ni             O   O  O           -
NiFe           O   O  O           -
SmCo5_noc             O            
SmCo5_oc       O   O  O           -
O: passed. X: failed, -: no reference
```

As an example, the block spectra $A(w,k)$ of NiFe (FCC $\mathrm{Ni}_{0.9}\mathrm{Fe}_{0.1}$) is generated under the NiFe directory as follows.
![](https://github.com/nim-hrkn/AkaiKKRPythonUtil/blob/cpa2021v01_supported/fig/NiFe_Awk_all.png?raw=true)

# ASE interface

`pyakaikkr.ase` runs AkaiKKR through an [ASE](https://wiki.fysik.dtu.dk/ase/) calculator (energy, magmom, magmoms; no forces).
Install ase by `pip install ase` (or `pip install {PREFIX}/library/PyAkaiKKR[ase]`).

```python
from ase.build import bulk
from pyakaikkr.ase import AkaiKKR, set_occupancy

cu = bulk("Cu", "fcc", a=3.615)
cu.calc = AkaiKKR(directory="cu", command="/path/to/specx < PREFIX.in > PREFIX.out",
                  magtyp="nmag")
print(cu.get_potential_energy())          # eV per cell

# CPA alloy: partial occupancies in the ASE convention (atoms.info["occupancy"])
nife = set_occupancy(bulk("Ni", "fcc", a=3.571), {0: {"Fe": 0.1, "Ni": 0.9}})
nife.calc = AkaiKKR(directory="nife", command="/path/to/specx < PREFIX.in > PREFIX.out",
                    sdftyp="pbeasa", bzqlty=8)
print(nife.get_potential_energy(), nife.get_magnetic_moment())
print(nife.calc.results["component_moments"])   # moment of each component of each type
```

The structure is passed as `brvtyp=aux` with `a=|cell[0]|`. Sites are grouped into AkaiKKR types
when they have the same occupancy and are symmetrically equivalent (spglib).
CIF files with partial occupancies read by `ase.io.read` work as they are.
Partial occupancies cannot come through pymatgen's `AseAtomsAdaptor.get_atoms` (it rejects disordered structures); build them with `set_occupancy` or read the CIF with `ase.io.read`.
See [docs/ase_calculator_usage.md](docs/ase_calculator_usage.md) (usage) and [docs/ase_calculator_spec.md](docs/ase_calculator_spec.md) (design).

The test script can use the ASE backend: run `python testrun_ase.py <program_path> --create_ref`
in `tests/akaikkr` to make `reference/ifort_ase.json`, then `python testrun_ase.py <program_path>`.
`testrun.py` (pymatgen backend) is unchanged. See [docs/testscript_usage.md](docs/testscript_usage.md) and [docs/testscript_ase_spec.md](docs/testscript_ase_spec.md).

# AiiDA

[aiida-akaikkr](../aiida-akaikkr) wraps specx as AiiDA CalcJobs (go / fsm / dos / j3.0 / tc / spc31 / cnd). Its `example/run_examples.py` uses the same `_<material>_common_param` definitions as the test script and reproduces `tests/*/reference/ifort.json`; see `aiida-akaikkr/docs/`.

# BUG
- TEST FAILED is always shown at the end of testrun.py.
- Awk\_both.png is generated at the top directory of testrun.py.
