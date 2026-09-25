"""Atomic (initial-potential) core levels of the elements from specx go=dsp on a fresh potential.

For every element Z=1..zmax a single-site fcc cell (a=1000000: experimental atomic volume,
nmag, sra, sdftyp) is run with go=dsp and no pot.dat: specx generates the initial atomic
potential, does one iteration and prints the core configuration and the core levels.
Output CSV columns: Z, element, orbital, energy_Ry, energy_eV, core_electrons, star,
ef_up, ef_dn, a_bohr, status. `star` marks levels specx moved to valence (above the
contour bottom of ewidth); ef_up/ef_dn are the dummy initial values (0.8 / 0.6) in this mode.

Usage: python atomic_levels_dsp.py --exe /path/to/specx -o atomic_core_levels_dsp.csv [--zmax 83] [--ewidth 1.2]
"""
import argparse
import csv
import os
import re
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor

from pymatgen.core.periodic_table import Element

RY_EV = 13.605693

CARD = """#--- go potentialfile
dsp pot.dat
#--- brvtyp a c/a b/a alpha beta gamma
fcc {a} 1.0 1.0 90 90 90
#--- edelt ewidth reltyp sdftyp magtyp record
0.001 {ewidth} sra {sdftyp} nmag 2nd
#--- outtyp bzqlty maxitr pmix
quit 4 1 0.02
#--- ntyp
1
#--- type ncmp
  X 1
#- rmt field mxl
1.0 0.0 {mxl}
#- anclr conc
{z} 100
#--- natm
1
#--- atmicx atmtyp
0.0a 0.0b 0.0c X
"""
ORBITALS = ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "5s", "5p", "4f", "5d", "6s", "6p", "5f", "6d", "7s"]


def run_element(z, exe, ewidth, sdftyp, threads, lattice="1000000"):
    """one dsp run; retried with the MJW volume table (a=0) then a=9 bohr when the
    experimental-volume table has no entry, and with ewidth 0.8 / 2.0 when reconf stops."""
    rows = _run_element(z, exe, ewidth, sdftyp, threads, lattice)
    st = rows[0]["status"]
    if st == "ok" and (rows[0]["a_bohr"] or 0) <= 0 and lattice == "1000000":
        rows = _run_element(z, exe, ewidth, sdftyp, threads, "0")
        if (rows[0]["a_bohr"] or 0) <= 0 or rows[0]["status"] != "ok":
            rows = _run_element(z, exe, ewidth, sdftyp, threads, "9.0")
    elif "reconf" in st:
        for ew in (0.8, 2.0, 0.5):
            rows = _run_element(z, exe, ew, sdftyp, threads, lattice)
            if rows[0]["status"] == "ok":
                break
    return rows


def _run_element(z, exe, ewidth, sdftyp, threads, lattice):
    mxl = 3 if z >= 57 else 2
    with tempfile.TemporaryDirectory() as d:
        with open(os.path.join(d, "inputcard"), "w") as f:
            f.write(CARD.format(z=z, ewidth=ewidth, sdftyp=sdftyp, mxl=mxl, a=lattice))
        env = dict(os.environ, OMP_NUM_THREADS=str(threads))
        with open(os.path.join(d, "inputcard")) as fin, open(os.path.join(d, "out.log"), "w") as fout:
            rc = subprocess.call([exe], stdin=fin, stdout=fout, stderr=subprocess.STDOUT, cwd=d, env=env)
        txt = open(os.path.join(d, "out.log")).read()
    el = str(Element.from_Z(z))
    rows = []
    err = re.search(r"\*\*\*err in (\S+)\.\.\.(.*)", txt)
    status = "ok" if rc == 0 and err is None and ("core level" in txt or z == 1) else ("err: " + (err.group(0).strip() if err else "rc=%d" % rc))
    m = re.search(r"core configuration for Z= *%d\n.*\n\s*up\s+(.*)\n" % z, txt)
    config = dict(zip(ORBITALS, m.group(1).split())) if m else {}
    m = re.search(r"\n\s+ef=\s*([-\d.]+)\s+([-\d.]+)", txt)
    ef_up, ef_dn = (float(m.group(1)), float(m.group(2))) if m else (None, None)
    m = re.search(r"bravais=\S+\s+a=\s*([\d.]+)", txt)
    a = float(m.group(1)) if m else None
    blk = re.search(r"core level  \(spin up  \)\n(.*?)(?:core level  \(spin down\)|\n\s*\n)", txt, re.S)
    if blk:
        for v, orb, star in re.findall(r"(-?\d+\.\d+) Ry\((\d[spdf])\)(\*?)", blk.group(1)):
            e = float(v)
            rows.append(dict(Z=z, element=el, orbital=orb, energy_Ry=e, energy_eV=round(e * RY_EV, 3),
                             core_electrons=config.get(orb, ""), star=int(bool(star)), ef_up=ef_up, ef_dn=ef_dn,
                             a_bohr=a, ewidth=ewidth, status=status))
    if not rows:
        rows.append(dict(Z=z, element=el, orbital="", energy_Ry="", energy_eV="", core_electrons="", star="",
                         ef_up=ef_up, ef_dn=ef_dn, a_bohr=a, ewidth=ewidth, status=status))
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--exe", required=True)
    p.add_argument("-o", "--output", default="atomic_core_levels_dsp.csv")
    p.add_argument("--zmin", type=int, default=1)
    p.add_argument("--zmax", type=int, default=83)
    p.add_argument("--ewidth", type=float, default=1.2)
    p.add_argument("--sdftyp", default="pbe")
    p.add_argument("--parallel", type=int, default=8)
    p.add_argument("--threads", type=int, default=2)
    a = p.parse_args()
    with ThreadPoolExecutor(a.parallel) as ex:
        results = list(ex.map(lambda z: run_element(z, a.exe, a.ewidth, a.sdftyp, a.threads), range(a.zmin, a.zmax + 1)))
    fields = ["Z", "element", "orbital", "energy_Ry", "energy_eV", "core_electrons", "star", "ef_up", "ef_dn", "a_bohr", "ewidth", "status"]
    with open(a.output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for rows in results:
            for r in rows:
                w.writerow(r)
    bad = [rows[0] for rows in results if rows[0]["status"] != "ok"]
    print("wrote", a.output, "elements:", len(results), "levels:", sum(len(r) for r in results if r[0]["status"] == "ok"),
          "failed:", [(r["Z"], r["element"], r["status"]) for r in bad])


if __name__ == "__main__":
    main()
