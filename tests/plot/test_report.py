"""pyakaikkr.report: HTML report from arrays (no specx) and from the test outputs when they exist."""
import os

import numpy as np
import pytest

from pyakaikkr.plot import figure_to_svg, figure_to_png, data_uri, save_figure
from pyakaikkr.report import (ReportData, Component, render_html, add_dos_figure, add_awk_figures, add_jij_figures,
                              add_gaes_figure, formula_of, symmetry_of, components_from_outputs, report_from_directory, T, LANGS)

AKAIKKR_TESTS = os.path.join(os.path.dirname(__file__), "..", "akaikkr")


def _data():
    d = ReportData(title="FeRhPt test", formula="Fe(Rh0.5Pt0.5)", formula_full="Fe1 Rh0.5 Pt0.5",
                   structure_source={"cif": "FeRh0.5Pt0.5.cif", "preset": "FeRh05Pt05"},
                   spacegroup={"number": 123, "symbol": "P4/mmm", "symprec": 1e-3, "crystal_system": "tetragonal"},
                   lattice={"brvtyp": "aux", "a_bohr": 5.2, "a_angstrom": 2.75}, natm=2,
                   calc={"code": "specx", "magtyp": "mag", "ewidth": 1.0, "edelt": 1e-3},
                   scf={"converged": True, "n_iter": 40, "rms_last": -6.1, "fermi_level": 0.56, "total_energy_Ry": -25743.01, "total_moment": 3.45},
                   components=[Component("Fe_0", "Fe", 26, 1.0, 2.98, 0.0, 26.2), Component("Rh0.5Pt0.5_1", "Rh", 45, 0.5, 0.61, None, 44.8)],
                   tc=217.1, cnd={"resistivity": 1.2e-6, "conductivity": [1.0, 2.0]}, provenance=[{"step": "go", "pk": 12, "where": "x"}])
    return d


def test_render_html_both_languages_with_svg_and_png_embedding():
    d = _data()
    e = np.linspace(-1, 0.5, 30)
    add_dos_figure(d, e, np.vstack([np.ones(30), np.ones(30)]), ewidth_go=1.0, lang="en")
    assert d.figures[0].svg.startswith("<svg") and d.figures[0].png[:4] == b"\x89PNG"
    for lang in LANGS:
        page = render_html(d, lang=lang)
        assert '<html lang="{}"'.format(lang) in page and T[lang]["formula"] in page and "P4/mmm" in page
        assert "Fe(Rh0.5Pt0.5)" in page and "217.1" in page and "<svg" in page and "Rh0.5Pt0.5_1" in page
    page_ja = render_html(d, lang="ja")
    assert "全エネルギー" in page_ja and "キュリー温度" in page_ja
    page_png = render_html(d, lang="en", embed="png")
    assert "data:image/png;base64," in page_png and "<svg" not in page_png
    with pytest.raises(ValueError):
        render_html(d, lang="fr")
    assert d.summary()["formula"] == "Fe(Rh0.5Pt0.5)" and d.summary()["figures"] == ["dos"]


def test_awk_mesh_is_rasterized_in_svg_and_jij_gaes_figures_are_added():
    d = _data()
    kdist = np.linspace(0, 2, 20)
    en = np.linspace(-1, 1, 25)
    awk = np.abs(np.random.default_rng(0).random((20, 25)))
    add_awk_figures(d, {"up": awk, "dn": None}, kdist, en, [0, 10, 19], ["G", "X", "L"])
    assert [f.name for f in d.figures] == ["awk_up"] and d.figures[0].svg.count("<image") >= 1   # mesh rasterized
    import pandas as pd
    df = pd.DataFrame({"type1": ["A"] * 3, "type2": ["A"] * 3, "comp1": [1, 1, 1], "comp2": [1, 1, 1],
                       "distance": [1.0, 1.4, 2.0], "J_ij(meV)": [5.0, 1.0, -0.5]})
    add_jij_figures(d, df, {"A": ["Fe"]}, tc=100.0)
    e = np.linspace(-2.2, 0.7, 100)
    dos = np.full(100, 5e-4)
    dos[e > -0.5] = 5.0
    add_gaes_figure(d, e, dos, 1.2, coarse=[(-2.2, -0.5)], fine=[(-2.2, -0.7)], bounds=(1.0, 2.0), levels={"Bi6s": [-0.9, True]},
                    natm=1, dosth=2e-2, dosth2=1e-3, final=1.4)
    assert [f.name for f in d.figures] == ["awk_up", "jij_A-A", "gaes"]
    page = render_html(d, lang="ja")
    assert page.count("<figure") == 3 and "GAES" in page


def test_figure_export_helpers(tmp_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.pcolormesh(np.random.default_rng(1).random((5, 6)))
    svg = figure_to_svg(fig)
    assert svg.startswith("<svg") and "<image" in svg
    assert data_uri(figure_to_png(fig)).startswith("data:image/png;base64,")
    paths = save_figure(fig, str(tmp_path / "f"))
    assert set(paths) == {"png", "svg"} and all(os.path.getsize(p) > 0 for p in paths.values())
    plt.close(fig)


def test_formula_symmetry_and_components():
    from pymatgen.core import Structure, Lattice
    s = Structure(Lattice.cubic(3.0), [{"Fe": 1.0}, {"Rh": 0.5, "Pt": 0.5}], [[0, 0, 0], [0.5, 0.5, 0.5]])
    assert formula_of(s)[0] == "Fe(Rh0.5Pt0.5)"
    sg = symmetry_of(s)
    assert sg and sg["number"] == 221 and sg["symbol"] == "Pm-3m"
    tos = [{"type": "HEA", "component": [{"anclr": 29.0, "conc": 0.5}, {"anclr": 30.0, "conc": 0.5}], "comp_shortname": ["HEA_Cu", "HEA_Zn"]}]
    comps = components_from_outputs(tos, [0.1, 0.2], None, [28.0, 29.0])
    assert [(c.element, c.conc, c.spin_moment, c.charge) for c in comps] == [("Cu", 0.5, 0.1, 28.0), ("Zn", 0.5, 0.2, 29.0)]


@pytest.mark.parametrize("name,expect", [("FeRh05Pt05", {"formula": "Fe(Rh0.5Pt0.5)", "sg": "P4/mmm", "figs": {"dos", "pdos", "awk_up", "awk_dn"}, "jij": 3, "tc": True}),
                                         ("SmCo5_oc", {"formula": "SmCo5", "sg": "P6/mmm", "figs": {"dos", "pdos", "awk_up", "awk_dn"}, "jij": 6, "tc": True}),
                                         ("Cu", {"formula": "Cu", "sg": "Fm-3m", "figs": {"dos", "pdos", "awk_up"}, "jij": 0, "tc": False})])
def test_report_from_test_outputs(name, expect, tmp_path):
    """the testrun outputs (tests/akaikkr/<name>/) when present: go + dos + spc31 (+ j3.0). SmCo5_oc has a type
    with f states (unequal l per type, NaN-padded PDOS)."""
    d = os.path.join(AKAIKKR_TESTS, name)
    if not os.path.isfile(os.path.join(d, "out_go.log")):
        pytest.skip("no testrun output in " + d)
    data = report_from_directory(d, structure_files={"cif": name + ".cif"}, lang="ja")
    assert data.formula == expect["formula"] and data.spacegroup["symbol"] == expect["sg"] and data.scf["converged"]
    names = [f.name for f in data.figures]
    assert expect["figs"] <= set(names) and sum(n.startswith("jij_") for n in names) == expect["jij"]
    assert (data.tc is not None) == expect["tc"] and data.components and data.scf["total_energy_Ry"] < 0
    out = tmp_path / (name + ".html")
    from pyakaikkr.report import write_report
    write_report(data, str(out), lang="ja", figure_dir=str(tmp_path / "figs"))
    page = out.read_text(encoding="utf-8")
    assert expect["formula"] in page and "<svg" in page and (tmp_path / "figs" / "dos.svg").exists()
