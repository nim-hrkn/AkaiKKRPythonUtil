# HTML レポート（`pyakaikkr.report`、`kkr-report`）の使い方

仕様: aiida-akaikkr/docs/report_spec.md（実装 2026-09-26）。1 つの計算（go とその後続の dos / spc / j / cnd）の式、構造の出典、空間群、格子、計算条件、SCF の結果、成分ごとのモーメントと電荷、図（DOS、PDOS、A(ω,k)、J_ij、GAES）、輸送の値を **1 枚の自己完結な HTML** にする。図は SVG（inline、A(ω,k) のメッシュだけ rasterize）と PNG（data URI、`--embed png`）。言語は英語と日本語。

## ファイルから

```
kkr-report tests/akaikkr/FeRh05Pt05 --lang ja -o report.html [--embed png] [--cif FeRh0.5Pt0.5.cif] [--poscar POSCAR] [--figure-dir figs]
```

- ディレクトリの `out_go.log`（必須）、`out_dos.log`、`out_spc*.log` + `pot.dat_up.spc` / `pot.dat_dn.spc` + `klabel.json`、`out_j*.log`、`out_cnd.log` を見つけたものだけ使う（`--go / --dos / --spc / --j / --cnd` で名前を指定できる）。
- 標準出力に JSON（`html` のパスと `summary`: 式、空間群、全エネルギー、モーメント、Tc、cnd、図の名前）。
- `--figure-dir` で各図の PNG と SVG も別ファイルに書く。

例（testrun の出力、docs/data に同梱）: [report_FeRh05Pt05_en.html](data/report_FeRh05Pt05_en.html)、[report_FeRh05Pt05_ja.html](data/report_FeRh05Pt05_ja.html)、[report_SmCo5_oc_ja.html](data/report_SmCo5_oc_ja.html)（Sm は f まであり PDOS の l 数が type で違う例）。

## Python から

```python
from pyakaikkr.report import report_from_directory, write_report, render_html
data = report_from_directory("tests/akaikkr/FeRh05Pt05", structure_files={"cif": "FeRh0.5Pt0.5.cif"}, lang="ja")
data.summary()                       # 式、空間群、全エネルギー、全モーメント、Tc、cnd、図の名前
write_report(data, "report.html", lang="ja", embed="svg", figure_dir="figs")
html_text = render_html(data, lang="en")
```

材料を自分で組む場合は `ReportData` / `Component` を作り、`add_dos_figure(data, energy, dos, ewidth_go=)`、`add_pdos_figure`、`add_awk_figures`、`add_jij_figures`、`add_gaes_figure` で図を足す（どれも配列を受け、`pyakaikkr.plot` で描く）。aiida-akaikkr の `akaikkr-aiida report --pk` / MCP `kkr_report` はこの経路でノードから作る。

## 構成

1. 構造: 組成式（約分、CPA は `Fe(Rh0.5Pt0.5)`）と胞全体の式、出典（CIF / POSCAR / preset / comp）、空間群（pymatgen SpacegroupAnalyzer、symprec 1e-3）、格子（brvtyp、a、c/a、角）、原子数。
2. 計算条件（code、build、magtyp、sdftyp、reltyp、ewidth、edelt、bzqlty、maxitr、pmix。ファイルから分かるものだけ）。
3. SCF の結果（収束、反復数、log10 rms、E_F、全エネルギー、全モーメント）。
4. 成分の表（type、元素、Z、濃度、スピン / 軌道モーメント、電荷）。Tc（J_ij があれば）。
5. 図と 1 行の説明。
6. 輸送（cnd の抵抗率・伝導度）。GAES（状態、ewidth、試した値、ギャップ）。
7. 来歴（段階ごとのファイル、または AiiDA の pk と作業ディレクトリ）。

テスト: `tests/plot/test_report.py`（合成配列と、tests/akaikkr の FeRh05Pt05 / SmCo5_oc / Cu の出力があればそれも）。
