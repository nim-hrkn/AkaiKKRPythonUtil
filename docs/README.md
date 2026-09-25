# AkaiKKRPythonUtil ドキュメント

| 文書 | 種類 | 内容 |
|---|---|---|
| [ase_calculator_usage.md](ase_calculator_usage.md) | 使い方 | `pyakaikkr.ase.AkaiKKR`（ASE calculator）。純物質・混晶（CPA）の指定、単位、パラメータ、ポテンシャル再利用、落とし穴 |
| [testscript_usage.md](testscript_usage.md) | 使い方 | `tests/akaikkr*/testrun.py` と `testrun_ase.py`。参照ファイル、backend、環境、既知の最終桁差 |
| [dos_plot_ewidth_line.md](dos_plot_ewidth_line.md) | 使い方 | DOS / PDOS 図の E − E_F = −|ewidth_go| の線。go の ewidth を使う理由、`ewidth_go` / `go_outfile` / `read_go_outfile` の優先順 |
| [option_access_usage.md](option_access_usage.md) | 使い方 | `begin_option` を pyakaikkr から書く・検証する・読む（`OPTION_KEYS`, `normalize_option`, `read_inputcard_option`, `AkaikkrJob.get_option / get_emesh_param`、ASE の `option=` / `code=`） |
| [akaikkr_option_keys.md](akaikkr_option_keys.md) | リファレンス | AkaiKKR の `begin_option` ブロックのキー（mse, tol, ng, dex, ...）の意味と既定値、読まれる場所 |
| [ase_calculator_spec.md](ase_calculator_spec.md) | 仕様書 | ASE calculator の設計仕様（2026-09-24 実装、末尾に実装メモと後日談） |
| [testscript_ase_spec.md](testscript_ase_spec.md) | 仕様書 | テストスクリプトの ASE backend の設計仕様（同上） |
| [option_access_spec.md](option_access_spec.md) | 仕様書 | `begin_option` を pyakaikkr から参照する設計仕様（2026-09-25 実装、末尾に実装メモ）: キー辞書 `OPTION_KEYS`、inputcard の option ブロックの読み取り、出力の `optnwrt:` echo と `meshr mse ng mxl` の実効値の読み取り、書く前の検証 |
| [gaes_usage.md](gaes_usage.md) | 使い方 | GAES（`pyakaikkr.gaes`, CLI `kkr-gaes`）: 既存 dos への判定、スキームの実行、互換モード、2022.0721 の各ビルドの注意 |
| [data/atomic_core_levels_dsp.csv](data/atomic_core_levels_dsp.csv) | データ | H〜Bi の単体 fcc（a=実験原子体積、nmag、pbe、sra）に specx `go=dsp` を新規（pot.dat 無し）で掛けた初期原子ポテンシャルの core 準位 646 本（Ry / eV、core 電子数、`*`）。図 [data/atomic_core_levels_dsp.png](data/atomic_core_levels_dsp.png)、スクリプト `tests/gaes/tools/atomic_levels_dsp.py`, `plot_atomic_levels.py` |
| [core_levels_vs_dos.md](core_levels_vs_dos.md) | 解析記録 | out_go.log の E_F（`ef=` 行）と成分ごとの core 準位（`*` = valence に切替）の読み方、2019 年 RUN 7,505 系での core 準位と total DOS の semicore ピークの一致（Ge 3d, Sn 4d, Bi 5d, Pb 5d, Sc 3p, Y 4p, Zr 4p, In 4d ...）、GAES への含意 |
| [ewidth_tuning_scheme.md](ewidth_tuning_scheme.md) | 仕様書 | ewidth 自動調整スキーム GAES（Gap-Anchored Ewidth Search: DOS < threshold の連続 mesh 区間 = バンドギャップに E_F − ewidth_go が入るまで go/dos を回す。2019 年の HEA 網羅計算 run_scheme2 を一般化。微分は使わない）と GAES-Committee（複数 ewidth_go の投機的並列実行と vote）の pyakaikkr.gaes 移植仕様。移植元の誤り 12 件、PDOS 拡張、SiteComposition、type 名 40 文字制限 |

## 関連プロジェクト

- **aiida-akaikkr**（`../aiida-akaikkr/`）: AkaiKKR を AiiDA の CalcJob として流すプラグイン。構造パラメータの生成に `akaikkr_testscript.get_kkr_struc_from_cif` を、出力解析に `pyakaikkr.AkaikkrJob` を使う。文書は `aiida-akaikkr/docs/README.md`。

## ライブラリの構成

```
library/PyAkaiKKR/src/pyakaikkr/         AkaikkrJob（入力生成・出力解析）、GoGo（go/dos/spc/... の実行クラス）、
                                          Cif2Kkr、HighsymmetryKpath、AwkReader/AwkPlotter、Fmg、ase/（ASE calculator）
library/AkaiKKRTestScript/src/akaikkr_testscript/
                                          testrun_class（物質ごとのパラメータと Go* の呼び出し）、asestruc（ASE backend）、
                                          OutputAnalyzer / resultutil（参照との比較）、exeutil
tests/akaikkr, tests/akaikkr_cnd, tests/akaikkr_cpa2021v01
                                          テストセットと reference/*.json
tests/ase/                                pytest（ASE calculator、backend 一致）
```
