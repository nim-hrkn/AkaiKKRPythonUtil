# AkaiKKRPythonUtil ドキュメント

| 文書 | 種類 | 内容 |
|---|---|---|
| [ase_calculator_usage.md](ase_calculator_usage.md) | 使い方 | `pyakaikkr.ase.AkaiKKR`（ASE calculator）。純物質・混晶（CPA）の指定、単位、パラメータ、ポテンシャル再利用、落とし穴 |
| [testscript_usage.md](testscript_usage.md) | 使い方 | `tests/akaikkr*/testrun.py` と `testrun_ase.py`。参照ファイル、backend、環境、既知の最終桁差 |
| [akaikkr_option_keys.md](akaikkr_option_keys.md) | リファレンス | AkaiKKR の `begin_option` ブロックのキー（mse, tol, ng, dex, ...）の意味と既定値、読まれる場所 |
| [ase_calculator_spec.md](ase_calculator_spec.md) | 仕様書 | ASE calculator の設計仕様（2026-09-24 実装、末尾に実装メモと後日談） |
| [testscript_ase_spec.md](testscript_ase_spec.md) | 仕様書 | テストスクリプトの ASE backend の設計仕様（同上） |

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
