# テストスクリプトの使い方

`tests/testrun/akaikkr`（標準ビルド）、`tests/testrun/akaikkr_cnd`（伝導度ビルド）、`tests/testrun/akaikkr_cpa2021v01`（CPA2021V01）にそれぞれ `testrun.py` があります。物質ごとのパラメータは `akaikkr_testscript.testrun_class` の `_<物質>_common_param` と `<物質>_<mode>` 関数、実行は `pyakaikkr.GoGo` の `GoGo / GoDos / GoSpc / Goj30 / GoTc / GoFsm / GoCnd / GoFmg` が担います。

## 1. 実行

```bash
conda activate akaikkr
cd tests/testrun/akaikkr
export OMP_NUM_THREADS=24
python testrun.py /path/to/AkaiKKRprogram.2022.0721.ifort [--set Cu] [--create_ref] [--compiler ifort]
python testrun_ase.py /path/to/AkaiKKRprogram.2022.0721.ifort            # ASE backend
```

- `<program_path>` は `akaikkr/specx`、`akaikkr_cnd/specx` などを含むディレクトリです。
- `--set` は `testsets.py` の `make_exe()` で定義した集合名（`all`, `Cu`, `Fe`, `Co`, `Ni`, `NiFe`, `AlMnFeCo`, `Fe2RhPt`, `Fe_lmd`, ...）。
- `--create_ref` は結果を `reference/<compiler>[_<backend>].json` に書きます。通常実行は同じファイルと比較して `SHORT SUMMARY` を出します（O 合格、X 不合格、- 参照無し）。
- 各物質のディレクトリ（`Cu/`, `Fe/`, ...）に inputcard、出力、図（dos.png、Awk_up.png ...）が残ります。`result.json` に今回の値、`Awk_both.png` がトップに出るのは既知の不便です。 dos.png / pdos_*.png の赤い一点鎖線は go の ewidth（同じディレクトリの `out_go.log` から読む）で、SCF の積分路の下端です（[dos_plot_ewidth_line.md](dos_plot_ewidth_line.md)）。

## 2. backend

| backend | 構造の作り方 | 参照ファイル |
|---|---|---|
| cif（既定） | pymatgen `Cif2Kkr`。慣用胞の `a` と `brvtyp=fcc` など | `reference/ifort.json` |
| ase | `ase.io.read` → `to_primitive` → `brvtyp=aux` | `reference/ifort_ase.json` |

両経路の全エネルギーは 1e-11 程度で一致しますが、Cu の aux と fcc で 1.8e-11 と閾値 1e-11 を超えるため、参照ファイルは分けてあります。type 名は経路で異なります（pymatgen は `Fe0.1Ni0.9_4a_0` のような Wyckoff 付き、ASE は出現順の番号）。type 順に依存する設定（SmCo5 の `mxl`）は type 名で決めています。

## 3. 環境（2026-09-24 時点）

- pandas 3 / scipy 1.17 / pymatgen 2026 で動くよう修正済み（`iteritems` → `items`、DataFrame の dtype=object など）。
- `akaikkr_testscript` は `pip install -e`。`pyakaikkr` も編集可能インストールに切り替えました。
- ifort ビルドはランタイムライブラリを oneAPI の `setvars.sh` で通しておく必要があります。

## 4. 既知の最終桁差（128 スレッド EPYC 7702 の参照と 24 スレッドのこのマシン）

`tests/testrun/akaikkr_cnd` の全 61 件のうち 6 件が X になります。いずれも閾値のすぐ外です。

| 計算 | 差 |
|---|---|
| Fe_lmd go / dos / spc31 | te rdiff 3e-10、moment 4e-5 |
| SmCo5_oc go | moment 差 1e-5 ちょうど |
| AlMnFeCo_bcc cnd | cnd rdiff 2.6e-5 |
| AlMnFeCo_bcc spc31 | dn の A(w,k) diff_max 5e-3（閾値 3e-3） |

閾値は変えていません。

## 5. 2026-09 に直した不具合

- `OutputAnalyzer._go_fix_localmoment` が参照行を結果側から作っていて、go / fsm / tc / cnd の比較が常に合格になっていた。
- `*_j30` の後処理で type 名を `"Rh0.5Pt0.5_1d_1"` と決め打ちしていた（pymatgen 2026 は `Pt0.5Rh0.5_1d_1`）。`typeofsite` から取るようにした。
- `OutputAnalyzer._canonical_type_name` で Jij の `pair` キーの元素順を正規化し、古い pymatgen で作った参照と突き合わせられるようにした（以前は不一致ペアが黙って落ちていた）。
- `GoFmg` の flip 名が一致しなくても無視されていた（`AlMnFeCo_bcc_gofmg` は Mn を反転していなかった）。未知の名前は例外、`*suffix` の末尾一致に対応。
- `tests/testrun/akaikkr_cpa2021v01/reference/ifort.json` は `tests/testrun/akaikkr` の参照のコピーで、cpa2021v01 の te（例 Cu_go −3304.747251823）とは合いません。cpa2021v01 用の参照は作り直しが必要です。

## 6. AiiDA から同じテストセットを流す

`../aiida-akaikkr/example/run_examples.py` は同じ `_<物質>_common_param` を使って全物質を AiiDA の CalcJob として投入し、`compare_reference.py` でここの `reference/ifort.json` と比べます。詳細は `aiida-akaikkr/docs/examples.md`。
