# GAES（ewidth 自動調整）の使い方

対象: pyakaikkr の `pyakaikkr.gaes`（2026-09-25 実装）。設計は [ewidth_tuning_scheme.md](ewidth_tuning_scheme.md)。

GAES は T. Fukushima, H. Akai, T. Chikyow, H. Kino, *Phys. Rev. Materials* **6**, 023802 (2022), [doi:10.1103/PhysRevMaterials.6.023802](https://doi.org/10.1103/PhysRevMaterials.6.023802) の HEA 網羅計算で使われた手法だが、手法そのものは未発表（同論文に記述は無い）。手法の記述は [ewidth_tuning_scheme.md](ewidth_tuning_scheme.md) が最初のもの。

## 何をするか

go の ewidth（E_F − ewidth_go が SCF の積分路の下端）が、valence 帯と semicore / core の間の**バンドギャップ**（連続した energy mesh 区間で DOS(E) < threshold）に入るまで、ewidth を置き直して go → dos を繰り返す。同時に STEP2 で edelt を 1e-2 → 1e-3 → 1e-4 と締め、各 edelt で pmix を段階的に落として SCF を追い込む。判定は DOS の値の threshold 比較だけで、微分は使わない。

## 既存の dos 出力に判定だけ掛ける

```
kkr-gaes check --dos bcc/out_dos.log --dos fcc/out_dos.log --ewidth 1.2 [--dosth 1e-3] [--eth 0.3] [--ediff 0.2]
```

Python では

```python
from pyakaikkr import AkaikkrJob
from pyakaikkr.gaes import dos_curves_from_outputs, gap_regions, choose_ewidth
energy, curves, _ = dos_curves_from_outputs([(AkaikkrJob("bcc"), "out_dos.log"), (AkaikkrJob("fcc"), "out_dos.log")])
regions = gap_regions(energy, curves, dosth=1e-3)      # 全曲線の AND、[GapRegion(e1, e2, i1, i2), ...]
flag, ewidth, candidates = choose_ewidth(regions, 1.2)   # ("old", 1.2, [...]) | ("new", 0.94, [...]) | ("fail", None, [])
```

## スキームを回す

```python
from pyakaikkr.gaes import Gaes, Layout, SiteComposition, make_single_site_param
comp = SiteComposition.from_type_name("Rh0.5Pt0.5")          # 任意の元素と比。"AlSiScTi" は等比
params = {pt: make_single_site_param(comp, pt) for pt in ("bcc", "fcc")}   # 2019 年と同じ既定値、a=1000000
g = Gaes("/path/to/specx", Layout("RUN", version=2), ewidth_init=1.2, ewidth_dos=3.0, ref=0.75, dosth=1e-3)
res = g.run(comp.key(), params)      # KeyResult。RUN/key_<key>.json にも書く
print(res.status, res.ewidth_final, res.gap_used, [j.flag for j in res.judgements])
```

- `status`: `finished`（E_F − ewidth がギャップ内で全 polytyp 収束）、`ewidth_fail`（幅 eth 以上のギャップが無い）、`ewidth_exhausted`（`max_ew` 回試した）、`not_converged`、`error`。
- `params` の値は `AkaikkrJob.make_inputcard` に渡せる辞書なら何でもよい（CIF や ASE から作ったものも可）。`ewidth`, `edelt`, `pmix`, `maxitr` はスキームが上書きする。
- ディレクトリ名は `key_<key>,ew_<n>-<ewidth>,ed_<n>-<edelt>,polytyp_<p>,pm_<n>-<pmix>`。同じ inputcard で正常終了した出力があれば再実行しない。
- `ref` は使うビルドの dos の窓に合わせる（akaikkr 0.75、cpa2021v01 / cnd 0.5）。ずれていれば `RUN/kkr-gaes.log` に警告が出る。
- `ewidth_dos` は E_F − ewidth_go − eth − ediff まで届くよう自動で広げる（`ewidth_dos_auto=True`）。specx が dos で止まれば狭めて再試行する。

CLI:

```
kkr-gaes site --exe /path/to/specx --comp Rh0.5Pt0.5 AlSiScTi --polytyp bcc fcc [--ref 0.5] [--threads 20]
kkr-gaes hea  --exe /path/to/specx --keys heakeylist0.csv --start 0 --stop 10      # 2019 年のキーリスト
kkr-gaes run  --exe /path/to/specx --input params.json                              # {key: {polytyp: param_go}}
kkr-gaes collect --prefix RUN -o result.csv        # key_*.json を表に
kkr-gaes collect --prefix run0/RUN --legacy -o legacy.csv   # 2019 年の RUN を読む
kkr-gaes comp Rh0.5Pt0.5_1 13142122               # 組成・type 名・heakey の相互変換と 40 文字検査
kkr-gaes site --exe ... --comp RbMnFeCo --polytyp fcc --orbital Rb4p=core        # Rb 4p を積分路の外（core）に
kkr-gaes site --exe ... --comp BiMnFeCo --polytyp fcc --orbital Bi6s=occupied --orbital Bi5d=core
kkr-gaes check --dos RUN/.../out_dos.log --go RUN/.../out_go.log --ewidth 1.2 --orbital Rb4p=core   # 範囲と準位を表示
```

## 軌道の valence / core 指定（`--orbital`、`Gaes(orbitals=[...])`）

`Rb4p=valence`（別名 `occupied`: 積分路に入れる）、`Rb4p=core`（別名 `unoccupied`: 積分路から外す）のように、元素・軌道ごとに指定できる（仕様: ewidth_tuning_scheme.md §15）。

- 範囲: valence → `min_ewidth = |E − E_F| + ediff`、core → `max_ewidth = |E − E_F| − ediff`。準位は各 go の `out_go.log` の成分ブロック（`AkaikkrJob.get_core_levels_by_component`）から判定のたびに読む。core 扱いの準位は valence 扱いより 0.7〜0.9 Ry 深く出るので、その key で見た値のうち core 指定には最も浅い値、valence 指定には最も深い値を使う。
- 最初の go の前は `converged_core_levels_2019.csv`（2019 年の収束値）だけで範囲を決め、`ewidth_init` が外れていれば範囲の中へ動かす。表に無い元素は `ewidth_init` のまま始め、最初の判定で直す。
- 指定があると既定の [1.0, 2.0] は使わない。`--min-ewidth` / `--max-ewidth` を明示すれば共通部分を取る。矛盾（min > max）や core 配置に無い軌道の `core` 指定は `GaesError`（status `error`）。
- ギャップの判定は範囲に依らず窓全体で行い、範囲は ewidth の選択だけに使う。範囲を外れた ewidth は `old` にならない（2026-09-25 修正）。go の `*` が指定と食い違えば（準位が積分路をまたいだ）、範囲を決め直して次の候補へ。
- `KeyResult.judgements[i].orbital_levels / orbital_bounds / orbital_mismatch`、`parameters["orbital_bounds_step0"]` に記録される。

## 2019 年の RUN との互換

`Gaes(compat=True)` は 2019 年の挙動（ディレクトリ名 v1、dosth 2e-2 = スピン平均 1e-2、dos / j を常に実行、`fail` で即終了、新 ewidth でも同じディレクトリを再利用）を再現する。検証用で、新規計算には使わない。

## 注意（AkaiKKR 2022.0721）

- akaikkr ビルドは Hf を含む go が `reconf` で止まる。cpa2021v01 ビルド（`supprs=.true.`）を使う。詳細は設計書 §13.1。
- cpa2021v01 の dos は ref=0.5 固定なので、2019 年と同じ窓（下端 −2.25 Ry）には `ewidth_dos=4.5` が要る。

## テスト

```
cd tests/gaes && python -m pytest -q                                   # specx 不要 17 件（2019 年の RUN があれば含む）
AKAIKKR_PROGRAM_PATH=/path/to/AkaiKKRprogram.2022.0721.ifort python -m pytest -q test_gaes_run.py   # 3 系、cpa2021v01
```
