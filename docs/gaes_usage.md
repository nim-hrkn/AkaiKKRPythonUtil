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
- ギャップ区間の判定は既定で **原子あたり** の DOS（total DOS / natm）で行う（`dos_per_atom=True`、CLI `--dos-per-cell` で生の胞あたり）。閾値 2e-2 / 1e-3 は 1 原子胞で決めた値なので、多原子胞ではこの正規化が要る（Bi2Se3 の 5 原子胞では胞あたりだと fail）。図の DOS は胞あたりのままで、閾値の線だけ natm 倍される。

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

## DOS 図の共通部品（`pyakaikkr.gaes.plot`）

GAES の DOS 図はどれも同じ要素（緑 = 粗いギャップ区間、青 = 細かい部分区間、斜線 = [min_ewidth, max_ewidth] の帯、赤 = −ewidth_go、青線 = 注目する core 準位（破線 `*`）、灰線 = 他の準位）を描くので、`draw_gaes_dos` に集約した。

```python
from pyakaikkr.gaes import dos_curves_from_outputs, levels_from_go
from pyakaikkr.gaes.ewidth import choose_ewidth2
from pyakaikkr.gaes.plot import draw_gaes_dos, legend_text
e, c, _ = dos_curves_from_outputs([(AkaikkrJob(d), "out_dos.log")])
dec = choose_ewidth2(e, c, 1.2, dosth=2e-2, dosth2=1e-3, min_ewidth=1.0, max_ewidth=2.0)
draw_gaes_dos(ax, e, c[0], ewidth=1.2, decision=dec, bounds=(1.0, 2.0), levels=levels_from_go(d + "/out_go.log"),
              highlight=["Bi6s"], xlim=(-2.4, 0.8))
fig.suptitle(legend_text())
```

`tests/gaes/tools/plot_hea_XMnFeCo.py`, `plot_orbital_rules.py`, `plot_conv_pairs.py` がこれを使う（docs/data の図）。aiida-akaikkr の `plot --gaes-pk` は自前の配色で同じ要素を描く。

## 軌道の valence / core 指定（`--orbital`、`Gaes(orbitals=[...])`）

`Rb4p=valence`（別名 `occupied`: 積分路に入れる）、`Rb4p=core`（別名 `unoccupied`: 積分路から外す）のように、元素・軌道ごとに指定できる（仕様: ewidth_tuning_scheme.md §15）。

- 範囲: valence → `min_ewidth = |E − E_F| + ediff`、core → `max_ewidth = |E − E_F| − ediff`。準位は各 go の `out_go.log` の成分ブロック（`AkaikkrJob.get_core_levels_by_component`）から判定のたびに読む。core 扱いの準位は valence 扱いより 0.7〜0.9 Ry 深く出るので、その key で見た値のうち core 指定には最も浅い値、valence 指定には最も深い値を使う。
- 最初の go の前は `converged_core_levels_2019.csv`（2019 年の収束値）だけで範囲を決め、valence 指定の下限より `ewidth_init` が浅ければ下限 + ediff へ深くする。core 指定のために浅くはしない（最初の go の DOS から範囲内の候補を取る）。表に無い元素は `ewidth_init` のまま始め、最初の判定で直す。valence 指定の下限が深いときは dos の窓をそこまで広げる。
- 指定があると既定の [1.0, 2.0] は使わない。`--min-ewidth` / `--max-ewidth` を明示すれば共通部分を取る。矛盾（min > max）や core 配置に無い軌道の `core` 指定は `GaesError`（status `error`）。
- ギャップの判定は範囲に依らず窓全体で行い、範囲は ewidth の選択だけに使う。範囲を外れた ewidth は `old` にならない（2026-09-25 修正）。go の `*` が指定と食い違えば（準位が積分路をまたいだ）、範囲を決め直して次の候補へ。
- `KeyResult.judgements[i].orbital_levels / orbital_bounds / orbital_mismatch`、`parameters["orbital_bounds_step0"]` に記録される。

Python から:

```python
from pyakaikkr.gaes import Gaes, Layout, SiteComposition, make_single_site_param
comp = SiteComposition.from_elements(["Rb", "Mn", "Fe", "Co"])
params = {"fcc": make_single_site_param(comp, "fcc", type_name="HEA")}
g = Gaes("/path/to/specx", Layout("RUN_orb", version=2), ewidth_init=1.2, ref=0.75,
         orbitals=["Rb4p=core"],          # 複数可: ["Bi6s=occupied", "Bi5d=core"]。occupied = valence、unoccupied = core
         with_j=False)                    # min_ewidth / max_ewidth は指定すれば共通部分、指定しなければ規則が決める
res = g.run("RbMnFeCo", params)
print(res.status, res.ewidth_final)      # finished / ewidth_fail / error（規則の矛盾、core 配置に無い軌道の core 指定）
print(res.message)                       # ewidth_fail のとき: 範囲と、区間ごとに候補が出なかった理由
for j in res.judgements:
    print(j.ewidth, j.flag, j.orbital_bounds,            # この判定に使った [min_ewidth, max_ewidth]
          j.orbital_levels["Rb4p"],                      # [E - E_F, star]  star=True は valence（*）
          j.orbital_mismatch,                            # go の扱いが指定と違った軌道
          j.reasons)                                     # fail のとき、区間ごとの理由
print(res.parameters["orbitals"], res.parameters.get("orbital_bounds_step0"))
```

`kkr-gaes check` は既存の go / dos 出力に対して同じ判定だけを行い、`note:` 行に理由を出す:

```
$ kkr-gaes check --dos RUN_orb/.../out_dos.log --go RUN_orb/.../out_go.log --ewidth 0.641 --orbital In4d=core
core levels (E - E_F, star): {'In4d': [-1.3473, False], ...}
  In4d=core: level -1.3473 -> max 1.1473
bounds from orbital rules: [None, 1.1473]
...
ewidth 0.6410: fail
candidates:
  note: sub-region [-1.103, -0.938] excluded by ediff 0.20 (E_F - ewidth must be below -1.123 = coarse upper edge -0.922 - ediff); ediff < 0.170 would give a candidate near 1.0925
  note: sub-region [-2.243, -1.657]: candidate 1.6675 is outside [min_ewidth, max_ewidth] = [None, 1.1473] and the bound is not inside the sub-region
```

図: `python tests/gaes/tools/plot_orbital_rules.py <tag> ...`（RUN_orb*/<tag>/key_*.json を読み、最終判定の DOS に規則の準位・範囲・理由を重ねる。`OUT=file.png` で出力名）。例は docs/data の `hea_XMnFeCo_fcc_orbital_rules_dos.png`（17 例）と `hea_XMnFeCo_fcc_orbital_pairs_dos.png`（Rb 4p / Se 4s / Bi 6s の valence と core の対）。

## 2019 年の RUN との互換

`Gaes(compat=True)` は 2019 年の挙動（ディレクトリ名 v1、dosth 2e-2 = スピン平均 1e-2、dos / j を常に実行、`fail` で即終了、新 ewidth でも同じディレクトリを再利用）を再現する。検証用で、新規計算には使わない。

## 注意（AkaiKKR 2022.0721）

- akaikkr ビルドは Hf を含む go が `reconf` で止まる。cpa2021v01 ビルド（`supprs=.true.`）を使う。詳細は設計書 §13.1。
- cpa2021v01 の dos は ref=0.5 固定なので、2019 年と同じ窓（下端 −2.25 Ry）には `ewidth_dos=4.5` が要る。
- ewidth を微少動かしただけで SCF が収束する / しないことがある（下端が準位や帯の端をまたぐため）。収束の可否は ewidth がギャップにある証拠にも反証にもならないので、判定は DOS で行い、収束しない go の DOS も判定に使う（設計書 §12）。

## テスト

```
cd tests/gaes && python -m pytest -q                                   # specx 不要 17 件（2019 年の RUN があれば含む）
AKAIKKR_PROGRAM_PATH=/path/to/AkaiKKRprogram.2022.0721.ifort python -m pytest -q test_gaes_run.py   # 3 系、cpa2021v01
```
