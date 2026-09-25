# E_F、原子（成分）の core 準位と DOS の semicore ピークの対応

作成日: 2026-09-25
データ: `fukushima_HEA_run_exprlattice/production_run/run0/RUN/`（2019 年、4 元等比 HEA の単一サイト CPA、fcc の 7,505 ディレクトリ）と、2026-09-25 に AkaiKKR 2022.0721 で走らせた GAES のテスト（[ewidth_tuning_scheme.md](ewidth_tuning_scheme.md) §13）
スクリプト: `tests/gaes/tools/scan_semicore.py`（DOS からギャップと semicore 区間を抽出）、`tests/gaes/tools/corelevel_vs_dos.py`（core 準位と DOS ピークの突き合わせ）

## 1. out_go.log のどこに何が書いてあるか

### E_F

SCF ループが終わった直後（最後の反復の `rms err=` 行の次）に 1 行だけ出る。反復途中には出ない。

```
   itr=212  neu  0.0000  chr,spn  7.7500  1.5151  intc,ints  0.7663 -0.0121
   rms err= -6.254 -6.055 -6.237 -6.426 -6.257 -6.053
            -6.238 -6.427
   ef=   0.6006938   0.6049192  def=   2.6920550   9.0556286
   total energy=  -4277.933100646
```

- `ef=` の 2 つの値はスピン up / down の Fermi 準位（Ry、絶対値、muffin-tin ゼロ基準）。`magtyp=mag` では僅かに違い、`nmag` では同じ値が 2 回出る。
- `def=` は E_F での DOS（states/Ry、スピンごと）。
- `out_dos.log`、`out_j.log` にも同じ形式で 1 行出る。
- pyakaikkr: `AkaikkrJob.get_fermi_level(outfile)` が up / down の平均を返す。

### 成分ごとの core 準位

出力末尾の、成分ごとのブロック `*** type-<type>  <元素> (z= <Z>) ***` の中に `core level (spin up)` / `(spin down)` として出る（`get_component_moment` が読む `spin moment=` 行と同じブロック）。

```
                             *** type-HEA      In (z= 49.0) ***
   ...
   core level  (spin up  )
     -2037.8783 Ry(1s)      -304.3372 Ry(2s)      -273.0079 Ry(2p)
       -57.0653 Ry(3s)       -46.8159 Ry(3p)       -30.9243 Ry(3d)
        -7.9472 Ry(4s)        -4.8964 Ry(4p)        -0.5365 Ry(4d)*
```

- 値は `ef=` と同じ絶対エネルギー基準。**E − E_F = 準位 − ef** で DOS の横軸に直せる。
- どの状態が core かは、その前の `core configuration for Z= 49` 表（各軌道の core 電子数）で決まる。表に無い軌道（例: Cd の 4d）は最初から valence で、core level には出ない。
- 準位の後ろの `*` は「core 配置ではあるが、積分路の下端 ebtm = E_F − ewidth より上にあるため占有を valence 側に移した」印（`cstate.f` の 186〜191 行、`reconf` の結果）。`*` の準位は DOS にピークとして現れる。
- pyakaikkr: `AkaikkrJob.get_core_level(outfile, core_state=[...])` は指定軌道の値を返すが、成分の区別と `*` は返さない。GAES で使うときはブロック単位で読む（§5）。

## 2. 元素ごとの semicore の位置（2019 年の RUN、E_F 基準）

`corelevel_vs_dos.py` で、fcc の 7,505 ディレクトリの各成分について −2.3 Ry より浅い core 準位を取り、total DOS（スピン和、edelt 1e-4、窓 [−2.24, 0.74] Ry）の極大（高さ > 1、E < −0.85）のうち最も近いものと比べた。

| 元素 | 準位 | 扱い | 件数 | 準位（中央値, Ry） | DOS ピーク − 準位（中央値, Ry） | 0.1 Ry 以上ずれた割合 | ピーク高さ（中央値, states/Ry） |
|---|---|---|---|---|---|---|---|
| Ge | 3d | core | 1,270 | −1.81 | −0.001 | 0% | 164 |
| Sn | 4d | core | 1,110 | −1.62 | 0.000 | 0% | 105 |
| Bi | 5d | core | 1,102 | −1.73 | 0.001 | 0% | 103 |
| Pb | 5d | core | 1,090 | −1.29 | 0.000 | 0% | 62 |
| Sc | 3p | core | 1,282 | −2.06 | 0.002 | 0% | 52 |
| Y | 4p | core | 1,276 | −1.68 | 0.013 | 0% | 26 |
| Zr | 4p | core | 1,272 | −1.98 | 0.011 | 0% | 29 |
| In | 4d | valence（`*`） | 1,118 | −1.10 | −0.001 | 0% | 69 |
| Ga | 3d | valence（`*`） | 1,276 | −1.10 | −0.002 | 0% | 125 |
| Tl | 5d | valence（`*`） | 1,124 | −0.90 | −0.001 | 5% | 38 |
| Bi | 6s | valence（`*`） | 1,102 | −0.88 | −0.84 | 76% | （独立ピーク無し） |
| Hf | 4f | core | 744 | −1.84 | 0.045 | 90% | 42 |
| Hf | 5p | core | 744 | −2.25 | 0.45 | 64% | （窓の下端） |
| Nb | 4p | core | 622 | −2.28 | 0.63 | 94% | （窓の下端） |

読み方:

- **core 扱いの semicore（Ge 3d、Sn 4d、Bi 5d、Pb 5d、Sc 3p、Y 4p、Zr 4p）は、DOS のピークが core 準位と数 mRy 以内で一致する**。例外は 1 件も無い。d 状態（Ge、Sn、Bi、Pb）は鋭く高く（60〜165）、p 状態（Sc、Y、Zr）は幅があり低い（25〜50）。
- `*` 付きで valence に切り替えられた準位も、In 4d、Ga 3d、Tl 5d のように鋭いまま残るものは一致する。Bi 6s は valence 帯に溶けて独立ピークにならない。
- Hf 4f は 2019 年版では core（−1.84 Ry）でピークが 0.05 Ry ずれる。2022.0721 版では 4f が −1.2 Ry に来て積分路の下端と衝突する（[ewidth_tuning_scheme.md](ewidth_tuning_scheme.md) §13.1）。
- Hf 5p、Nb 4p は dos の窓の下端（−2.24 Ry）に掛かり、ピークが切れているので一致しない。窓の外の準位（Pt / Au / Hg の 5p が −3.9〜−4.8 Ry、4f が −6.2〜−8.3 Ry、In 4p −5.5 Ry など）は DOS には現れないが core level 表には出る。

### 2.1 同族での比較（4d 系列）

| 元素 | 4d の扱い | 4d の位置 E − E_F | DOS での見え方 |
|---|---|---|---|
| Cd | valence（core 配置の 4d 列が 0） | −0.73 Ry | valence 帯の底の一部（高さ 27）。ギャップ無し |
| In | core 配置だが ebtm より上（`*`） | −1.10 Ry（範囲 −1.17〜−1.04） | valence 帯から 0.3 Ry 離れた鋭いピーク（高さ 50〜70）。2019 年の ewidth 1.2 はその下側の裾を切る |
| Sn | core | −1.62 Ry | 上にギャップ [−1.24, −0.98]（幅 0.26 Ry） |

同様に Tl 5d（−0.90、`*`）→ Pb 5d（−1.29、core）→ Bi 5d（−1.73、core）、Ga 3d（−1.10、`*`）→ Ge 3d（−1.81、core）と、周期表で右へ行くほど深くなる。

### 2.2 semicore ピークを作らない元素

Pt、Au、Hg は 5d が valence 帯の中（Hg 5d は −0.7 Ry 付近の山）で、最も浅い core 準位は 5p（−3.9〜−4.8 Ry）と 4f（−6.2〜−8.3 Ry）。RUN の Pt / Au / Hg 系で「ギャップの下の semicore」と判定された例は多数あるが、その semicore はすべて同居する Ge / Bi / Sn / Y / Zr 側のものだった。Cd も同じ。Al と 3d 遷移金属だけの組では、valence 帯の底（−0.8 Ry）から窓の下端まで DOS は単調な裾で、ギャップは検出されない。

## 3. 4 元 HEA の DOS の semicore ピークを元素に割り当てる

上の一致から、**total DOS のピークは out_go.log の core level 表を引くだけで元素と軌道に割り当てられる**。PDOS は要らない。手順:

1. `ef=` 行から E_F（up / down の平均、または片方）を取る。
2. 成分ブロックごとに core level を読み、E − E_F に直す。`*` の有無も記録する。
3. DOS の窓の中にある準位（E − E_F > 窓の下端）ごとに、最も近い DOS の極大を探す。差が 0.03 Ry 以内なら同一。
4. 窓の下端から 0.1 Ry 以内の準位はピークが切れているので割り当てない。

`corelevel_vs_dos.py` がこれを行い、`corelevel_vs_dos.json` に (成分, 軌道, `*`, E − E_F, ピーク位置, 高さ) を書く。

## 4. GAES への含意

- ギャップ判定は DOS の値だけで行う（[ewidth_tuning_scheme.md](ewidth_tuning_scheme.md) §0.2）が、**E_F − ewidth_go が core 準位に衝突していないか**は go の出力だけで、dos を回す前に検査できる。Hf 4f（2022 版で −1.2 Ry）や Pb 5d（−1.29 Ry）のように、2019 年の初期値 1.2 がちょうど準位の上に乗る場合は SCF が収束しない、または収束しても積分路が semicore を切っている。
- 実装案（未実装）: 成分ごとの core 準位（E − E_F）を幅 0.1 Ry 程度の擬似ピークとして `gap_regions` の曲線列に加える。窓の外の準位や、`*` で valence に切り替わる準位の位置も判定に入る。
- `*` 付きの準位は ewidth を変えると 0.05〜0.1 Ry 動く（In 4d: ewidth 1.2 → 1.46 → 1.53 で −1.10 → −1.17 → −1.13 Ry）。判定のたびに読み直す。

## 5. pyakaikkr での読み方（現状）

```python
import re
from pyakaikkr import AkaikkrJob
job = AkaikkrJob("RUN/key_13252649,ew_000,ed_000,polytyp_fcc,pm_000")
ef = job.get_fermi_level("out_go.log")                  # up/down の平均
txt = "\n".join(job._read("out_go.log"))
for block in re.split(r"\*\*\* type-", txt)[1:]:
    m = re.match(r"(\S+)\s+(\S+)\s+\(z=\s*([\d.]+)\)", block)   # type, element, Z
    levels = re.findall(r"(-?\d+\.\d+) Ry\((\d[spdf])\)(\*?)", block)
    for value, orbital, star in levels:
        print(m.group(2), orbital, float(value) - ef, "valence" if star else "core")
```

成分ごとの core 準位を `*` 付きで返す `AkaikkrJob.get_core_levels_by_component(outfile)` は未実装（GAES の core 準位判定を入れるときに追加する）。
