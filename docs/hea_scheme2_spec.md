# HEA 単一サイト CPA スクリーニング（run_scheme2）の pyakaikkr 移植仕様

作成日: 2026-09-25
移植元: `fukushima_HEA_run_exprlattice/production_run/run0/`（2019-11、`1306.run_scheme2.py`, `HEARun.py`, `HEAPathSearch.py`, `hea_util.py`, `kkrinput_brvtyp.py`, `akaikkrio2.py`）
移植先: AkaiKKRPythonUtil（pyakaikkr 2023.2.1）、AkaiKKR 2022.0721（`akaikkr/specx`）

## 0. 目的と範囲

- 4 元等比 HEA（各成分 25 at.%）を bcc / fcc の単一サイト CPA で計算し、go → dos → j を流して全エネルギー・モーメント・DOS・J_ij・Tc を得るスクリーニングを、pyakaikkr の `AkaikkrJob` の上に再実装する。
- 移植元のスキーム（ewidth の自動選定、edelt と pmix を段階的に落とす SCF 追い込み）はそのまま残す。ただし §1 の誤りは修正し、移植元の結果を再現したいときのために「互換モード」を用意する。
- 独自パーサ（`akaikkrio2.py`）と独自入力生成（`kkrinput_brvtyp.py`）は捨て、`AkaikkrJob.make_inputcard / run / get_*` に一本化する。
- 対象は akaikkr/specx のみ。akaikkr_cnd（displc 必須）と cpa2021v01 は対象外。
- 4 元・等比・単一サイトという移植元の制約は外し、任意の成分数と濃度（合計 100）、任意の brvtyp を受け付ける。

## 1. 移植元の誤りと修正方針

| # | 箇所 | 内容 | 修正 |
|---|---|---|---|
| 1 | `1306.run_scheme2.py` STEP1 | 新しい ewidth（`nextplan=="new"`）が返っても `iew` を進めないため、同じディレクトリ `ew_000` の既存 out_go.log を再利用して再計算しない。str.out で `('new', ...)` は 2053 回出ているが、そのすべてが ewidth=1.2 の結果を「新 ewidth の結果」として再判定している（2 回目は必ず `old` になる構造）。 | 実行ディレクトリ名に実際の ewidth 値を含め、inputcard の ewidth と一致しなければ再計算する（§4）。互換モードでも既存結果の再利用は inputcard の内容が一致する場合に限る。 |
| 2 | `HEARun.analyzeDOS.__init__` | 低 DOS 領域が上端まで続くとき、未定義名 `smalldos` を参照して NameError。 | `len(mask)` を使う（§5）。 |
| 3 | `akaikkrio2.OutputGo.converging` | `falg = True` のタイポ。当該分岐（`r2moment>0.8 and r2err>0.8`）は直前の `r2err>0.8` に含まれるので削除する。 | §6 の判定式に置き換える。 |
| 4 | `akaikkrio2.OutputDOS.get` | スピン和ではなく平均 `(up+dn)/2` を返しているので、`dosth=1e-2` は平均 DOS への閾値になっている。 | スピン和 DOS（states/Ry/cell）を使い、既定閾値は `dosth=2e-2` にする。互換モードでは同じ判定になる。 |
| 5 | `HEArun.make_heainput` | `self.dic` を破壊的に更新するため、go → dos → j と呼ぶうちに `go`, `ewidth`, `record` が前の呼び出しの値を引きずる。 | 入力辞書は毎回 `deepcopy` して作る。 |
| 6 | `HEArun.run_all` | go が未収束でも dos と j を実行する（`stopifnotconverged=False` 固定）。 | 既定では収束した go の後だけ dos / j を実行する。互換モードでは移植元どおり常に実行する。 |
| 7 | `1306.run_scheme2.py` | `nextplan=="fail"`（新 ewidth を作れない）を `keyconverged=True` にして「収束」と同じ扱いで終了する。 | 状態を `ewidth_fail` として記録し、収束とは区別する（§7）。 |
| 8 | `calculate_pm` | `ipm` をローカルで進めるので、呼び出し元の `ipm` と実際に使ったディレクトリ番号がずれる。 | 使った pm 段の番号を戻り値で返す。 |
| 9 | `analyzeDOS`, `make_new_ewidth` | 裸の `raise`（RuntimeError）。 | 専用例外 `HEASchemeError` を送出する。 |
| 10 | `HEADos` | 未定義の `OutFirstDOS` を使う。 | 削除し、プロットは `pyakaikkr.DosPlotter` に任せる。 |
| 11 | `hea_util` | `from pymatgen import Element` は pymatgen 2022 以降で動かない。 | `pymatgen.core.periodic_table.Element` を使う。 |
| 12 | `OutputGo.__init__` | 履歴が空だと `h_err[-1]` で IndexError。 | `AkaikkrJob.get_convergence`（`*** no convergence` の有無）で判定する。specx は `tol=1e-6`（log10 err ≤ -6）で収束、`maxitr` 到達で `*** no convergence` を出すので、移植元の `err <= -6` 判定と等価。 |

移植元で仕様どおりだったもの（誤りではない）:

- inputcard の `a=1000000` は標準 specx の機能で、`chklat.f` が `a > 0.99e6` のとき `qvolum(anclr,1)`（純物質の実験原子体積、cc/mol）を濃度平均して格子定数を決める。`a ≈ 0` なら MJW 値の表（`qvolum(anclr,2)`）を使う。合金化による体積変化は入らない。
- `record=2nd` は pm_000 でも指定されるが、pot.dat が無いときは specx が `eof detected; data generated` と警告して初期化から始めるので害はない。

## 2. ファイル構成

| 区分 | パス | 内容 |
|---|---|---|
| 追加 | `library/PyAkaiKKR/src/pyakaikkr/hea/__init__.py` | 公開 API（`HeaKey`, `make_hea_param`, `HeaRunner`, `Scheme2`, `choose_ewidth`, `is_converging`） |
| 追加 | `pyakaikkr/hea/keys.py` | heakey ⇄ 元素列（§3.1） |
| 追加 | `pyakaikkr/hea/inputs.py` | `make_hea_param`（§3.2） |
| 追加 | `pyakaikkr/hea/dosgap.py` | 低 DOS 領域の検出と ewidth 選定（§5） |
| 追加 | `pyakaikkr/hea/convergence.py` | SCF 履歴の「まだ落ちている」判定（§6） |
| 追加 | `pyakaikkr/hea/runner.py` | 1 ディレクトリで go / dos / j を実行する `HeaRunner`（§4） |
| 追加 | `pyakaikkr/hea/layout.py` | ディレクトリ命名と検索（§4.1） |
| 追加 | `pyakaikkr/hea/scheme2.py` | ewidth / edelt / pmix の状態機械 `Scheme2` と結果集計（§7, §8） |
| 追加 | `pyakaikkr/hea/cli.py` | `hea-scheme2` コマンド（§9） |
| 修正 | `pyakaikkr/Error.py` | `HEASchemeError(Exception)` を追加 |
| 修正 | `library/PyAkaiKKR/setup.cfg` | `console_scripts: hea-scheme2 = pyakaikkr.hea.cli:main`、`extras_require: hea = scikit-learn`（§6 で使う場合。numpy だけで書けるので既定は不要、後述） |
| 追加 | `tests/hea/` | pytest（§10） |
| 追加 | `docs/hea_scheme2_usage.md` | 使い方（実装後） |

`pyakaikkr.hea` は `import pyakaikkr` からは読み込まない（明示 import）。

## 3. 入力

### 3.1 heakey（`keys.py`）

```python
@dataclass(frozen=True)
class HeaKey:
    z: tuple[int, ...]              # 原子番号、入力順を保つ
    @classmethod
    def from_key(cls, key: str | int) -> "HeaKey"   # "13142122" → (13,14,21,22)
    @classmethod
    def from_elements(cls, elements: Sequence[str] | str) -> "HeaKey"  # ["Al","Si","Sc","Ti"] または "AlSiScTi"
    @property
    def key(self) -> str            # 2 桁ゼロ埋め連結。桁数は len(z)*2
    @property
    def elements(self) -> tuple[str, ...]
```

- 移植元の `heakey2zlist / heakey2elements / material2heakey` に対応する。桁数が奇数なら `ValueError`。
- `heakeylist0.csv`（列 `heakey, elements`）は `load_keylist(path) -> list[HeaKey]` で読む。`heakey` 列は文字列として読む（先頭ゼロの保護。pandas の `dtype={"heakey": str}`）。

### 3.2 inputcard 用辞書（`inputs.py`）

```python
def make_hea_param(key: HeaKey, brvtyp: str, *, conc: Sequence[float] | None = None,
                   lattice: str | float = "expr", go: str = "go", ewidth: float = 1.2,
                   edelt: float = 1e-4, pmix: float = 0.005, maxitr: int = 500,
                   bzqlty: int = 10, sdftyp: str = "pbe", reltyp: str = "sra",
                   magtyp: str = "mag", record: str = "2nd", outtyp: str = "update",
                   rmt: float = 1.0, mxl: int = 2, field: float = 0.0,
                   type_name: str = "HEA", potentialfile: str = "pot.dat",
                   option: dict | None = None) -> dict
```

- 戻り値は `AkaikkrJob.default` を `deepcopy` して更新した辞書で、`AkaikkrJob.make_inputcard` にそのまま渡せる。`ntyp=1, type=[type_name], ncmp=[len(z)], anclr=[list(z)], conc=[conc]`, `natm=1, atmicx=[["0.0a","0.0b","0.0c",type_name]]`, `c/a=b/a=1.0`, `alpha=beta=gamma=90`。
- `conc` 省略時は等比（`100/len(z)`）。`conc` の合計が 100 から 1e-6 以上ずれれば `ValueError`。等比のとき移植元と同じく整数 `25` を書く（4 元の場合）。5 元では `20`、3 元では `33.333333` のように小数になる。
- `lattice`:
  - `"expr"`（既定）: `a=1000000`。specx が実験原子体積の濃度平均から格子定数を決める（移植元の exprlattice）。
  - `"mjw"`: `a=0`。MJW 値の表を使う。
  - `float`: その値（bohr）を書く。
- `option` を渡すと `begin_option` ブロックを書く（`akaikkr_option_keys.md` 参照）。移植元には無い。
- 移植元との差: `AkaikkrJob.make_inputcard` は `#--- type ncmp / #- rmt field mxl / #- anclr conc` の短い形式を出すので、行の並びは `kkrinput_brvtyp.py` と同じになる。`a` は `1000000` を整数で書く（移植元と同じ文字列。specx は `f` 形式で読むので `1e6` でも動くが、互換のため整数）。

### 3.3 ewidth と DOS 計算

- go と j は同じ `ewidth`（スキームが決めた値）を使う。dos は `ewidth_dos=3.0`（固定、`Scheme2` の引数で変更可）で、specx は DOS を `[EF-0.75*ewidth_dos, EF+0.25*ewidth_dos]` に 201 点で出す（例: 3.0 のとき -2.2425 〜 0.7425 Ry、刻み 0.015 Ry）。§5 の判定はこのエネルギー軸上で行う。

## 4. 1 ディレクトリの実行（`runner.py`）

```python
class HeaRunner:
    def __init__(self, akaikkr_exe: str, directory: str, param_go: dict, *,
                 ewidth_dos: float = 3.0, files: dict | None = None, compat: bool = False)
    def run_go(self, copy_potential_from: str | None = None, force: bool = False) -> bool  # converged
    def run_dos(self, force: bool = False) -> None
    def run_j(self, force: bool = False) -> None
    def run_all(self, copy_potential_from=None) -> bool
    def result(self) -> dict     # §8 の 1 行分
```

- `files` の既定は `{"inputcard_go": "inputcard_go", "out_go": "out_go.log", "inputcard_dos": "inputcard_dos", "out_dos": "out_dos.log", "inputcard_j": "inputcard_j", "out_j": "out_j.log", "potential": "pot.dat"}`（移植元と同名）。
- `copy_potential_from` が与えられれば、実行前に `AkaikkrJob.copy_potential` で pot.dat をコピーする（移植元の `copyfiles_from`）。
- 再利用の判定（移植元の誤り #1 の修正）: 既存の inputcard を読んで、これから書く内容と**文字列として一致**し、かつ出力が正常終了（`sbtime report` あり、`***err` 無し）していれば実行しない。inputcard が異なるときは、既存の出力を `out_go.log.bak-<n>` に退避してから再実行する。`force=True` なら無条件に実行する。
- 実行は `AkaikkrJob.run`。`KKRFailedExecutionError`（return code ≠ 0、または `***err`）はそのまま上げる。
- 収束は `AkaikkrJob.get_convergence(out_go)`。
- `run_all` は go の後、`compat=False`（既定）なら収束したときだけ dos と j を実行し、`compat=True` なら常に実行する。
- dos と j の inputcard は `param_go` を `deepcopy` し `go` と（dos なら）`ewidth` だけ差し替えて作る。`record=2nd` で go の pot.dat を読む。
- `result()` は `AkaikkrJob` の `get_total_energy, get_total_moment, get_component_moment, get_lattice_constant, get_unitcell_volume, get_rms_error, get_te_history, get_moment_history, get_curie_temperature(out_j), get_jij_as_dataframe(out_j), get_dos(out_dos)` を使って §8 の項目を返す。
  - `get_lattice_constant` は `bravais=fcc a= 8.01512` の行から読むので、ヘッダの `a=*********` には影響されない。
  - `get_total_moment / get_component_moment` は `spin moment=` 行を読む。J は `get_jij_as_dataframe`。

### 4.1 ディレクトリ命名（`layout.py`）

```python
@dataclass
class RunPoint:
    key: str; iew: int; ewidth: float; ied: int; edelt: float; polytyp: str; ipm: int; pmix: float

class Layout:
    def __init__(self, prefix: str = "RUN", version: int = 2)
    def dirname(self, p: RunPoint) -> str
    def parse(self, dirname: str) -> RunPoint
    def find(self, **fixed) -> list[RunPoint]
```

- `version=1`（互換）: `key_{key},ew_{iew:03},ed_{ied:03},polytyp_{polytyp},pm_{ipm:03}`。移植元の `ordertype=2` と同じ。既存の RUN を読むときに使う。
- `version=2`（既定）: `key_{key},ew_{iew:03}-{ewidth:.4f},ed_{ied:03}-{edelt:.0e},polytyp_{polytyp},pm_{ipm:03}-{pmix:.0e}`。ewidth / edelt / pmix の実値を名前に含めるので、同じ `iew` で ewidth が変わっても別ディレクトリになる。
- `find` は glob ではなく `os.scandir` で列挙し、`parse` で戻す（名前にカンマを含むため `,` で区切り、各項目は最初の `_` で分ける）。
- 区切り文字のカンマは移植元の資産と互換にするため残す。新規運用で嫌うときは `Layout(sep=";")` を許す（既定はカンマ）。

## 5. 低 DOS 領域と ewidth の選定（`dosgap.py`）

```python
def low_dos_regions(energy: np.ndarray, dos_list: Sequence[np.ndarray], dosth: float = 2e-2
                    ) -> list[tuple[float, float]]
def choose_ewidth(regions, ewidth: float, *, eth: float = 0.30, ediff: float = 0.20,
                  margin: float = 0.01) -> tuple[str, float | None, list[float]]
```

- `energy` は E - EF（Ry）、`dos_list` は polytyp ごとのスピン和 total DOS（`get_dos` の `up+dn`）。全 polytyp のエネルギー軸が一致しなければ `HEASchemeError`。
- `low_dos_regions`: 全 polytyp で `dos < dosth` となる点の AND を取り、連続区間 `[e_start, e_end]` の列にする。区間の上端は「最後に条件を満たした点の次の点」（移植元と同じ、`i` を使う）。末尾まで続くときは最後の点（移植元の NameError 箇所）。
- `choose_ewidth` は移植元の `calc_new_ewidth` と同じ判定:
  1. 区間を高エネルギー側から見る。幅 `e2-e1 > eth` かつ `e1 < -ewidth < e2` かつ `-ewidth < e2 - ediff` を満たす区間があれば `("old", ewidth, candidates)`。
  2. 無ければ幅 `> eth` の各区間から `-(e2 - ediff - margin)` を候補にして `("new", candidates[0], candidates)`。候補が無ければ `("fail", None, [])`。
- 既定 `dosth=2e-2` は移植元の `1e-2`（スピン平均に対する閾値）と同じ判定。互換モードで `dosth_compat=1e-2` を平均 DOS に使っても同じ結果になるので、内部は常にスピン和で持つ。
- `eth`, `ediff`, `margin` は移植元のハードコード（0.30, 0.20, 0.01 Ry）を引数にしたもの。

## 6. 「まだ収束に向かっている」判定（`convergence.py`）

```python
def is_converging(err_history: Sequence[float], moment_history: Sequence[float], *,
                  last: int = 100, r2_th: float = 0.80, mae_err_th: float = 5e-4,
                  mae_moment_th: float = 5e-5, mae_err_loose_th: float = 6e-3) -> bool
```

- 直近 `last` 点の err（log10 rms）と moment に対して、反復番号を説明変数にした最小二乗直線を引き、標準化した値で r² を、生の値で MAE（残差の平均絶対値）を求める。
- 判定は移植元の生きている条件だけ: `r2_err > r2_th`、または `mae_err < mae_err_th`、または `mae_moment < mae_moment_th and mae_err < mae_err_loose_th`。タイポで死んでいた `r2moment and r2err` 条件は入れない。
- scikit-learn は使わず `numpy.polyfit` で書く（標準化は平均・標準偏差で自前。定数列のときは r²=0 とする）。
- 点数が 3 未満なら `False`。

## 7. スキーム（`scheme2.py`）

```python
class Scheme2:
    def __init__(self, akaikkr_exe: str, layout: Layout, *, polytyps=("bcc", "fcc"),
                 ewidth_init=1.2, ewidth_dos=3.0, edelt_steps=(1e-4, 1e-3, 1e-2),
                 pmix_steps=(1e-2, 5e-3, 1e-3, 5e-4, 1e-4), pmix_init=0.005,
                 maxitr_init=500, maxitr_2nd=200, maxitr_pm=300, max_pm_iter=20, max_ew=10,
                 dosth=2e-2, eth=0.30, ediff=0.20, lattice="expr", compat=False,
                 param_override: dict | None = None, logger=None)
    def run_key(self, key: HeaKey) -> KeyResult
    def run_keys(self, keys: Iterable[HeaKey]) -> list[KeyResult]
```

移植元の 2 段階を状態機械として書き直す。1 キーの処理:

1. **STEP1（ewidth の粗い決定）**: `iew=0, ewidth=ewidth_init, edelt=edelt_steps[0], pmix=pmix_init, maxitr=maxitr_init` で各 polytyp を 1 回だけ go → dos → j（`only_once`）。DOS から `choose_ewidth`。
   - `old` かつ全 polytyp 収束 → **finished**。
   - `old` だが未収束 → STEP2 へ。
   - `new` → `iew += 1`、ewidth を候補に置き換えて STEP1 をやり直す（誤り #1 の修正。`version=2` の layout なら別ディレクトリになる）。同じ ewidth を二度試さない（`ewidth_tried`）。`iew >= max_ew` で **ewidth_exhausted**。
   - `fail` → **ewidth_fail**（収束扱いにしない）。ただし go が全 polytyp 収束していれば結果は残す。
2. **STEP2（SCF の追い込み）**: `maxitr=maxitr_2nd`。`edelt_steps` を順に、各 edelt で `pmix_steps` を順に試す。各 pmix では前段の pot.dat をコピーして `maxitr=maxitr_pm` で最大 `max_pm_iter` 回まで継続し、収束したら次の polytyp、`is_converging` が偽なら次の pmix。全 pmix を使い切ったら次の edelt。
   - 各 edelt の後に DOS を見て `choose_ewidth`。`old` かつ全 polytyp 収束 → **finished**。`new` → 新 ewidth で STEP1 からやり直し。`fail` → **ewidth_fail**。
3. `iew` が `max_ew` に達する、または候補が尽きたら **not_converged**。

- `KeyResult` は `status`（`finished | ewidth_fail | ewidth_exhausted | not_converged | error`）、最終 `RunPoint`（polytyp ごと）、試した ewidth 列、§8 の結果行、例外メッセージを持ち、`RUN/key_<key>.json` に書く。1 キーの例外は捕まえて `error` にし、次のキーへ進む（移植元は `raise` で全体が止まる）。
- `compat=True` のときの差: `Layout(version=1)`、dos / j を常に実行、STEP1 の `new` で `iew` を進めずに既存結果を再利用（移植元の挙動の再現。データベース検証用で、新規計算には使わない）。
- 移植元の `pmix` の意味: STEP1 では `pmix_init=0.005`、STEP2 の pm 段では `pmix_steps` の値を `maxitr_pm=300` と組で使う。edelt は移植元どおり `param0["edelt"]` を上書きして inputcard に反映する。
- `param_override` は `make_hea_param` に渡す追加キー（`bzqlty`, `sdftyp`, `mxl`, `option` など）。

## 8. 出力（結果行）

`HeaRunner.result()` と `KeyResult` の 1 polytyp 分:

| 列 | 取得元 |
|---|---|
| key, elements, conc, polytyp, ewidth, edelt, pmix, iew, ied, ipm, directory | `RunPoint` |
| lattice_mode, a_bohr, volume_bohr3 | `get_lattice_constant`, `get_unitcell_volume`（out_go） |
| converged, n_iter, last_err | `get_convergence`, `len(get_rms_error)`, `get_rms_error()[-1]` |
| total_energy_Ry | `get_total_energy`（out_go） |
| total_moment, component_moment[] | `get_total_moment`, `get_component_moment`（out_go） |
| Tc_K | `get_curie_temperature`（out_j） |
| jij_csv | `get_jij_as_dataframe(out_j)` を `jij.csv` に保存したパス |
| dos_csv | `get_dos(out_dos, "dataframe")` を `dos.csv` に保存したパス |
| low_dos_regions | §5 の区間 |

集計は `Scheme2.collect(prefix) -> pd.DataFrame`（全 `key_*.json` を読む）。既存の 2019 年の RUN を読むときは `Layout(version=1)` と `HeaRunner.result()` だけを使う `collect_legacy(prefix)` を用意する（inputcard から ewidth 等を読み直す。ディレクトリ名の `iew` は信用しない）。

## 9. CLI

```
hea-scheme2 run   --exe /path/to/specx --keys heakeylist0.csv [--prefix RUN] [--polytyp bcc fcc]
                  [--lattice expr|mjw|<bohr>] [--ewidth-init 1.2] [--compat] [--threads 24] [--start N --stop M]
hea-scheme2 collect --prefix RUN [--legacy] -o result.csv
hea-scheme2 key   AlSiScTi          # heakey の相互変換
```

- `--threads` は `OMP_NUM_THREADS` を子プロセスにだけ設定する（移植元は `MKL_NUM_THREADS=1` を固定していたが specx は MKL を使わない）。
- ログは `logging` で `prefix/hea-scheme2.log` に書く。移植元の `_1306.run_scheme2.warning.log` の「Failed to have a new ewidth」に相当するものは `status=ewidth_fail` として JSON にも残す。

## 10. テスト（`tests/hea/`）

specx 不要:

- `test_keys.py`: `"13142122"` ⇄ `("Al","Si","Sc","Ti")`、奇数桁で `ValueError`、5 元。
- `test_inputs.py`: `make_hea_param` の inputcard が `run0/RUN/key_13487580,ew_000,ed_000,polytyp_fcc,pm_000/inputcard_go`（AlCdReHg fcc）と、`#` 行と空白の正規化を除いて一致すること。`lattice="mjw"` で `a=0`、float で bohr 値。
- `test_dosgap.py`: 上記ディレクトリの `out_dos.log` から `get_dos` で作った DOS に対し、`low_dos_regions` が str.out に記録されたギャップ領域（例: key 13142122 で `[-2.2425,-2.1075],[-2.0025,-0.7275]`）を返し、`choose_ewidth(1.2)` が `("old", 1.2, [0.9375])` を返すこと。合成 DOS で `new` と `fail`、末尾まで低 DOS が続くケース（旧 NameError）。
- `test_convergence.py`: 単調減少列で `True`、雑音のみで `False`、`last` 未満の長さ。
- `test_layout.py`: v1 / v2 の往復、旧ディレクトリ名の `parse`。
- `test_result_legacy.py`: 上記ディレクトリの out_go / out_dos / out_j を `HeaRunner.result()` で読み、`total_energy=-21075.683771993`、`a=8.01512`、`Tc=0.0`、全成分モーメント 0 を得ること。

specx が必要（`AKAIKKR_PROGRAM_PATH` があるときだけ）:

- `test_run.py`: Cu 単体（`HeaKey.from_elements("Cu")`, fcc, `lattice=6.82`）で `HeaRunner.run_all` が収束し、`tests/akaikkr/reference/ifort.json` の Cu_go の te と一致すること（ewidth と bzqlty を参照と揃える）。
- `test_scheme2_small.py`: 1 キー（AlSiScTi、bcc のみ、`max_pm_iter=2`, `maxitr_init=50`）で `Scheme2.run_key` が例外なく `status` を返すこと。

## 11. 実装時の注意

- `AkaikkrJob.make_inputcard` は `atmicx` を `["0.0a","0.0b","0.0c","HEA"]` の形で受ける。移植元の `0 0 0 HEA` と等価。
- `record=2nd` で pot.dat が無いときの `***wrn in spmain...eof detected` は正常。`check_stopped_by_errtrp` は最終行の `***err in` だけを見るので影響しない。
- specx は `a=1000000` をヘッダで `a=*********` と印字する。`get_struc_param` も `bravais=` の行から読むことを確認する。
- DOS の閾値判定は E - EF の軸で行う。`get_dos` が返すエネルギーは out_dos.log の 1 列目そのもの（EF 基準）。
- 移植元の `pot.dat.info`（specx が go / dos / j ごとに 1 行追記する a, te, moment）は読まない。同じ情報は out_*.log から取る。
- 2019 年の RUN を `collect_legacy` で読むときは、誤り #1 のため `ew_000` ディレクトリに「ewidth=1.2 で計算した結果」しか無い。ewidth は必ず inputcard_go から読み直す。
