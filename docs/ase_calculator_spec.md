# AkaiKKR ASE calculator ソース修正仕様（混晶 / CPA 対応を含む）

作成日: 2026-09-24
対象: AkaiKKRPythonUtil（pyakaikkr 2023.2.1）、ase 3.29.0、AkaiKKR 2022.0721（akaikkr / akaikkr_cpa2021v01）

## 0. 目的と範囲

- ASE の `Atoms` から AkaiKKR の inputcard を生成し、`specx` を実行し、結果を ASE の `Calculator` の流儀で返す。
- 混晶（CPA）に対応する。ASE 標準の部分占有表現（`atoms.info['occupancy']` と `atoms.arrays['spacegroup_kinds']`）を AkaiKKR の `type / ncmp / anclr / conc` に写像する。
- energy-only calculator とする。`forces` と `stress` は実装しない（AkaiKKR は力を出力しないため）。
- 既存の `AkaikkrJob`（`make_inputcard`, `run`, `get_*`）と `ElementKKR` を再利用し、pymatgen 経由の `Cif2Kkr` には依存しない。

以下、「追加」は新規ファイル、「修正」は既存ファイルの変更を指す。

## 1. ファイル構成

| 区分 | パス | 内容 |
|---|---|---|
| 追加 | `src/pyakaikkr/ase/__init__.py` | `AkaiKKR` calculator と補助関数を公開する |
| 追加 | `src/pyakaikkr/ase/structure.py` | `Atoms` → 格子・原子座標（`brvtyp=aux`, `a`, `r1..r3`, `natm`, `atmicx`） |
| 追加 | `src/pyakaikkr/ase/occupancy.py` | 部分占有の検査・正規化、kind と対称性による type 決定、`ntyp / type / ncmp / anclr / conc / rmt / field / mxl` の生成 |
| 追加 | `src/pyakaikkr/ase/calculator.py` | `class AkaiKKR(FileIOCalculator)` |
| 修正 | `src/pyakaikkr/AkaiKkr.py` | デバッグ print の除去、成分別モーメント取得の追加（§5） |
| 修正 | `src/pyakaikkr/__init__.py` | 変更しない（`pyakaikkr.ase` は明示 import。ase 未導入環境で `import pyakaikkr` が壊れないようにするため） |
| 修正 | `pyproject.toml` | `[options.extras_require] ase = ase>=3.23` を追加 |
| 追加 | `tests/ase/test_inputcard.py`, `tests/ase/test_occupancy.py`, `tests/ase/test_run.py` | §7 |

`pyakaikkr/ase/` はパッケージ名 `ase` と衝突しない（相対 import は `from . import ...`、本家は `import ase` で絶対 import する）。`from __future__ import annotations` を付け、モジュール先頭で `import ase` する。

## 2. 構造の写像（`ase/structure.py`）

関数: `atoms_to_kkr_lattice(atoms) -> dict`

| ASE | AkaiKKR の dict キー | 規則 |
|---|---|---|
| `atoms.cell` (Å) | `brvtyp` | 常に `"aux"` |
| `|cell[0]|` | `a` | bohr 単位。`a = norm(cell[0]) / ase.units.Bohr` |
| `cell[i]` | `r1`, `r2`, `r3` | `cell[i] / norm(cell[0])`（無次元、`a` 単位）。回転はしない |
| `atoms.get_scaled_positions(wrap=True)` | `atmicx` | `["%12.8fa" % x, "%12.8fb" % y, "%12.8fc" % z, type名]`。`a/b/c` 接尾辞は r1, r2, r3 に対する分率座標（`Cif2Kkr` と同じ規約） |
| `len(atoms)` | `natm` | |
| | `c/a`, `b/a`, `alpha`, `beta`, `gamma` | `aux` では inputcard に書かれないが `make_inputcard` の dict には既定値（1.0, 1.0, 90, 90, 90）を入れておく |

検査（違反は `ase.calculators.calculator.InputError`）:

- `atoms.pbc` がすべて True。
- `cell` が 3 本とも非零で右手系（行列式 > 0）。左手系は AkaiKKR の体積計算が負になるので拒否する（自動で並べ替えない）。
- 原子数が 0 でない。

`brvtyp` に `fcc / bcc` などを使う経路は作らない。`a` が慣用胞の格子定数を指すため、ASE の基本胞（`bulk('Cu')` の `cell[0]` 長 2.556 Å）をそのまま渡すと誤るからである。

## 3. 混晶（部分占有）の写像（`ase/occupancy.py`）

### 3.1 入力表現

ASE 標準の 2 つを使う（`ase.io.read` の CIF 読み込みと `ase.spacegroup.crystal(occupancies=...)` が付ける形式）。

- `atoms.info['occupancy']`: `{ "<kind>": { "<元素記号>": 占有率, ... }, ... }`。キーは文字列化した kind 番号。
- `atoms.arrays['spacegroup_kinds']`: 各原子の kind 番号（int 配列）。無い場合は kind = 原子番号 i（0 始まり）とみなす（ASE の CIF リーダーの非対称経路と同じ）。
- `atoms.symbols[i]` は kind の占有 dict に含まれる元素でなければならない（ASE 自身の CIF 書き出しと同じ制約）。

`info['occupancy']` が無い場合は全原子を占有率 1.0 の純物質として扱う。

補助関数 `set_occupancy(atoms, occupancy: dict[int, dict[str, float]]) -> Atoms`
を追加する。原子 index をキーに占有 dict を受け取り、`info['occupancy']` と `spacegroup_kinds` を整合した形で付ける（同じ占有 dict を持つ原子は同じ kind にまとめる）。`bulk('Fe')` などから手で混晶を作るための入口である。

### 3.2 検査（違反は `InputError`）

1. 各原子の kind に対応するキーが `info['occupancy']` にある。
2. 各占有率は `0 < occ <= 1`。
3. kind ごとの占有率の和 s は `s <= 1 + 1e-6`。`s < 1 - 1e-6` のときは空孔成分を追加する（§3.4）。
4. 元素記号は `ElementKKR` の表にある。
5. `atoms.symbols[i]` が kind の占有 dict に含まれる。

### 3.3 type の決定

AkaiKKR では同じ `type` の原子は 1 つのポテンシャルを共有する。したがって type は「kind が同じ」だけでは不十分で、「kind が同じ ∧ 結晶対称性で等価」でなければならない。

手順:

1. kind 番号を 0..K-1 に付け直す（`spacegroup_kinds` が無ければ、占有 dict が同一の原子をまとめて kind とする）。
2. `spglib.get_symmetry_dataset((cell, scaled_positions, kind_ids), symprec=symprec)` を呼ぶ。原子番号の代わりに kind_id を渡すので、占有の違いと対称性の両方が同時に区別される。
3. `dataset.equivalent_atoms` の等価類を type とする。type の並び順は等価類の代表原子の index 昇順。
4. 各原子の `atmicx` 末尾に所属 type 名を書く。

パラメータ `type_mode`:

| 値 | 動作 |
|---|---|
| `"symmetry"`（既定） | 上記手順 |
| `"kind"` | kind をそのまま type にする（対称性を見ない。ユーザーが等価性を保証する場合） |
| `"atom"` | 原子ごとに別 type（`ntyp = natm`。歪んだ超格子など） |

`symprec` の既定値は 1e-5（Å）。spglib が失敗したら `InputError`。

### 3.4 type ごとのパラメータ

| キー | 規則 |
|---|---|
| `type` | 名前。§3.5 |
| `ncmp` | 占有 dict の元素数。占有和 < 1 なら +1（空孔） |
| `anclr` | `ElementKKR(Vc).getAtomicNumber(sym)`。並びは原子番号 Z の昇順（既存テスト `AlMnFeCo_bcc` の inputcard `13, 25, 26, 27` と同じ）。空孔は最後に Z=0 |
| `conc` | `occ * 100.0`。空孔は `(1 - s) * 100.0`。合計が 100 になるよう最後の成分で丸め誤差を吸収する |
| `rmt` | 成分の `getAtomicRMT` の算術平均（`Cif2Kkr` と同じ。表は全元素 0.0 なので実質 0.0 = 自動決定） |
| `field` | 成分の `getAtomicField` の絶対値の最大 |
| `mxl` | 成分の `getAtomicLMax` の最大（表は 2。ユーザー上書きで 3 を指定できる） |

空孔は `Cif2Kkr` と同じく `Vc="Og"` を Z=0 の元素として扱う。ユーザーが `Og` を実元素として使うことは想定しない。

ユーザー上書き: calculator パラメータ `type_params: dict[str, dict]` で type 名または元素記号をキーに `rmt / field / mxl` を上書きできる。適用順は「既定 → 元素記号キー → type 名キー」。

### 3.5 type 名

- 形式: `<組成文字列>_<通し番号>`。組成文字列は Z 昇順に `記号+占有率`（純物質は記号のみ）。例: `Fe_0`, `Fe0.1Ni0.9_0`, `Al0.25Mn0.25Fe0.25Co0.25_0`（この例は Z 順で `Al0.25Mn0.25Fe0.25Co0.25`）。空孔は `Vc0.05` と書く。
- 長さ上限 `len_type = 40`（`akaikkr/source/m_param.f`。cpa2021v01 は README のパッチ適用後に 40）。上限を超えるときは `T<通し番号>` に落とす。
- 空白と `#` を含めない。名前は type 間で一意。

### 3.6 結果の戻し方

- `energy`: 単位胞あたりの全エネルギー。`te [Ry] * ase.units.Rydberg` で eV。
- `magmom`: `get_total_moment` の値（単位胞あたり、μB）。
- `magmoms`: 原子ごとの配列。原子 i の値は所属 type の成分スピンモーメントの濃度加重平均 `Σ conc_c * m_c`。純物質ではその原子のスピンモーメントそのもの。
- `results['component_moments']`: `{type名: {元素記号: {"spin": m_s, "orbital": m_o, "conc": c}}}`。混晶の成分分解値はここで返す（§5 の新関数を使う）。
- `results['converged']`: bool。
- `results['type_of_site']`: `get_type_of_site` の戻り値（inputcard の写像が正しく解釈されたかの検証用）。

## 4. Calculator（`ase/calculator.py`）

```python
class AkaiKKR(FileIOCalculator):
    name = "akaikkr"
    implemented_properties = ["energy", "magmom", "magmoms"]
    _legacy_default_command = "specx < PREFIX.in > PREFIX.out"
    default_parameters = dict(
        go="go", potentialfile="pot.dat",
        edelt=1e-3, ewidth=1.0, reltyp="sra", sdftyp="mjw", magtyp="mag",
        record="init", outtyp="update", bzqlty=6, maxitr=200, pmix=0.02,
        option=None,            # begin_option ... end_option（dict）
        type_mode="symmetry", symprec=1e-5, type_params=None, Vc="Og",
        reuse_potential=True, allow_unconverged=False,
    )
```

- コマンド: ASE の慣例どおり `command=` 引数、環境変数 `ASE_AKAIKKR_COMMAND`、または `~/.config/ase/config.ini` の `[akaikkr] command=`。`PREFIX` は `label` に置換される（`FileIOCalculator.OldShellProfile`）。既定 label は `akaikkr`（ファイルは `akaikkr.in`, `akaikkr.out`）。
- `write_input(atoms, properties, system_changes)`:
  1. `super().write_input` でディレクトリ作成。
  2. §2, §3 で dict を作り、`self.parameters` をかぶせる（`AkaikkrJob.default` を出発点にする）。
  3. `record` の決定（§4.1）。
  4. `AkaikkrJob(self.directory).make_inputcard(dic, open(<prefix>.in, "w"))`。`make_inputcard` は `io.IOBase` を受け付けるのでファイルハンドルで渡す。
  5. 生成 dict を `self._last_dic` に保持する（§4.1 と結果検証に使う）。
- `execute()`: 親クラスのまま（`profile.execute`）。戻り値が非零なら ASE が `CalculationFailed` を出す。
- `read_results()`:
  1. `job = AkaikkrJob(self.directory)`; `job.check_stopped_by_errtrp(<prefix>.out)` が True なら `CalculationFailed`（末尾行 ` ***err in ...` を含める）。
  2. `job.get_convergence` が False で `allow_unconverged=False` なら `SCFError`。
  3. §3.6 の値を `self.results` に入れる。
  4. `job.get_prim_vec` と `job.get_atom_coord` を `_last_dic` の `r1..r3`, `atmicx` と 1e-5 で突き合わせ、不一致なら `CalculationFailed`（構造の解釈ミスを黙って通さない）。
  5. `self.job = job` を残し、`get_dos` などの既存 API に到達できるようにする。

### 4.1 ポテンシャルの再利用（`record`）

- `reuse_potential=True` かつ `<directory>/<potentialfile>` が存在し、かつ前回の `_last_dic` と今回の dict で `ntyp / type / ncmp / anclr / conc` が一致するときだけ `record="2nd"`。それ以外は `"init"`。
- 一致しないポテンシャルで `2nd` を指定すると specx が errtrp で止まるため、構造（格子や座標）は変わっても type 集合が同じときのみ引き継ぐ。EOS の体積スキャンで収束が速くなる。
- ユーザーが `record` を明示した場合はそれを優先する。

### 4.2 `check_state` の上書き

`Calculator.check_state`（`compare_atoms`）は `positions / numbers / cell / pbc / initial_charges / initial_magmoms` しか比べない。`info['occupancy']` と `spacegroup_kinds` は比較対象外なので、占有率だけを変えて `get_potential_energy()` を呼ぶと再計算されない。

```python
def check_state(self, atoms, tol=1e-15):
    changes = super().check_state(atoms, tol)
    if self.atoms is not None and not _same_occupancy(self.atoms, atoms):
        changes = all_changes[:]
    return changes
```

`_same_occupancy` は `info.get('occupancy')` の dict 比較（占有率は 1e-8 で比較）と `spacegroup_kinds` の配列比較。`Atoms.copy()` は `info` を deepcopy し `arrays` も複製するので、`self.atoms` 側に占有情報は残る。

### 4.3 `magtyp` と初期磁気モーメント

- 既定 `magtyp="mag"`（スピン分極）。非磁性計算は `magtyp="nmag"` を明示する。
- `atoms.get_initial_magnetic_moments()` は inputcard に写像しない（AkaiKKR の `field` は初期磁場 [Ry] であり μB ではない）。`field` は `type_params` で与える。ただし初期モーメントが非零の原子があり `magtyp="nmag"` なら警告を出す。

### 4.4 `go` モード

- `go="go"` のみ結果（energy, magmom）を保証する。`fsm`（`fspin` が必要）、`tc`、`dos`、`spc*`（`kpath_raw` が必要）は inputcard の生成だけ行い、結果の解釈は `self.job` の既存 API に任せる。`spc*` は v1 では `InputError` とする。
- `option`（`begin_option` ブロック）は dict をそのまま `make_inputcard` に渡す。

### 4.5 未実装（意図的）

- `forces`, `stress`: `get_forces()` は ASE が `PropertyNotImplementedError` を出す。`calculate_numerical_forces` は既定で使わない（原子数 × 6 回の SCF が必要で実用的でない）。
- `displc`（akaikkr_cnd の格子変位）: v1 では非対応。`dic["displc"]` が無ければ `make_inputcard` は通常形式を書く。
- `brvtyp` の自動判定（fcc/bcc 形式の inputcard 出力）。

## 5. `AkaiKkr.py` の修正

1. `make_inputcard` 内の `print("dic", dic)` を削除する（calculator は inputcard を繰り返し生成するので標準出力が汚れる）。
2. 成分別モーメント取得を追加する。

```python
def get_component_moment(self, outfile) -> list[dict]:
    """*** type-<type>  <El> (z= Z) *** ブロックごとの spin/orbital moment を
    [{"type": str, "element": str, "Z": float, "spin": float, "orbital": float}, ...]
    で返す。出現順は type of site の順、各 type 内は component の順。"""
```

  出力の該当行は `*** type-Mn0.25Al0.25Fe0.25Co0.25_2a_0     Al (z= 13.0) ***` と、その後の `spin moment=  -0.14303  orbital moment=   0.00000`。既存の `get_local_moment` は同じ順序の平坦なリストを返すだけで type や元素の対応を持たないため、ラベル付きの関数を新設する。`get_type_of_site` の `component` 順（`conc` と `anclr`）と zip して濃度加重平均を取る。

3. `ElementKKR.__init__` の `print("debug, _ElementKKR", ...)` も削除する（同じ理由）。
4. （任意）`make_inputcard` の `if "kpath_raw":` は常に真になるバグ。`if "kpath_raw" in dic:` に直す。v1 で `spc*` を拒否するなら影響しないが、直しておくのが自然。

## 6. 依存関係と import

- `setup.cfg`: `[options.extras_require]` に `ase = ase>=3.23`。spglib は既に `install_requires` にある。
- `pyakaikkr/__init__.py` は変更しない。利用側は `from pyakaikkr.ase import AkaiKKR, set_occupancy` と書く。
- `pyakaikkr/ase/__init__.py` は `import ase` の失敗時に、`pip install ase` を促す `ImportError` を出す。

## 7. テスト仕様（`tests/ase/`）

specx 不要（inputcard 生成のみ）:

1. 純物質: `bulk('Cu', 'fcc', a=3.615)` → `brvtyp aux`、`a = 2.5562/0.529177`、`r1..r3` が基本胞、`ntyp=1, ncmp=1, anclr=[[29]], conc=[[100.0]]`。
2. 二元混晶: `set_occupancy(bulk('Ni','fcc',a=3.571), {0: {"Fe":0.1,"Ni":0.9}})` → `ncmp=2, anclr=[[26,28]], conc=[[10.0,90.0]]`（`tests/testrun/akaikkr_cpa2021v01/NiFe/inputcard_go` と成分部が一致）。
3. 四元 HEA: `ase.io.read('tests/testrun/structure/AlMnFeCo-Im3m.cif')` → `ntyp=1, ncmp=4, anclr=[[13,25,26,27]], conc=[[25.0]*4]`（`AlMnFeCo_bcc/inputcard_go` と一致）。
4. 空孔: 占有和 0.95 → `ncmp` が +1、最後が `anclr=0, conc=5.0`。
5. 対称性による type 分割: 2×2×2 超格子で 1 原子だけ z 方向に 0.05 Å ずらす → `type_mode="symmetry"` で `ntyp` が 2 以上、`type_mode="kind"` で 1。
6. 検査: 占有率 1.2、和 1.1、未知元素、`symbols[i]` が占有 dict に無い、左手系 cell、`pbc=False` がそれぞれ `InputError`。
7. `check_state`: 同じ `Atoms` で占有率だけ変えると `all_changes` が返る。
8. `record`: 同じ type 集合で 2 回目は `2nd`、組成を変えると `init` に戻る。

specx 必要（`akaikkr` 環境、`tests/testrun/akaikkr` の Cu 入力と同じ条件）:

9. Cu の `get_potential_energy()` が `tests/testrun/akaikkr/reference/ifort.json` の `Cu_go` の `te` に対して `Rydberg` 換算で 1e-6 eV 以内。
10. NiFe の `magmom` が参照値と 1e-3 μB 以内、`component_moments` が Fe と Ni の 2 成分を持つ。
11. 出力の `primitive translation vectors` と `atoms in the unit cell` が入力と一致する（§4 の検証を通る）。
12. EOS: Cu の体積 5 点を `ase.eos.EquationOfState` で当てはめ、2 回目以降の呼び出しで `record=2nd` が使われて反復回数が減る。

## 8. 使用例

```python
from ase.build import bulk
from ase.io import read
from pyakaikkr.ase import AkaiKKR, set_occupancy

# 純物質
cu = bulk("Cu", "fcc", a=3.615)
cu.calc = AkaiKKR(directory="cu", command="specx < PREFIX.in > PREFIX.out",
                  sdftyp="pbeasa", bzqlty=8)
print(cu.get_potential_energy())        # eV / cell

# 二元混晶
nife = set_occupancy(bulk("Ni", "fcc", a=3.571), {0: {"Fe": 0.1, "Ni": 0.9}})
nife.calc = AkaiKKR(directory="nife", magtyp="mag", type_params={"Fe0.1Ni0.9_0": {"mxl": 3}})
print(nife.get_magnetic_moment(), nife.calc.results["component_moments"])

# HEA（CIF の部分占有をそのまま使う）
hea = read("tests/testrun/structure/AlMnFeCo-Im3m.cif")   # info['occupancy'] 付き
hea.calc = AkaiKKR(directory="hea")
print(hea.get_potential_energy())
```

## 9. 実装メモ（2026-09-24 実装時の仕様からの変更点）

- `record` の既定値は `None`（自動判定: ポテンシャルファイルがあり `ntyp/type/ncmp/anclr/conc/magtyp` が前回と同じなら `2nd`、それ以外は `init`）。文字列を与えれば固定。
- bohr 換算は ASE の `units.Bohr` ではなく `pyakaikkr.Unit().length_au2ang`（0.529177）を使う。`Cif2Kkr` と同じ値にして両経路の `a` を揃えるため。
- 基本胞への縮約 `to_primitive(atoms, symprec)` を追加した。`spglib.find_primitive` に kind 番号を原子種として渡すので占有の違うサイトは併合されない。縮約後の原子は「対応する入力原子の最小 index」順に並べる。
- type 順は「原子の出現順」。ASE は CIF の `_atom_site` 順を保つが、pymatgen の慣用胞は電気陰性度順に並べ替えるため、cif 経路と type 順が一致しない構造がある（FeRh0.5Pt0.5、SmCo5）。type 順に依存する設定は type 名で決めること（`testscript_ase_spec.md` §3.2 の SmCo5 の例）。
- 例外 `KKRStructureMismatchError` を `pyakaikkr.Error` に追加し、`check_kkr_output_structure(job, outfile, param)` を `pyakaikkr.ase.structure` に置いた（calculator とテストスクリプトの geom 検証で共用）。
- `pyakaikkr` は編集可能インストール（`pip install -e .`）に切り替えた。
- pytest は `tests/ase/` に置いた。specx を使うテストは環境変数 `AKAIKKR_PROGRAM_PATH`（`akaikkr/specx` を含むディレクトリ）が無ければ skip する。

- **type 順（追記）**: `rmt=0`（自動決定）のマフィンティン半径は type の並び順に依存する（FeB1.95 で Fe→B と B→Fe で Fe の rmt が 0.381 と 0.474、全エネルギーが 0.02 Ry 違った）。そのため `atoms_to_kkr_types` に `type_order` を追加し、既定を `"electronegativity"`（占有率加重の Pauling 電気陰性度の昇順、同点は出現順、Og などの NaN は最後）にした。これは pymatgen がサイトを並べる規則と同じで、cif 経路と type 順が一致する（consistency テストで順序一致も検査する）。`"appearance"` で出現順にできる。
- **原子の並び（追記）**: 自動 rmt は原子の並びにも依存する（SmCo5 で type 順を揃えても原子が Co, Co, Sm の順だと rmt が変わった）。`atoms_to_kkr_atmicx` は原子を type 順にまとめて並べる（Cif2Kkr と同じ）。`param["atom_order"]` に行ごとの元の原子 index を持つ。calculator の `magmoms` は元の原子順で返す。

## 10. 後日談（2026-09-24、AiiDA プラグイン整備時に分かったこと）

- **pymatgen の ASE 変換は部分占有を受け付けない。** `pymatgen.io.ase.AseAtomsAdaptor.get_atoms` は部分占有サイトがあると `ValueError: ASE Atoms only supports ordered structures` を投げる（pymatgen 2026.9.24）。本仕様 §3.1 の ASE 表現（`info["occupancy"]` + `spacegroup_kinds`）に pymatgen の `Structure` から変換する経路は無いので、CPA 構造は CIF を `ase.io.read` で読むか `set_occupancy` で作る。AiiDA では `StructureData(pymatgen=...)` が Kind の weights として占有を保持する（aiida-akaikkr のパーサーはこれを使う）。
- 利用者向けの要点は [ase_calculator_usage.md](ase_calculator_usage.md) にまとめた。
- 参照ファイルは `reference/ifort_ase.json` として cif 経路と分けた（Cu の aux/fcc の差 1.8e-11 が閾値 1e-11 を超えるため）。
