# AkaiKKRTestScript と tests/testrun/akaikkr/testrun.py の ASE 対応 修正仕様

作成日: 2026-09-24
前提: `docs/ase_calculator_spec.md`（pyakaikkr 側の `pyakaikkr.ase` パッケージ）が実装済みであること。本仕様はその上に載るテストスクリプト側の変更を定める。

## 0. 方針

- 既存の動作を変えない。`python testrun.py <program_path> --set Cu` は今までどおり pymatgen（`Cif2Kkr` / `CompareCifKkr`）経路で動き、既存の参照 `reference/ifort.json` と比較できる。
- ASE 経路は「構造パラメータ dict を作る部分」だけを差し替える。`GoGo` 系クラス、`OutputAnalyzer`、参照 JSON の形式、`ExeUtil` / `ResultUtil` は共通のまま使う。
- ASE 版の実行スクリプトは `tests/testrun/akaikkr/testrun_ase.py` とし、`testrun.py` は残す。
- バグ修正（§5）は両経路に効く。修正により既存参照と一致しなくなる項目は明記する。
- ase は必須依存にしない。`akaikkr_testscript` は ase 未導入でも import・実行できる。

## 1. 変更ファイル一覧

| 区分 | パス | 内容 |
|---|---|---|
| 追加 | `src/akaikkr_testscript/asestruc.py` | `get_kkr_struc_from_ase` と geom 検証（§2） |
| 修正 | `src/akaikkr_testscript/testrun_class.py` | 構造取得の振り分け関数 `get_kkr_struc`、各 `_X_common_param` の呼び出し置換、`change_atomic_type` の index キー対応、`all_go` の `backend` 引数（§3） |
| 修正 | `src/akaikkr_testscript/__init__.py` | 変更なし（`asestruc` は `get_kkr_struc` 内で遅延 import） |
| 修正 | `pyproject.toml` | `[options.extras_require] ase = pyakaikkr[ase]` |
| 修正 | `src/pyakaikkr/Fmg.py` | flip 名不一致の検出（§5.1） |
| 修正 | `src/pyakaikkr/GoGo.py` | `GoFmg` の flip 名解決（§5.1） |
| 追加 | `tests/testrun/akaikkr/testsets.py` | `make_exe()` を `testrun.py` から移す（§4） |
| 修正 | `tests/testrun/akaikkr/testrun.py` | `testsets.make_exe` を使う。動作は同じ（§4） |
| 追加 | `tests/testrun/akaikkr/testrun_ase.py` | ASE 経路の実行スクリプト（§4） |
| 追加 | `tests/testrun/akaikkr/reference/<compiler>_ase.json` | ASE 経路の参照（`--create_ref` で生成、§4.3） |
| 追加 | `tests/ase/test_backend_consistency.py` | cif 経路と ASE 経路の dict 一致テスト（§6） |

`tests/testrun/akaikkr_cpa2021v01/` と `tests/testrun/akaikkr_cnd/` の `testrun.py` は本仕様の範囲外。同じ手順（§4）でそれぞれ `testrun_ase.py` を足せる。

## 2. `asestruc.py`

```python
def get_kkr_struc_from_ase(structure, akaikkr_exe: str, displc: bool,
                           use_bravais=False, remove_temperaryfiles=True,
                           Vc: str = "Og", directory: str = "temporary",
                           fmt=None, type_mode="symmetry", symprec=1e-5,
                           check_geom=True) -> dict
```

戻り値は `get_kkr_struc_from_cif` と同じキーを持つ dict:
`brvtyp, a, c/a, b/a, alpha, beta, gamma, r1, r2, r3, ntyp, type, ncmp, rmt, field, mxl, anclr, conc, natm, atmicx`。
`displc` は受け取るだけで使わない（`GoGo.execute` が `anclr` から `displc` を作るため。cif 経路と同じ）。

### 2.1 入力

- `structure` は `ase.Atoms` またはファイルパス。パスのときは `ase.io.read(structure, format=fmt)` で読む。CIF は ASE 既定で `fractional_occupancies=True` なので `info['occupancy']` と `spacegroup_kinds` が付く。
- `tests/testrun/structure/GaAsVc-F43m.cif` の `Og`（空孔）は ASE が元素 Og として読む。`ElementKKR(Vc="Og")` が Z=0 に写すので追加処理は不要。
- `use_bravais` は署名互換のためだけに受け取る。ASE 経路は常に `brvtyp=aux`（`ase_calculator_spec.md` §2）。`True` が渡されたら一度だけ警告を出して無視する。

### 2.2 変換

`pyakaikkr.ase.structure.atoms_to_kkr_lattice(atoms)` と `pyakaikkr.ase.occupancy.atoms_to_kkr_types(atoms, Vc=Vc, type_mode=type_mode, symprec=symprec)` を呼び、結果を 1 つの dict にまとめる。この 2 関数は calculator と共有する（テストスクリプト専用の変換ロジックを持たない）。

順序の規則（既存参照の `spinlocalmoment` などはリストの並びで比較されるため、cif 経路と揃える）:

| 対象 | 規則 | cif 経路との対応 |
|---|---|---|
| type の並び | `Atoms` 内で最初に現れる原子の index 昇順 | `Cif2Kkr` は慣用胞のサイト出現順で type を作る。CIF を読んだ `Atoms` は CIF の `_atom_site` 順なので一致する |
| type 内の成分の並び | 原子番号 Z 昇順 | 既存 inputcard は `AlMnFeCo_bcc` が 13,25,26,27、`NiFe` が 26,28 でいずれも Z 昇順 |
| `atmicx` の並び | `Atoms` の原子順 | |

`SmCo5_oc` の `param["mxl"] = [3, 2, 2]` のように type 順に依存する設定があるので、§6 のテストで type 順の一致を確認する。

### 2.3 geom 検証（`check_geom=True`）

cif 経路の `CompareCifKkr.convert_and_compare` に相当する検査を行う。

1. `directory` に `AkaikkrJob` で `go="geom"` の inputcard を書き、`akaikkr_exe` を実行する（`inputcard_geom`, `out_geom.log`。cif 経路と同じファイル名）。
2. `get_prim_vec(out, unitof="relative")` を `r1..r3` と、`get_atom_coord(out)` を `atmicx` の分率座標 × (r1,r2,r3) と比較する。許容差 1e-6（a 単位）。原子座標は格子ベクトルの整数倍を法として比べる。type 名も一致すること。
3. 不一致なら `CIF2KKRNsiteInconsistentError` か `CIF2KKRCellShapeError`（既存の例外を流用）を投げる。
4. `remove_temperaryfiles=True` なら inputcard と出力を削除し、空なら `directory` も消す（cif 経路と同じ後始末）。

`check_geom=False` は specx を持たない環境で dict だけ作るための逃げ道（§6 のテストで使う）。

## 3. `testrun_class.py` の修正

### 3.1 振り分け関数

```python
def get_kkr_struc(ciffilepath: str, akaikkr_exe: dict, displc: bool,
                  use_bravais=True, remove_temperaryfiles=True,
                  Vc: str = "Og", directory: str = "temporary") -> dict:
    backend = akaikkr_exe.get("backend", "cif")
    if backend == "cif":
        return get_kkr_struc_from_cif(ciffilepath, akaikkr_exe["specx"], displc,
                                      use_bravais=use_bravais,
                                      remove_temperaryfiles=remove_temperaryfiles,
                                      Vc=Vc, directory=directory)
    elif backend == "ase":
        from .asestruc import get_kkr_struc_from_ase   # ase はここで初めて import
        return get_kkr_struc_from_ase(ciffilepath, akaikkr_exe["specx"], displc,
                                      use_bravais=use_bravais,
                                      remove_temperaryfiles=remove_temperaryfiles,
                                      Vc=Vc, directory=directory)
    raise ValueError(f"unknown backend {backend}")
```

- `akaikkr_exe` は `all_go` が作る `prog` dict（`specx`, `fmg`, `args`）。ここに `backend` を足す。キーが無ければ `"cif"` なので、`prog` を自前で作る既存の利用者コードも今までどおり動く。
- `get_kkr_struc_from_cif` は削除・変更しない（公開名として残す）。

### 3.2 `_X_common_param` の置換

13 個の `_X_common_param`（Cu, Fe, Co, Ni, AlMnFeCo_bcc, FeRh05Pt05, NiFe, Fe_lmd, FeB195, GaAs, Co2MnSi, SmCo5_oc, SmCo5_noc）で

```python
param = get_kkr_struc_from_cif(
    ciffilepath=ciffilepath, akaikkr_exe=akaikkr_exe["specx"], displc=displc, ...)
```

を

```python
param = get_kkr_struc(
    ciffilepath=ciffilepath, akaikkr_exe=akaikkr_exe, displc=displc, ...)
```

に置き換える。他の行（`magtyp`, `rmt`, `mxl` などの上書き、`Fe_lmd` の `ncmp/anclr/conc` の二重化、`GaAs` の `ncmp` 再設定）は触らない。これらは dict のキーだけに依存し、backend に依存しない。

### 3.3 `change_atomic_type`

`type_rep` のキーに type 名（str）に加えて type の index（int）を許す。

```python
def change_atomic_type(param, type_rep):
    name_rep = {}
    for k, v in type_rep.items():
        if isinstance(k, int):
            name_rep[param["type"][k]] = v
        else:
            name_rep[k] = v
    ...  # 以降は name_rep で従来どおり置換。未知の名前は従来どおり KeyError
```

`_Cu_common_param` の `change_atomic_type(param, {"Cu_4a_0": "Cu"})` を `{0: "Cu"}` に変える。cif 経路では `Cu_4a_0` が唯一の type なので結果は同じ。ASE 経路の type 名（`Cu_0`）でも動く。

### 3.4 `all_go`

```python
def all_go(akaikkr_exe, fmg_exe, exe_dic, displc=False, backend="cif"):
```

- `prog = {"specx": ..., "fmg": ..., "args": args, "backend": backend}`。
- 参照ファイル名と結果ファイル名に backend の接尾辞を付ける。

| backend | 参照 | 結果 |
|---|---|---|
| `cif` | `reference/<compiler>.json`（従来どおり） | `result.json`（従来どおり） |
| `ase` | `reference/<compiler>_ase.json` | `result_ase.json` |

- ASE 経路の参照を cif 経路と分ける理由: ASE 経路は `aux` + `r1..r3` で格子を渡すため、specx 内部の対称性処理が `fcc/bcc` 指定と厳密には同じ経路を通らず、`rdiff_te = 1e-11` の閾値で一致する保証がない。両経路の突き合わせは §6 の consistency テストで dict のレベルで行う。
- `meta` には `"backend"` を追加する（`make_meta` の戻り値に `meta["backend"] = backend`）。
- `parse_args` は変更しない。backend の切り替えはスクリプト名（`testrun.py` / `testrun_ase.py`）で行う。

## 4. `tests/testrun/akaikkr/` のスクリプト

### 4.1 `testsets.py`（追加）

現在の `testrun.py` の `make_exe()` をそのまま移す。`from akaikkr_testscript import *` を先頭に置く。`exe_dic` の内容は変えない。

### 4.2 `testrun.py`（修正）

```python
from akaikkr_testscript import all_go
from testsets import make_exe

if __name__ == "__main__":
    all_go("specx", fmg_exe="fmg", exe_dic=make_exe(), displc=False)
```

コマンドラインと出力は従来と同じ。`if "pyakaikkr" not in sys.modules: from pyakaikkr import *` の行は不要になるので削除する（`akaikkr_testscript` が `pyakaikkr` を import する）。

### 4.3 `testrun_ase.py`（追加）

```python
from akaikkr_testscript import all_go
from testsets import make_exe

if __name__ == "__main__":
    all_go("specx", fmg_exe="fmg", exe_dic=make_exe(), displc=False, backend="ase")
```

- 実行は `python testrun_ase.py <program_path> --set Cu` 。初回は `--create_ref` で `reference/ifort_ase.json` を作る（`meta.json` の扱いは従来どおり）。
- `--set` のキーと関数は `testsets.py` 経由で `testrun.py` と共通。
- ase 未導入なら `get_kkr_struc` の遅延 import で `ImportError` が出る。メッセージに `pip install ase` を含める。
- 実行ディレクトリの基底名（`akaikkr`）から `specx` の場所を決める仕組み（`kkrtype = basename(cwd)`）は共通なので、`testrun_ase.py` も `tests/testrun/akaikkr/` で実行する。

## 5. バグ修正（両経路に効く）

### 5.1 `GoFmg` の flip が黙って無視される

`AlMnFeCo_bcc_gofmg` は `flip_list = ["HEA_Mn_25.0%"]` を渡すが、`typeofsites` の shortname は `Mn0.25Al0.25Fe0.25Co0.25_2a_0_Mn_25.0%` なので一致せず、`Fmg.make_inputfile` は `all_flip_dic["HEA_Mn_25.0%"] = 1` を追加するだけで何も反転しない。実際に `tests/testrun/akaikkr_cpa2021v01/AlMnFeCo_bcc/fmg.input` は `pot.dat 1 2 3 4 / pot_fmg.dat 1 2 3 4` で反転無し。

修正:

1. `Fmg.make_inputfile`: `flip_list` の要素が `comp_shortname_list` に無ければ `ValueError`（候補一覧をメッセージに含める）。
2. `GoFmg.__init__` の `flip_list` に「元素と濃度だけの指定」を許す。要素が `"*"` で始まる場合は shortname の末尾一致（`shortname.endswith(entry[1:])`）で解決し、複数一致すれば全部反転する。`prescript` で `typeofsites` を得た後に解決して `Fmg` に渡す。
3. `AlMnFeCo_bcc_gofmg` の `flip_list` を `["*_Mn_25.0%"]` にする。type 名が backend で変わっても効く。

影響: 修正後は Mn が実際に反転するので `AlMnFeCo_bcc_gofmg` の `te`, `tm`, `spinlocalmoment` は既存参照 `reference/ifort.json` と一致しなくなる。cif 経路の参照はこの項目だけ作り直す（`--create_ref` は参照ファイル全体を作るので、既存 JSON の `AlMnFeCo_bcc_gofmg` を新しい結果で差し替える手順を README に書く）。

### 5.2 `make_inputcard` の `if "kpath_raw":`

`ase_calculator_spec.md` §5 の項目。常に真になるため `spc*` 以外でも `dic["kpath_raw"]` があれば書き出す。`if "kpath_raw" in dic:` に直す。既存テストは `spc*` のときだけ `kpath_raw` を入れるので結果は変わらない。

### 5.3 デバッグ print

`AkaikkrJob.make_inputcard` の `print("dic", dic)` と `ElementKKR.__init__` の `print("debug, ...")` を削除する。テスト出力が読みやすくなる以外の影響は無い。

修正しないもの: README の既知バグ「`TEST FAILED` が常に表示される」「`Awk_both.png` がトップに生成される」は本仕様の範囲外（再現条件を確認してから別途）。

## 6. テスト（`tests/ase/test_backend_consistency.py`）

specx が必要なので環境変数 `AKAIKKR_PROGRAM_PATH`（`testrun.py` の `program_path` と同じ意味）が無ければ skip する。

各 CIF（`tests/testrun/structure/*.cif` のうち `testrun.py` で使う 12 件。`Fe_lmd` は `Fe-Im3m.cif` を共用）について:

1. `get_kkr_struc_from_cif(...)` と `get_kkr_struc_from_ase(...)` を同じ `directory` 設定で呼ぶ。
2. 次を比較する。

| キー | 判定 |
|---|---|
| `ntyp`, `natm` | 等しい |
| `ncmp` | 等しい |
| `anclr`, `conc` | type ごとに等しい（成分は Z 昇順に並べてから比較。cif 経路が Z 昇順でない CIF があればここで判明する） |
| `type` の順 | cif 経路の type i と ASE 経路の type i が同じ原子集合を指す（`atmicx` の type 名を index に直して比較） |
| 格子 | `a * r_i`（bohr）で張る格子が同じ格子である（cif 経路は `fcc` などなので `_TranslationKKR.getMatrix(brvtyp)` で `r1..r3` に展開してから比較。格子ベクトルの整数結合で一致すればよい） |
| `atmicx` | 同じ格子上で同じ位置集合（並びは type 順の一致を確認したうえで比較） |

3. `SmCo5_oc` は `mxl=[3,2,2]` が Sm, Co(2c), Co(3g) の順に対応することを個別に確認する。

specx 実行を伴う確認（手動）:

4. `python testrun_ase.py <program_path> --set Cu --create_ref` → `reference/ifort_ase.json` ができる。
5. `python testrun_ase.py <program_path> --set Cu` → `ALL TESTS PASSED`。
6. `python testrun.py <program_path> --set Cu` → 従来どおり `ALL TESTS PASSED`（既存参照、回帰確認）。
7. `ifort.json` と `ifort_ase.json` の `Cu_go` の `te` の差を記録する。`rdiff_te = 1e-11` を超えるなら §3.4 の「参照を分ける」判断が正しかったことになり、超えないなら将来 参照を共有する検討材料になる。
8. `--set AlMnFeCo` を両経路で回し、§5.1 の修正後に `gofmg` の `tm` が `go` と異なる（Mn が反転している）ことを確認する。

## 7. 実装順

1. §5.3, §5.2（無害な修正）
2. `asestruc.py` と §3.1〜3.3（cif 経路の回帰: `testrun.py --set Cu` が通ること）
3. §3.4 と §4（`testrun_ase.py --set Cu --create_ref` → 比較）
4. §6 の consistency テスト、12 件の CIF で通す
5. §5.1（参照の差し替えを伴うので最後）

## 8. 実装メモ（2026-09-24 実装時の仕様からの変更点）

- §2.2 の「type 順は CIF の順で cif 経路と一致する」は誤りだった。pymatgen の慣用胞はサイトを電気陰性度順に並べ替えるため、FeRh0.5Pt0.5（cif: Fe, RhPt / ase: RhPt, Fe）と SmCo5（cif: Sm, Co, Co / ase: Co, Co, Sm）で順が異なる。対処として `_SmCo5_oc_common_param` と `_SmCo5_noc_common_param` の `mxl = [3, 2, 2]` を「type 名が `Sm` で始まれば 3、それ以外 2」に変えた。両経路で同じ結果になる。
- `get_kkr_struc_from_ase` に `use_primitive=True` を追加した。ASE の CIF 読み込みは慣用胞を返すため、`pyakaikkr.ase.to_primitive` で基本胞に落としてから変換する（cif 経路の `use_primitive=True` に相当）。
- §6 の consistency テストは type の対応を「組成と位置集合の一致」で取り、順序には依存しない。
- 参照ファイル名の接尾辞は `_ref_suffix(backend)`（`cif` → 無し、それ以外 → `_<backend>`）。`meta["backend"]` を記録する。
- `pyakaikkr.Fmg.make_inputfile` は未知の flip 名で `ValueError`、`GoGo._resolve_flip_list` が `"*..."` を末尾一致で解決する。`AlMnFeCo_bcc_gofmg` の `flip_list` は `["*_Mn_25.0%"]`。

- **type 順（追記）**: `rmt=0`（自動決定）のマフィンティン半径は type の並び順に依存する（FeB1.95 で Fe→B と B→Fe で Fe の rmt が 0.381 と 0.474、全エネルギーが 0.02 Ry 違った）。そのため `atoms_to_kkr_types` に `type_order` を追加し、既定を `"electronegativity"`（占有率加重の Pauling 電気陰性度の昇順、同点は出現順、Og などの NaN は最後）にした。これは pymatgen がサイトを並べる規則と同じで、cif 経路と type 順が一致する（consistency テストで順序一致も検査する）。`"appearance"` で出現順にできる。
- **原子の並び（追記）**: 自動 rmt は原子の並びにも依存する（SmCo5 で type 順を揃えても原子が Co, Co, Sm の順だと rmt が変わった）。`atoms_to_kkr_atmicx` は原子を type 順にまとめて並べる（Cif2Kkr と同じ）。`param["atom_order"]` に行ごとの元の原子 index を持つ。calculator の `magmoms` は元の原子順で返す。

## 9. 後日談（2026-09-24）

- `tests/testrun/akaikkr_cnd` を全件流した際に直した不具合（`*_j30` 後処理の type 名決め打ち、`_canonical_type_name` による Jij `pair` キーの正規化）と、残る 6 件の最終桁差は [testscript_usage.md](testscript_usage.md) §4–5 にまとめた。
- `AlMnFeCo_bcc_gofmg` は flip 名の修正で本当に Mn を反転するようになったため、旧 cif 参照とは一致しない。参照を作り直すこと。
- 同じ `_<物質>_common_param` を AiiDA から使う経路を `aiida-akaikkr/example/run_examples.py` に置いた。共通パラメータを構造キー（`brvtyp, a, ..., atmicx, displc`）とそれ以外に分けて CalcJob の `structure` / `parameters` に渡す。akaikkr_cnd では `displc` が無いと specx が "illegal input" で止まるので、`GoGo.execute` と同様に `make_displc_list(anclr)` を付ける。
