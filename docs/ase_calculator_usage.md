# pyakaikkr.ase の使い方（ASE calculator）

設計の詳細は [ase_calculator_spec.md](ase_calculator_spec.md) にあります。ここは利用者向けの要点です。

## 1. インストール

```bash
pip install -e library/PyAkaiKKR[ase]      # ase >= 3.23 を一緒に入れる
```

specx の起動コマンドは次のいずれかで与えます。`PREFIX` は calculator の label に置き換わります。

- `AkaiKKR(command="/path/to/specx < PREFIX.in > PREFIX.out")`
- 環境変数 `ASE_AKAIKKR_COMMAND`
- ASE の設定ファイルの `[akaikkr]` セクション

OpenMP スレッド数は `OMP_NUM_THREADS` で指定します。ifort ビルドはログインシェル以外だとランタイムライブラリが見つからないことがあるので、必要なら `source /opt/intel/oneapi/setvars.sh` をコマンドに含めてください。

## 2. 純物質

```python
from ase.build import bulk
from pyakaikkr.ase import AkaiKKR

cu = bulk("Cu", "fcc", a=3.615)
cu.calc = AkaiKKR(directory="cu", command="/path/to/specx < PREFIX.in > PREFIX.out", magtyp="nmag")
e = cu.get_potential_energy()        # eV / cell
```

得られる量は `energy`（全エネルギー、eV / セル）、`magmom`（セルの全モーメント）、`magmoms`（原子ごと。CPA では成分の濃度加重平均）です。力と応力はありません。`calc.results["component_moments"]` に type ごと・成分ごとのモーメントが入ります。

## 3. 混晶（CPA）

ASE 標準の部分占有表現（`atoms.info["occupancy"]` と `atoms.arrays["spacegroup_kinds"]`）を使います。`set_occupancy` で作るか、部分占有を含む CIF を `ase.io.read` で読めばそのまま使えます。

```python
from pyakaikkr.ase import AkaiKKR, set_occupancy

nife = set_occupancy(bulk("Ni", "fcc", a=3.571), {0: {"Fe": 0.1, "Ni": 0.9}})
nife.calc = AkaiKKR(directory="nife", command=..., sdftyp="pbeasa", bzqlty=8)
nife.get_potential_energy(); nife.get_magnetic_moment()
```

- type は「占有が同じ、かつ対称的に等価（spglib）」なサイトをまとめたものです（`type_mode="symmetry"`）。`type_mode="kind"` なら占有だけで分けます。
- 空孔は `Vc="Og"`（既定）のダミー元素で表します。
- 慣用胞を渡した場合は `to_primitive(atoms)` で基本胞に落とせます。占有の違うサイトは併合されません。

## 4. 構造の渡し方と単位

- 格子は常に `brvtyp=aux`、`a=|cell[0]|`（bohr）、`r1..r3 = cell / a` で渡します。bohr 換算は `pyakaikkr.Unit().length_au2ang`（0.529177）で、`Cif2Kkr` と同じ値です。
- `rmt` を `type_params` で与えるときは、この `a` を単位にした値です。出力の `rmt` も `a` 単位なので、慣用胞 `a` で走らせた cif 経路とは fcc で √2 倍などの見かけの差が出ます（物理的な半径は同じ）。
- `rmt=0`（自動）のときマフィンティン半径は **type の並び順と原子の並び順の両方**に依存します。既定の `type_order="electronegativity"`（占有加重 Pauling 電気陰性度の昇順）と、原子を type ごとにまとめる並びは pymatgen の `Cif2Kkr` と同じ規則で、両経路の全エネルギーが一致します。`type_order="appearance"` にすると出現順になり、結果が変わることがあります。

## 5. パラメータ

`AkaiKKR(**kwargs)` の kwargs は inputcard のキーがそのまま使えます。既定値:

| キー | 既定 | 備考 |
|---|---|---|
| go | go | `dos` などは未対応（エネルギー計算のみ） |
| edelt / ewidth | 1e-3 / 1.0 | |
| reltyp / sdftyp | sra / mjw | |
| magtyp | mag | 初期磁気モーメントは AkaiKKR 側の既定 |
| record | None | 自動判定。ポテンシャルがあり type 構成が前回と同じなら `2nd`、それ以外は `init`。文字列で固定可 |
| bzqlty / maxitr / pmix | 6 / 200 / 0.02 | |
| option | None | `begin_option` ブロック（[akaikkr_option_keys.md](akaikkr_option_keys.md)） |
| type_params | None | type 名 → `{"rmt":..., "mxl":..., "field":...}` |
| reuse_potential | True | `pot.dat` を作業ディレクトリに残して再利用 |
| allow_unconverged | False | 未収束なら例外 |

## 6. 落とし穴

- **pymatgen の ASE 変換は部分占有を受け付けない。** `AseAtomsAdaptor.get_atoms` は `ValueError: ASE Atoms only supports ordered structures` を出します（pymatgen 2026.9 で確認）。CPA 構造を pymatgen から持ってくるときは `Structure` を経由せず、CIF を `ase.io.read` で読むか `set_occupancy` を使ってください。AiiDA では `StructureData(pymatgen=...)` が Kind の weights として占有を保持します。
- `check_state` は上書きしてあり、`info["occupancy"]` や `spacegroup_kinds` の変化も再計算の引き金になります。ASE 標準の `compare_atoms` はこれらを見ません。
- 出力構造が入力と食い違うと `KKRStructureMismatchError` になります（`check_kkr_output_structure`）。
- 例外 `KKRFailedExecutionError` は specx が異常終了したとき、`KKRValueAquisitionError` は出力から値が取れないときです。

## 7. テスト

```bash
export AKAIKKR_PROGRAM_PATH=/path/to/AkaiKKRprogram.2022.0721.ifort   # akaikkr/specx を含む
pytest tests/ase
```

specx を使うテストは `AKAIKKR_PROGRAM_PATH` が無いと skip されます。`test_backend_consistency.py` は cif 経路と ase 経路の inputcard と全エネルギーの一致を検査します。
