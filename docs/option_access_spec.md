# `begin_option` を pyakaikkr から参照する仕様

作成日: 2026-09-25
対象: pyakaikkr 2023.2.1（`library/PyAkaiKKR`）、AkaiKKR 2022.0721（akaikkr, akaikkr_cnd, akaikkr_cpa2021v01）
根拠: `akaikkr_common/source/m_optn.f`（`optnrd_rd`, `optnwrt_*`）、`docs/akaikkr_option_keys.md`、2026-09-25 に Cu で option 付き go を 1 回走らせて採った出力

## 0. 目的と範囲

いま pyakaikkr は `begin_option` を**書く**ことしかできない（`AkaikkrJob.make_inputcard` の `make_option_card` が `dic["option"]` を書き出す）。次の 4 つの「参照」をできるようにする。

| 記号 | 参照したいもの | いまの状態 |
|---|---|---|
| A | キーの一覧（正式名・別名・型・既定値・効くコード・意味）を Python から引く | `docs/akaikkr_option_keys.md` の表にしか無い |
| B | 既存の inputcard の `begin_option` ブロックを dict に戻す | 無い |
| C | 出力ログから、specx が実際に受け取った値（`optnwrt:` 行）と、実効値（ヘッダの `meshr mse ng mxl` 表）を読む | 無い |
| D | 書く前に検証する（未知キー、別名の正規化、`key= value` の書式、論理値の表記） | 無い。未知キーは specx が `unknown token` で停止して初めて分かる |

対象外: `klabel` の意味付け（この版の Fortran では保持と表示だけ）、inputcard 全体のパーサ（option ブロックだけを読む）、aiida-akaikkr の parser への組み込み（§7 に接続点だけ書く）。

## 1. specx 側の事実（実装が前提にするもの）

`m_optn.f` と実測から。

- 文法: `begin_option` から `end_option` までを空白区切りトークンで読む。スカラは `key=` トークンの**次のトークン**が値（`mse= 5` は可、`mse=5` は `unknown token: mse=5` を出して `stop 100`）。配列は `begin_<name>` から `end_<name>` までのトークン列。`begin_option` が無ければ何もしない。
- スカラキー（`t_optparam` の 15 個、値はすべて `character(80)` で保持、使う側で `read(...,*)` する）: `mse`, `tol`, `ng`, `ie`, `critic`, `rkick`, `dex`, `cemesr_ref`, `spmain_bnd2`, `cpaitr_show`, `cpaitr_tol`, `spckkr_dmpc0`, `spckkr_itrmx`, `ddos`, `tempmu`。別名: `number_emesh`→`mse`、`thresh_scf`/`thresh_go`→`tol`、`ndegree_cheb`→`ng`、`ndirection_cnd`→`ie`。
- 配列キー: `klabel` だけ処理される。他の `begin_xxx` は `unknown array labelxxx, but ignore it.` を出して続行する（`label` と名前の間に空白が無いのは Fortran 側の書式）。
- 出力への echo: 値が与えられたキーが**使われた時点**で ` optnwrt:<name> <value>` が 1 行出る（`optnwrt_wt` で全部をまとめて出す機能は呼ばれていない）。したがって
  - `tol`, `ng`, `dex`, `ie`, `mse` は `specx.f` の冒頭で読まれ、ヘッダの前に出る。
  - `critic` は SCF 反復の直前（`spmain.f`、`show="once"` なので 1 回だけ）。`rkick` は kick のとき。`spmain_bnd2` は与えたときだけ ` optnwrt:spmain_bnd2 b` と出る（`debug bnd2_type a` の行は option の有無によらず常に出る別物）。
  - `cemesr_ref` は `cemesr.f` が呼ばれる `dos`/`spc` 等でしか出ない。go の出力には出ない。**echo 名は `ref`**（` optnwrt:ref  0.750000000000000`）で、変数名 `cemesr_ref` とは違う（実装後に akaikkr_cnd の dos で判明。`OptionKey.echo_name` で対応）。他のキーは変数名と同じ。
  - `klabel` は echo されない。
  - 書式は list-directed 出力: ` optnwrt:mse           5`, ` optnwrt:tol  1.000000000000000E-005`, ` optnwrt:critic  -2.00000000000000     `, ` optnwrt:cpaitr_show  T`。名前と値の間の空白数は不定。
- 実効値: ヘッダの `meshr   mse    ng   mxl` の次の行に、option の有無によらず実際に使う `mse` と `ng` が出る（Cu go、`mse= 5` で `400     5    21     3`）。
- エラー: 未知の `key=` は `stop 100`（`AkaikkrJob.run` は return code ≠ 0 で `KKRFailedExecutionError("return_code=100")` を上げるが、`check_stopped_by_errtrp` は最終行 `***err in` しか見ないので理由は伝わらない）。`end_option` 前の EOF は `failed to read "end_option", but found EOF`。
- 値は 80 文字で切られる。論理値は list-directed 読み（`T`/`F`/`.true.`/`.false.`、ifort は `True` も受ける）。

## 2. ファイル構成

| 区分 | パス | 内容 |
|---|---|---|
| 追加 | `library/PyAkaiKKR/src/pyakaikkr/option.py` | `OptionKey`, `OPTION_KEYS`, `canonical_name`, `list_option_keys`, `normalize_option`, `format_option_card`, `parse_option_block`, `read_inputcard_option`, `parse_option_echo`（§3〜§5） |
| 修正 | `pyakaikkr/AkaiKkr.py` | `make_option_card` を `option.format_option_card` の呼び出しにする。`get_option`, `get_emesh_param`, `check_option_error` を追加（§5）。`run` の失敗時メッセージに `unknown token` を含める |
| 修正 | `pyakaikkr/Error.py` | `KKRUnknownOptionError(ValueError)`, `KKROptionValueError(ValueError)` |
| 修正 | `pyakaikkr/__init__.py` | `from .option import *`（`__all__` で公開名を絞る） |
| 修正 | `pyakaikkr/ase/calculator.py` | `option=` を `normalize_option` に通す。`get_option()` を追加（§6） |
| 修正 | `docs/akaikkr_option_keys.md` | §1 末尾に「Python からは `pyakaikkr.OPTION_KEYS` で引ける」と、echo の事実（§1 のこと）を追記 |
| 追加 | `tests/option/test_option.py`, `tests/option/data/out_go_option.log` | pytest と、上の Cu 実測出力（162 行）を fixture として保存 |
| 追加 | `docs/option_access_usage.md` | 使い方（実装後） |

## 3. キーの辞書（A）

```python
@dataclass(frozen=True)
class OptionKey:
    name: str                       # 正式名（Fortran の変数名）
    aliases: tuple[str, ...]        # specx が受ける別名
    kind: str                       # "int" | "float" | "str" | "bool" | "list"
    default: object | None          # specx 内部の既定値。ewidth 依存など固定で書けないものは None
    codes: tuple[str, ...]          # 効くコード: "akaikkr", "akaikkr_cnd", "cpa2021v01"
    echo: bool                      # 出力に optnwrt:<echo_name> が出るか
    echo_name: str                  # echo 行での名前。cemesr_ref だけ "ref"、他は name と同じ
    description: str                # 1 行

OPTION_KEYS: dict[str, OptionKey]   # 正式名 → OptionKey、16 個（15 スカラ + klabel）
```

- 内容は `docs/akaikkr_option_keys.md` §2 の表をそのまま写す。表と辞書が食い違わないよう、テストで表の行数とキー名を突き合わせる（表の `| \`name=\`` を正規表現で拾う）。
- `canonical_name(key: str) -> str`: 末尾の `=` は取り、別名なら正式名に。未知なら `KKRUnknownOptionError`。
- `list_option_keys(code: str | None = None) -> list[OptionKey]`。`code` を渡すとそのコードで効くものだけ。
- `default`: `mse` は None（`2·int(ewidth·35/2)+1`）、`tol` 1e-6、`ng` 21、`dex` 0.05、`rkick` 1、`critic` -1.0、`cemesr_ref` 0.75（akaikkr）/ 0.5（akaikkr_cnd。コード依存なので `default` は dict `{"akaikkr": 0.75, "akaikkr_cnd": 0.5}` を許す）、`spmain_bnd2` "a"、`spckkr_dmpc0` 1.0、`spckkr_itrmx` 100、`ie` 3、`cpaitr_show` False、`cpaitr_tol` 1e-8、`ddos` False、`tempmu` 1e-6、`klabel` None。

## 4. 検証と書き出し（D）

```python
def normalize_option(option: dict, *, code: str | None = None, strict: bool = True) -> dict
def format_option_card(option: dict) -> list[str]
```

- `normalize_option` は正式名に揃えた新しい dict を返す（入力は変更しない）。
  - 値は inputcard に書く文字列に変換する。`bool` → `"T"`/`"F"`（ifort/gfortran どちらの list-directed 読みでも確実な形）。`int`/`float` は `str()`（`1e-05` のような指数表記は Fortran が読める）。`list` は要素を `str()` にしたリスト。`str` はそのまま。
  - `kind` と合わない値（`mse` に `"abc"`、`spmain_bnd2` に `"z"` 以外の検査は最小限で、`int`/`float` は変換できるか、`bool` は真偽に解釈できるかだけ）は `KKROptionValueError`。
  - 未知キー: `strict=True` なら `KKRUnknownOptionError`、`False` なら `warnings.warn` して**そのまま残す**（specx 側が新しいキーを持つ場合の逃げ道）。
  - `code` を渡し、そのコードで効かないキー（akaikkr に `cpaitr_show` など）があれば `warnings.warn`（止めない。§1 のとおり specx は停止しない）。
  - 値が 80 文字を超えれば `KKROptionValueError`。
- `format_option_card` は現在の `make_option_card` と同じ行を返す（`begin_option` の前後の空行、先頭 1 空白、`key= value`、配列は `begin_key`/1 行/`end_key`）。`make_inputcard` は `dic["option"]` があれば `format_option_card(normalize_option(dic["option"], strict=True))` を書く。既存の inputcard との差分は、別名で渡していた場合に正式名になることだけ。
- `make_inputcard` の `dic["option"]` が空 dict のときはブロックを書かない（いまと同じ）。

## 5. 読み取り（B, C）

```python
def parse_option_block(lines: Iterable[str], *, typed: bool = True, strict: bool = True) -> dict
def read_inputcard_option(inputcard: str | Iterable[str], *, typed: bool = True) -> dict
def parse_option_echo(lines: Iterable[str], *, typed: bool = True) -> dict
```

- `parse_option_block`: specx と同じトークン規則で `begin_option` 〜 `end_option` を読む。`#` で始まる行はコメントとして捨てる（specx の `xtoken` と同じ）。
  - `key=` の次のトークンを値に。`typed=True` なら `OptionKey.kind` で `int`/`float`/`bool`/`str` に変換（`bool` は `T`/`F`/`.true.`/`.false.`/`true`/`false`/`True`/`False` を大小無視で受ける）。`typed=False` なら文字列のまま。
  - `begin_x` 〜 `end_x` は `list[str]`。`klabel` 以外の名前は `strict=True` なら `KKRUnknownOptionError`、`False` なら警告して dict に入れる（specx は無視するだけなので、既存 inputcard を読む用途では `False` が便利）。
  - `mse=5` のように `=` の後ろが続くトークン、`end_option` の無い入力は `KKRUnknownOptionError`。specx の停止と同じ条件で止める。
  - 戻りは正式名 → 値。別名は正式名に直す。
- `read_inputcard_option`: ファイルパスまたは行列を受け、`begin_option` トークンを探して `parse_option_block` に渡す。ブロックが無ければ `{}`。`spc` 入力では `atmicx` の後と k 点の前の 2 箇所に書けるので、複数あれば**後のブロックで上書き**した dict を返す（specx も 2 回目の `optnrd_rd` で同じ `optparam` に上書きする）。全ブロックが要るときは `read_inputcard_option_blocks(...) -> list[dict]`。
- `parse_option_echo`: 出力の ` optnwrt:<echo_name> <value>` 行を `re.compile(r"^\s*optnwrt:(\w+)\s+(\S+)")` で拾い、`echo_name`（`ref` → `cemesr_ref`）を正式名に直して 正式名 → 値。同じキーが複数回出れば最後を採る（`critic` は once なので 1 回、`cemesr_ref` は dos で 2 回）。`typed=True` なら `kind` で変換（Fortran の `T`/`F`、`1.000000000000000E-005` を受ける）。`debug bnd2_type a` の行は拾わない。
- `AkaikkrJob` への追加:
  - `get_option(outfile) -> dict`: `parse_option_echo(self._read(outfile))`。option を渡していない出力では `{}`（例外にしない。他の `get_*` と違い「無い」が正常だから）。
  - `get_emesh_param(outfile) -> dict`: `meshr mse ng mxl` の見出し行の次の行を読み `{"meshr": int, "mse": int, "ng": int, "mxl": int}`。無ければ `KKRValueAquisitionError`。option の有無によらず実効値が取れる。
  - `check_option_error(outfile) -> str | None`: `unknown token:` または `failed to read "end_option"` の行があればその行を返す。`run` は return code ≠ 0 のとき `check_option_error` を呼び、見つかればメッセージに付けて `KKRFailedExecutionError` を上げる（`return_code=100, unknown token: mse=5`）。
- `read_inputcard_option` は `AkaikkrJob.read_inputcard_option(inputcard)` としても呼べる（`self.path_dir` 相対）。

## 6. ASE calculator

- `AkaiKKR(option=...)` の値は `set` 時ではなく inputcard を書く直前に `normalize_option(option, code=<推定>, strict=True)` に通す。`code` は実行ファイル名から推定できないので `AkaiKKR(code="akaikkr")` 引数を追加（既定 `"akaikkr"`、警告の精度にしか使わない）。
- `calc.get_option()`: 直近の出力に対する `AkaikkrJob.get_option`。`calc.get_emesh_param()` も同様。
- `results` には入れない（ASE の `results` はエネルギー等の物理量に限る）。

## 7. 他プロジェクトとの接続（実装はここではしない）

- aiida-akaikkr: parser で `AkaikkrJob.get_option(out)` と `get_emesh_param(out)` を `results["option"]`, `results["emesh"]` に入れると、`contour_bottom` 等が `cemesr_ref` を provenance から取れる（いまは 0.75/0.5 をビルドで決め打ち）。dos 出力には `optnwrt:cemesr_ref` が出る。
- ewidth 自動調整スキーム（`docs/ewidth_tuning_scheme.md` §3.1 の `option`）は `normalize_option` を通すだけ。

## 8. テスト（`tests/option/`）

specx 不要:

- `test_keys.py`: `OPTION_KEYS` が 16 個、`canonical_name("number_emesh=") == "mse"`、`canonical_name("thresh_go") == "tol"`、未知で `KKRUnknownOptionError`、`docs/akaikkr_option_keys.md` の表のキー名集合と一致。
- `test_normalize.py`: `{"mse": 3, "cpaitr_show": True, "klabel": ["G", "X"]}` → `{"mse": "3", "cpaitr_show": "T", "klabel": ["G", "X"]}`。`"number_emesh"` が `"mse"` になる。`code="akaikkr"` で `cpaitr_show` に警告。`strict=False` で未知キーが警告付きで残る。80 文字超で例外。
- `test_format.py`: `make_inputcard` の出力が現行と同じ（`docs/akaikkr_option_keys.md` §1 の例と一致）。option 無し・空 dict でブロックが出ない。
- `test_parse.py`: §1 の例、実測に使った inputcard（`mse= 5`, `tol= 1e-5`, `critic= -2.0`, `cemesr_ref= 0.6`, `klabel`, `begin_foo`）を `strict=False` で読んで `{"mse": 5, "tol": 1e-5, "critic": -2.0, "cemesr_ref": 0.6, "klabel": ["G","X","W"], "foo": ["1","2"]}`、`strict=True` で `foo` により例外。`mse=5`（空白無し）で例外。`end_option` 無しで例外。`make_inputcard` → `read_inputcard_option` の往復で元の dict（正規化後）に戻る。2 ブロックの上書き。
- `test_echo.py`: fixture `out_go_option.log` に対し `get_option` が `{"tol": 1e-5, "mse": 5, "critic": -2.0}`（`cemesr_ref` は go では出ないので含まれない）、`get_emesh_param` が `{"meshr": 400, "mse": 5, "ng": 21, "mxl": 3}`。`tests/akaikkr/Cu/out_go.log`（option 無し）があれば `get_option` が `{}`。`check_option_error` が `unknown token: mse=5` の行を返す（fixture `out_go_badoption.log`、3 行）。

specx が必要（`AKAIKKR_PROGRAM_PATH` があるときだけ、`tests/ase/conftest.py` の `needs_specx` を流用）:

- `test_option_run.py`: Cu go に `option={"number_emesh": 5, "tol": 1e-5}` を渡して `get_option` が `{"mse": 5, "tol": 1e-5}`、`get_emesh_param()["mse"] == 5`。`option={"mse=5": ...}` のような壊れたキーは `normalize_option` で実行前に止まる。`strict=False` で `{"nosuchkey": 1}` を通すと `KKRFailedExecutionError` のメッセージに `unknown token: nosuchkey=` が入る。

## 9. 実装時の注意

- `make_option_card` の出力形式は変えない（先頭 1 空白、`key= value`）。既存ユーザの inputcard 差分は別名→正式名のみ。
- `_read` はファイル単位でキャッシュするので `get_option` と `get_emesh_param` を同じ outfile に続けて呼んでも 1 回しか読まない。
- echo 行は値が「使われた」時点で出るため、`get_option` の結果は「この実行で効いた option」と読める。渡したが効かなかったキー（go 実行の `cemesr_ref` など）は `read_inputcard_option(inputcard)` との差分で分かる。この差分を返す `AkaikkrJob.unused_option(inputcard, outfile) -> set[str]` を用意する（`klabel` は echo されないので除外する）。
- `parse_option_echo` の正規表現は `optnwrt:` の直後に名前が続く前提（`this="optnwrt:"` は 8 文字固定）。
- Python 側の `OPTION_KEYS` は AkaiKKR の版に依存する。冒頭に対象版（2022.0721）を定数 `OPTION_KEYS_VERSION` として持ち、docstring で書く。

## 10. 実装メモ（2026-09-25）

- 実装済み: `pyakaikkr/option.py`, `AkaikkrJob.get_option / get_emesh_param / check_option_error / read_inputcard_option / unused_option`, `run` のメッセージ, ASE の `code=` と `get_option() / get_emesh_param()`, `tests/option/`（specx 不要 14 件、specx 3 件）。
- 仕様との差: `spmain_bnd2` は echo される（§1 を訂正）。`cemesr_ref` の echo 名は `ref`（`OptionKey.echo_name` を追加）。
- akaikkr_cnd の dos で `cemesr_ref= 0.75` を入れると DOS の窓が [-1.0, 1.0] から [-1.5, 0.5]（ewidth=2）に移ることを `tests/option/test_option_run.py::test_cnd_dos_with_cemesr_ref` で確認。その出力を `tests/option/data/out_dos_cnd_ref.log` に置き、specx 無しでも echo の読み取りを検査する。
