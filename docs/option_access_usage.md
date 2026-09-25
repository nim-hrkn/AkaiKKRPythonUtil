# `begin_option` を pyakaikkr から使う

対象: pyakaikkr 2023.2.1（2026-09-25 実装）、AkaiKKR 2022.0721。キーの意味は [akaikkr_option_keys.md](akaikkr_option_keys.md)、設計は [option_access_spec.md](option_access_spec.md)。

## 書く

入力 dict の `"option"` に dict を入れる。`AkaikkrJob.make_inputcard` が書く前に検証する。

```python
from pyakaikkr import AkaikkrJob
job = AkaikkrJob("Cu")
dic = dict(job.default)
dic["option"] = {"number_emesh": 5, "tol": 1e-5, "cpaitr_show": True, "klabel": ["G", "X", "W"]}
job.make_inputcard(dic, "inputcard_go")
# begin_option
#  mse= 5          <- 別名 number_emesh は正式名 mse に
#  tol= 1e-05
#  cpaitr_show= T  <- bool は T/F に
#  begin_klabel
#  G X W
#  end_klabel
# end_option
```

- 未知キーは `KKRUnknownOptionError`、型に合わない値（`mse` に文字列、`mse=5` のようなキー）は `KKROptionValueError`。specx が `unknown token` で止まる前に分かる。
- 検証だけしたいときは `normalize_option(option, code="akaikkr_cnd")`。`code` を渡すと、そのビルドで効かないキー（akaikkr に `cpaitr_show` など）に警告が出る。`strict=False` で未知キーをそのまま通せる。
- `pyakaikkr.OPTION_KEYS["cemesr_ref"]` で正式名・別名・型・既定値（ビルド依存なら `default_for("akaikkr_cnd")`）・効くコード・説明が引ける。`list_option_keys("akaikkr_cnd")` で一覧。

## 読む（inputcard）

```python
from pyakaikkr import read_inputcard_option
read_inputcard_option("Cu/inputcard_go")            # {'mse': 5, 'tol': 1e-05, 'cpaitr_show': True, 'klabel': [...]}
job.read_inputcard_option("inputcard_go")           # 同じ（job.path_dir 相対）
read_inputcard_option(text, typed=False)            # 値を文字列のまま
read_inputcard_option(path, strict=False)           # 未知キー・未知配列を警告付きで残す（specx が無視する begin_foo など）
```

spc 入力のように 2 箇所に書いてあれば後のブロックで上書きした dict になる。全部要るなら `read_inputcard_option_blocks`。

## 読む（出力）

```python
job.get_option("out_go.log")        # {'mse': 5, 'tol': 1e-05}  specx が使った option だけ
job.get_emesh_param("out_go.log")   # {'meshr': 400, 'mse': 5, 'ng': 21, 'mxl': 3}  option の有無によらず実効値
job.unused_option("inputcard_go", "out_go.log")   # 与えたが使われなかったキー、例 {'cemesr_ref'}
job.check_option_error("out_go.log")              # 'unknown token: mse=5' など。無ければ None
```

- specx は与えた値を**使う時点**で ` optnwrt:<name> <value>` を 1 行出す。それを拾うので、`get_option` の結果は「この実行で効いた option」。go では `cemesr_ref` は出ない（dos/spc でだけ出る。echo 名は `ref` だが戻り値は `cemesr_ref`）。`klabel` は出ない。
- option を渡していない出力では `{}`（例外にならない）。
- `AkaikkrJob.run` は specx が option で止まったとき（return code 100）、`KKRFailedExecutionError("return_code=100, unknown token: nosuchkey=")` のように理由を付ける。

## ASE calculator

```python
from pyakaikkr.ase import AkaiKKR
calc = AkaiKKR(option={"number_emesh": 5}, code="akaikkr")   # code は警告の判定にだけ使う
atoms.calc = calc
atoms.get_potential_energy()
calc.get_option()        # {'mse': 5}
calc.get_emesh_param()   # {'meshr': ..., 'mse': 5, 'ng': 21, 'mxl': ...}
```

option の誤りは `write_input` で `ase.calculators.calculator.InputError` になる。

## akaikkr_cnd の dos で ref を変える例

akaikkr_cnd ビルドの dos は既定 `ref=0.5` で、窓は [E_F − 0.5·ewidth, E_F + 0.5·ewidth]。akaikkr と同じ [E_F − 0.75·ewidth, E_F + 0.25·ewidth] にするには `cemesr_ref= 0.75` を渡す。

```python
dic["go"] = "dos"; dic["ewidth"] = 2.0
dic["option"] = {"cemesr_ref": 0.75}
job.make_inputcard(dic, "inputcard_dos")
job.run(specx_cnd, "inputcard_dos", "out_dos.log")
job.get_option("out_dos.log")            # {'cemesr_ref': 0.75}
e, dos = job.get_dos_as_list("out_dos.log")   # e は -1.5 〜 0.5 Ry（ref 無しなら -1.0 〜 1.0）
```

`tests/option/test_option_run.py::test_cnd_dos_with_cemesr_ref` がこれを確かめる（`AKAIKKR_PROGRAM_PATH` と `tests/akaikkr_cnd/Cu/pot.dat` が要る）。

## テスト

```
cd tests/option && python -m pytest -q                      # specx 不要 14 件
AKAIKKR_PROGRAM_PATH=/path/to/AkaiKKRprogram.2022.0721.ifort python -m pytest -q   # + specx 3 件
```
