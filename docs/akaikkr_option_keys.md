# AkaiKKR の `begin_option` ブロックとキーの意味

作成日: 2026-09-24
対象: AkaiKKRprogram.2022.0721（akaikkr, akaikkr_cnd, akaikkr_cpa2021v01）
根拠: `akaikkr_common/source/m_optn.f`（読み取り）、`akaikkr/source/specx.f`, `spmain.f`, `cemesr.f`, `spckkr.f`、`akaikkr_cnd/source/specx.f`, `cpaitr.f`, `cndct.f`, `m_seebeck.f`（使用箇所）

## 1. 書き方と読まれる位置

inputcard の原子リスト（`natm` 個の `atmicx atmtyp` 行）の直後に置く。

```
#--- natm
1
#--- atmicx atmtyp
0.0a 0.0b 0.0c Cu
begin_option
 mse= 3
 cpaitr_show= True
 begin_klabel
  G X W K G L
 end_klabel
end_option
```

- 読むのは `m_optn.f` の `optnrd_rd`。`specx.f` が `readin` の直後に 1 回呼ぶ。`spc` モードでは `spmain.f` が k 点を読む直前にもう 1 回呼ぶので、そこにも置ける。
- 次のトークンが `begin_option` でなければ何もしない（ブロックは省略可能）。
- キーと値は**別トークン**にする。`mse= 3` は可、`mse=3` は「unknown token」で停止する。
- `begin_xxx ... end_xxx` は文字列配列として読む。処理されるのは `klabel` だけで、他の名前は表示して無視する。
- 未知のキーは `unknown token` を表示して停止する。`end_option` の前でファイルが終わると停止する。
- `fsm` の `fspin` は option ブロックの後に書く（`specx.f` の `optnrd_rd` が先に走るため）。

pyakaikkr から渡すときは入力 dict に `"option"` キーを入れる。`AkaikkrJob.make_inputcard` の `make_option_card` が上の形式に書き出す（値がリストなら `begin_<key> ... end_<key>`、それ以外は `<key>= <value>`）。ASE calculator では `AkaiKKR(option={...})` で同じ dict を渡す。

```python
param["option"] = {"mse": 3, "cpaitr_show": True, "klabel": ["G", "X", "W"]}
```

Python からの参照（2026-09-25 実装、[option_access_usage.md](option_access_usage.md)）: `pyakaikkr.OPTION_KEYS` にこの表の内容（正式名、別名、型、既定値、効くコード）があり、`AkaikkrJob.make_inputcard` は書く前に `normalize_option` で検証する（未知キーは実行前に `KKRUnknownOptionError`、別名は正式名に、論理値は `T`/`F` に）。既存 inputcard は `read_inputcard_option`、出力は `AkaikkrJob.get_option`（specx が値を使った時点で出す ` optnwrt:<name> <value>` 行。`cemesr_ref` は `ref` の名前で、dos/spc でしか出ない）と `get_emesh_param`（ヘッダの `meshr mse ng mxl`）で読める。

## 2. キーの一覧

値は文字列として保持され、使う側で `read(...,*)` により数値化される。「使用箇所」は値が実際に効く場所。

| キー（別名） | 意味 | 既定値 | 使用箇所 |
|---|---|---|---|
| `mse=`（`number_emesh=`） | 複素エネルギー経路の点数。既定は `ewidth` から `2·int(ewidth·35/2)+1`（上限 201）。`dos`, `spc`, `mcd` では 201 に固定されるが、このキーで上書きできる | ewidth 依存 | `specx.f` |
| `tol=`（`thresh_scf=`, `thresh_go=`） | SCF 収束判定。ポテンシャル残差 `cnvq` がこれを下回ると反復終了 | 1e-6 | `specx.f` → `spmain.f`（`if(cnvq .lt. tol)`） |
| `ng=`（`ndegree_cheb=`） | エネルギー積分の Chebyshev 展開の項数。上限 `ngmx=21` | 21 | `specx.f` → `gtchst`, `banden` |
| `dex=` | 初期反復でスピンを分裂させる対称性破りの場。Fermi 準位を `ef±dex` にずらし、`kick` の場の大きさにも使う | 0.05 | `specx.f` → `spmain.f` |
| `rkick=` | `kick` の反復回数。`magtyp=kick3` の数字と同じ役割 | 1 | `spmain.f` |
| `critic=` | 混合法の切り替え点。`log10(cnvq)` がこの値を下回ると Tchebyshev 混合から単純混合に切り替え、`pmix` を 10 倍（上限 0.1）にする | -1.0 | `spmain.f` |
| `cemesr_ref=` | 実軸エネルギー経路（`cemesr.f`）で Fermi 準位を窓のどの位置に置くかの比。`kef = (mse-1)·ref+1`。0.75 なら窓の 75% が EF より下 | 0.75 | `cemesr.f` |
| `spmain_bnd2=` | バンドエネルギーの端点補正式の選択。`a`: `Im(e·z)`、`b`: `ef·Im(z)`、`c`: `Re(e)·Im(z)`。全エネルギーに影響する | `a` | `spmain.f` |
| `spckkr_dmpc0=` | スペクトル関数計算での CPA 反復の初期減衰係数 | 1.0 | `spckkr.f` |
| `spckkr_itrmx=` | 同じく CPA 反復の最大回数 | 100 | `spckkr.f` |
| `ie=`（`ndirection_cnd=`） | 伝導度を計算する電流方向の index（1, 2, 3 = x, y, z） | 3 | akaikkr_cnd `specx.f` → `cndct.f`。akaikkr 本体では読むだけで未使用 |
| `cpaitr_show=` | CPA 反復の途中経過を表示する | False | akaikkr_cnd `cpaitr.f` |
| `cpaitr_tol=` | CPA 反復の収束しきい値 | 1e-8 | akaikkr_cnd `cpaitr.f` |
| `ddos=` | DOS のエネルギー微分を計算して `ddos.txt` に出す（Seebeck 係数用）。`true`/`.true.`/`True` 系の文字列で真 | false | akaikkr_cnd `spmain.f` |
| `tempmu=` | Seebeck 係数計算（`go=sbk`）の温度。Fermi 分布の微分の幅を決め、`ewidth` もこれから自動決定される | 1e-6 | akaikkr_cnd `specx.f`, `m_seebeck.f` |
| `begin_klabel ... end_klabel` | k 点ラベルの配列。保持と表示だけで、この版の Fortran では計算に使う箇所はない | なし | `m_optn.f` のみ |

## 3. 注意

- `ie`, `cpaitr_show`, `cpaitr_tol`, `ddos`, `tempmu` は akaikkr_cnd 向け。akaikkr 本体に書いても停止はしないが効かない（`ie` は読まれて捨てられる）。
- 論理値は Fortran の list-directed 読みなので `T`/`F`、`.true.`/`.false.` が確実。デモ `demo/akaikkr_cnd.small.devel/cnd_test/finite_temperature.py` は `cpaitr_show= True` を渡しており、ifort ではこれも受け付ける。
- `tests/` の回帰テストは option ブロックを使っていない（すべて既定値）。`testrun_class.py` の `Cu_go`, `Cu_spc` にコメントアウトされた例が残っている。
