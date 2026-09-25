# ewidth 自動調整スキーム: GAES（Gap-Anchored Ewidth Search）と GAES-Committee の pyakaikkr 移植仕様

作成日: 2026-09-25（旧 `hea_scheme2_spec.md` を、HEA に限らない一般のアルゴリズムとして書き直し、名前を GAES に定めたもの）
使われた論文: T. Fukushima, H. Akai, T. Chikyow, H. Kino, Phys. Rev. Materials 6, 023802 (2022), doi:10.1103/PhysRevMaterials.6.023802（手法自体は未発表）
移植元: `fukushima_HEA_run_exprlattice/production_run/run0/`（2019-11、`1306.run_scheme2.py`, `HEARun.py`, `HEAPathSearch.py`, `hea_util.py`, `kkrinput_brvtyp.py`, `akaikkrio2.py`）
移植先: AkaiKKRPythonUtil（pyakaikkr 2023.2.1）、AkaiKKR 2022.0721（`akaikkr/specx`）

## 0. 背景と目的

全電子計算の電子状態には valence、semicore、core の領域がある。AkaiKKR の go は [E_F − ewidth, E_F] の複素エネルギー経路で電荷を積分するので、`ewidth` は valence 帯と、その下の core または semicore 帯との間の、DOS が小さいエネルギー領域（以下**バンドギャップ**と呼ぶ）の中を指すように選ばなければならない。ewidth が valence 帯の中にあれば電荷を取りこぼし、semicore 帯の中にあれば semicore の一部だけを valence として数えるので、どちらも自己無撞着な解にならない。

- バンドギャップは「連続した energy mesh の区間で DOS(E) < threshold」と定義する。threshold の既定は 1e-3（states/Ry、後述の単位）で、変更できる。DOS = 0 で定義しないのは、AkaiKKR は小さな `edelt`（虚部）を入れて計算するため、原理的に DOS = 0 のエネルギー領域が存在しないからである（VASP 等の実軸の計算では存在する）。
- 手順は「ewidth を仮に決めて go を回し、dos で DOS を出し、E_F − ewidth_go がバンドギャップの中にあるか見る。無ければ ewidth を直して go からやり直す」である。これまでは人間がこれを行っていた。AkaiKKR 本体にも、pyakaikkr にも、これを検知するアルゴリズムは入っていない。
- 2019 年の `fukushima_HEA_run_exprlattice/` は、この ewidth 調整を自動化したアルゴリズムで HEA（4 元等比、bcc/fcc、単一サイト CPA）の網羅計算を行ったものである。アルゴリズム自体は HEA に限らない一般のものなので、本仕様では HEA 固有の部分（2 桁原子番号の連結キー、等比濃度、`a=1000000`）を切り離し、任意の inputcard 辞書に対する **ewidth tuning scheme** として `pyakaikkr` に実装する。単一サイト CPA の入力は、等比に限らず任意の比の組成（`Rh0.5Pt0.5`、`B0.975Vc0.025` のような AkaiKKR の type 名と同じ書き方）を受ける `SiteComposition` から作る（§3.1）。
- `run0/RUN/` にはその結果が置かれているが、残っているのは**最終的に採用された ewidth の計算だけ**であり、途中の試行は無い（§2 の誤り #1 のため、ディレクトリ名の `ew_000` は試行番号として信用できない）。

### 0.1 名前

- **GAES（Gap-Anchored Ewidth Search）**: DOS のバンドギャップ区間に E_F − ewidth_go を固定（anchor）するまで ewidth_go を置き直して go / dos を繰り返す探索。族の名前であり、同時に移植元の逐次方式（ewidth を 1 本置き、DOS を見て直す、を繰り返す。日本語ではギャップ追従法）を指す。この手法は T. Fukushima, H. Akai, T. Chikyow, H. Kino, *Phys. Rev. Materials* **6**, 023802 (2022), [doi:10.1103/PhysRevMaterials.6.023802](https://doi.org/10.1103/PhysRevMaterials.6.023802) の HEA 網羅計算で使われたものだが、**手法そのものは未発表**（同論文には記述が無い）。移植元 `fukushima_HEA_run_exprlattice/production_run/run0/` はその計算のスクリプトであり、本仕様書が手法の最初の記述になる。
- **GAES-Committee**: ewidth_go を投機的に複数置いて並列に go + dos を走らせ（ewidth_dos は固定）、各 DOS から得たギャップ区間の投票（vote）でギャップを推定し、その推定から GAES で仕上げる方式（§7.1）。複数の ewidth_go の計算を committee（委員会）、集約段を vote と呼ぶ。
- パッケージ名 `pyakaikkr.gaes`、CLI 名 `kkr-gaes`、例外 `GaesError`。文書はこのファイル（`ewidth_tuning_scheme.md`）にまとめる。

### 0.2 判定は DOS の値で行い、微分は使わない

バンドギャップの端を DOS の微分（勾配の急変）で見つける方式は採らない。AkaiKKR の DOS は mesh 点ごとにギザギザしていることがあり（CPA の共鳴、粗い k 点、小さな edelt）、微分はその凹凸に敏感で一般的なアルゴリズムにならない。判定は §4 のとおり「連続した mesh 区間で DOS(E) < threshold」という値の比較だけで行う。平滑化してから微分する案も、平滑化幅という新たなパラメータを増やすだけなので採らない。

目的は次の 4 つ。

1. DOS（将来は PDOS）からバンドギャップ区間を検出し、与えた ewidth_go がその中にあるか判定し、無ければ次の ewidth 候補を返す関数群（§4、§5）。人間が手で行っていた判定をそのまま関数にしたもので、単独でも使える。
2. 上の判定と、edelt・pmix を段階的に落として SCF を追い込む手順を組み合わせた状態機械 GAES（§7）。移植元の scheme2。
3. 投機的並列実行と投票で初期 ewidth を決める GAES-Committee（§7.1）。
4. 2019 年の RUN を読み直す互換モード（§8）。

対象は akaikkr/specx。akaikkr_cnd は `displc` 必須で、dos の窓が ref=0.5 なので、§3.3 の ref を差し替えれば動くが本仕様の検証対象には入れない。

### 0.3 将来の拡張（PDOS）— 設計で先に確保しておくこと

いまは total DOS だけで判定する。HEA のように各元素が等しい濃度で入っていれば、どの元素の semicore も total DOS に同じ重みで現れるので、total DOS で正しいバンドギャップを検知できた。しかし希薄な成分では成り立たない。例えば La0.999Ge0.001 では Ge の total DOS への寄与は濃度分（0.001 程度）しかなく、threshold 1e-3 と同程度なので、total DOS からは Ge の semicore を認識できない。この場合は Ge の **PDOS（成分ごとの DOS）**で semicore を検知しなければならない。

- したがって §4 の区間検出は「energy mesh と、任意本数の DOS 曲線の列」を入力に取り、全曲線の AND でバンドギャップを決める形にする。曲線が polytyp ごとの total DOS でも、成分ごとの PDOS でも同じ関数が使える。
- AkaiKKR の `DOS of component i` は成分 i の on-site Green 関数から出る成分あたりの DOS（`spmain.f`: `-Im wkc(l,i,k)/π`）で、濃度は掛かっていないと読める（実装時に希薄系で確認する）。掛かっていなければ、希薄成分の semicore は PDOS では本来の高さで見え、threshold をそのまま使える。
- PDOS 用の threshold は total DOS 用と別に持てるようにする（`dosth_pdos`）。PDOS は l 成分に分かれて出るので、成分ごとに l を足した値を曲線にする（`get_pdos_as_list` の `output_format="spin_separation"` の l 和）。

## 1. 用語と記号

| 記号 | 意味 |
|---|---|
| ewidth_go | go（と j, tc）の ewidth。SCF の積分路は [E_F − ewidth_go, E_F] |
| ewidth_dos | dos の ewidth。窓は [E_F − ref·ewidth_dos, E_F + (1 − ref)·ewidth_dos]、ref は akaikkr で 0.75、akaikkr_cnd で 0.5、`begin_option` の `cemesr_ref=` で変更（[akaikkr_option_keys.md](akaikkr_option_keys.md)、[option_access_usage.md](option_access_usage.md)） |
| E | dos 出力の横軸 E − E_F（Ry） |
| dosth | バンドギャップ判定の threshold。既定 1e-3、変更可 |
| バンドギャップ区間 | mesh 上で連続して DOS(E) < dosth となる区間 [e1, e2] |
| eth | 区間を採用する最小幅（Ry）。移植元 0.30 |
| ediff | 区間上端からの余裕（Ry）。移植元 0.20。新 ewidth は −(e2 − ediff − margin)、margin 0.01 |

DOS の単位: pyakaikkr の `get_dos` が返す total DOS はスピンごと（states/Ry/cell/spin）。判定は**スピン和**（up + dn、nmag なら up の 1 本）に対して行う。移植元は平均 (up+dn)/2 に 1e-2 を掛けていた（= スピン和に 2e-2）。ユーザー指定の既定 1e-3 はスピン和に対する値とし、互換モードは 2e-2 を使う。

## 2. 移植元の誤りと修正方針

| # | 箇所 | 内容 | 修正 |
|---|---|---|---|
| 1 | `1306.run_scheme2.py` STEP1 | 新しい ewidth（`nextplan=="new"`）が返っても `iew` を進めないため、同じディレクトリ `ew_000` の既存 out_go.log を再利用して再計算しない。str.out で `('new', ...)` は 2053 回出ているが、そのすべてが ewidth=1.2 の結果を「新 ewidth の結果」として再判定している（2 回目は必ず `old` になる構造）。 | 実行ディレクトリ名に実際の ewidth 値を含め、inputcard の ewidth と一致しなければ再計算する（§6）。互換モードでも既存結果の再利用は inputcard の内容が一致する場合に限る。 |
| 2 | `HEARun.analyzeDOS.__init__` | 低 DOS 領域が上端まで続くとき、未定義名 `smalldos` を参照して NameError。 | `len(mask)` を使う（§4）。 |
| 3 | `akaikkrio2.OutputGo.converging` | `falg = True` のタイポ。当該分岐（`r2moment>0.8 and r2err>0.8`）は直前の `r2err>0.8` に含まれるので削除する。 | §5 の判定式に置き換える。 |
| 4 | `akaikkrio2.OutputDOS.get` | スピン和ではなく平均 `(up+dn)/2` を返しているので、`dosth=1e-2` は平均 DOS への閾値になっている。 | スピン和で判定する（§1）。 |
| 5 | `HEArun.make_heainput` | `self.dic` を破壊的に更新するため、go → dos → j と呼ぶうちに `go`, `ewidth`, `record` が前の呼び出しの値を引きずる。 | 入力辞書は毎回 `deepcopy` して作る。 |
| 6 | `HEArun.run_all` | go が未収束でも dos と j を実行する（`stopifnotconverged=False` 固定）。 | 既定では収束した go の後だけ dos / j を実行する。互換モードでは移植元どおり常に実行する。 |
| 7 | `1306.run_scheme2.py` | `nextplan=="fail"`（新 ewidth を作れない）を `keyconverged=True` にして「収束」と同じ扱いで終了する。 | 状態を `ewidth_fail` として記録し、収束とは区別する（§7）。 |
| 8 | `calculate_pm` | `ipm` をローカルで進めるので、呼び出し元の `ipm` と実際に使ったディレクトリ番号がずれる。 | 使った pm 段の番号を戻り値で返す。 |
| 9 | `analyzeDOS`, `make_new_ewidth` | 裸の `raise`（RuntimeError）。 | 専用例外 `GaesError` を送出する。 |
| 10 | `HEADos` | 未定義の `OutFirstDOS` を使う。 | 削除し、プロットは `pyakaikkr.DosPlotter`（go の ewidth の線付き、[dos_plot_ewidth_line.md](dos_plot_ewidth_line.md)）に任せる。 |
| 11 | `hea_util` | `from pymatgen import Element` は pymatgen 2022 以降で動かない。 | `pymatgen.core.periodic_table.Element` を使う。 |
| 12 | `OutputGo.__init__` | 履歴が空だと `h_err[-1]` で IndexError。 | `AkaikkrJob.get_convergence`（`*** no convergence` の有無）で判定する。specx は `tol=1e-6`（log10 err ≤ -6）で収束、`maxitr` 到達で `*** no convergence` を出すので、移植元の `err <= -6` 判定と等価。 |

移植元で仕様どおりだったもの（誤りではない）:

- inputcard の `a=1000000` は標準 specx の機能で、`chklat.f` が `a > 0.99e6` のとき `qvolum(anclr,1)`（純物質の実験原子体積、cc/mol）を濃度平均して格子定数を決める。`a ≈ 0` なら MJW 値の表（`qvolum(anclr,2)`）を使う。合金化による体積変化は入らない。
- `record=2nd` は pm_000 でも指定されるが、pot.dat が無いときは specx が `eof detected; data generated` と警告して初期化から始めるので害はない。

## 3. ファイル構成

| 区分 | パス | 内容 |
|---|---|---|
| 追加 | `library/PyAkaiKKR/src/pyakaikkr/gaes/__init__.py` | 公開 API（`gap_regions`, `check_ewidth`, `choose_ewidth`, `is_converging`, `KkrRunner`, `Layout`, `Gaes`, `GaesCommittee`, `SiteComposition`） |
| 追加 | `pyakaikkr/gaes/gap.py` | バンドギャップ区間の検出（§4） |
| 追加 | `pyakaikkr/gaes/ewidth.py` | ewidth の判定と候補生成（§4.2） |
| 追加 | `pyakaikkr/gaes/convergence.py` | SCF 履歴の「まだ落ちている」判定（§5） |
| 追加 | `pyakaikkr/gaes/runner.py` | 1 ディレクトリで go / dos / j を実行する `KkrRunner`（§6） |
| 追加 | `pyakaikkr/gaes/layout.py` | ディレクトリ命名と検索（§6.1） |
| 追加 | `pyakaikkr/gaes/scheme.py` | 状態機械 `Gaes` と結果集計（§7, §9） |
| 追加 | `pyakaikkr/gaes/committee.py` | `GaesCommittee`: 投機的並列実行と vote（§7.1） |
| 追加 | `pyakaikkr/gaes/cli.py` | `kkr-gaes` コマンド（§10） |
| 追加 | `pyakaikkr/gaes/composition.py` | `SiteComposition`（任意の元素と比の単一サイト組成、AkaiKKR の type 名との相互変換、type 名の長さ検査）、`make_single_site_param`（§3.1） |
| 追加 | `pyakaikkr/gaes/legacy_hea.py` | 2019 年の 2 桁原子番号キー（`13142122`）⇄ `SiteComposition`、`heakeylist0.csv` の読み込み。互換モード専用 |
| 修正 | `pyakaikkr/Error.py` | `GaesError(Exception)` を追加 |
| 修正 | `library/PyAkaiKKR/setup.cfg` | `console_scripts: kkr-gaes = pyakaikkr.gaes.cli:main` |
| 追加 | `tests/gaes/` | pytest（§11） |
| 追加 | `docs/gaes_usage.md` | 使い方（実装後） |

`pyakaikkr.gaes` は `import pyakaikkr` からは読み込まない（明示 import）。

### 3.1 単一サイト組成と入力（`composition.py`, `legacy_hea.py`）

```python
@dataclass(frozen=True)
class SiteComposition:
    elements: tuple[str, ...]        # 元素記号（"Vc" = 空孔、Z=0 も可）、入力順を保つ
    fractions: tuple[float, ...]     # 合計 1（1e-6 以内）。等比なら 1/n
    @classmethod
    def from_elements(cls, elements: Sequence[str] | str) -> "SiteComposition"   # 等比。["Al","Si","Sc","Ti"] または "AlSiScTi"
    @classmethod
    def from_dict(cls, comp: Mapping[str, float]) -> "SiteComposition"           # {"Rh": 0.5, "Pt": 0.5}
    @classmethod
    def from_type_name(cls, name: str) -> "SiteComposition"                      # "Rh0.5Pt0.5"（末尾の "_1" は捨てる）、"Fe"、"B0.975Vc0.025"
    @property
    def z(self) -> tuple[int, ...]                                               # 原子番号（Vc は 0）
    @property
    def conc(self) -> tuple[float, ...]                                          # 百分率（AkaiKKR の conc）
    def type_name(self, suffix: str = "", ndigits: int = 3, max_len: int = 40) -> str
    def key(self) -> str                                                         # ディレクトリ名に使う安全な文字列（type_name から "." を "p" に）
```

- 組成は**任意の元素数と比**。HEA の 4 元等比は `from_elements("AlSiScTi")` の特殊な場合にすぎない。テストセットの `Rh0.5Pt0.5_1`（FeRh05Pt05）や `B0.975Vc0.025_1`（FeB195）と同じ書き方を `from_type_name` で読み、`type_name()` で書く。
- `type_name` は pyakaikkr の Cif2Kkr と同じ書式（元素記号 + 分率、分率 1 なら省略、`_<site>` の接尾辞は `suffix` で付ける）。`ndigits` は分率の桁数。**長さが `max_len`（既定 40）を超えれば `ValueError`**（§3.1.1）。空白とカンマを含めば `ValueError`。
- 2019 年のキーは `legacy_hea.py` の `heakey_to_composition("13142122") -> SiteComposition`（2 桁ずつ原子番号、等比。桁数が奇数なら `ValueError`）と `composition_to_heakey(comp)`（等比でなければ `ValueError`）、`load_heakeylist(path)`（列 `heakey, elements`、`heakey` は文字列で読む。先頭ゼロ保護）で扱う。

```python
def make_single_site_param(comp: SiteComposition, brvtyp: str, *, lattice: str | float = "expr",
                           go="go", ewidth=1.2, edelt=1e-4, pmix=0.005, maxitr=500, bzqlty=10,
                           sdftyp="pbe", reltyp="sra", magtyp="mag", record="2nd", outtyp="update",
                           rmt=1.0, mxl=2, field=0.0, type_name: str | None = None,
                           potentialfile="pot.dat", option: dict | None = None) -> dict
```

- `AkaikkrJob.default` を `deepcopy` して更新した辞書を返す。`ntyp=1, type=[type_name], ncmp=[len(comp.z)], anclr=[list(comp.z)], conc=[comp.conc]`, `natm=1, atmicx=[["0.0a","0.0b","0.0c",type_name]]`。`type_name` 省略時は `comp.type_name()`。移植元と同じ inputcard を出したいとき（互換モード）は `type_name="HEA"` と整数の conc（等比 4 元で `25`）を使う。
- `lattice`: `"expr"`（既定、`a=1000000`、実験原子体積の濃度平均）、`"mjw"`（`a=0`）、`float`（bohr）。
- `option` は `normalize_option` を通して `dic["option"]` に入れる。

#### 3.1.1 AkaiKKR の type 名・atmtyp の制限（ソース確認、2026-09-25）

`AkaiKKRprogram.2022.0721` のソースと、Cu で 40 / 41 / 45 文字の type 名を実際に走らせて確かめた。

| 項目 | 値 | 根拠 |
|---|---|---|
| type 名と atmtyp の最大長 | **40 文字** | `akaikkr_common/source/m_param.f`: `len_type = 40`。`akaikkr/source/specx.f` と `akaikkr_cnd/source/specx.f` は `type(:)*(len_type), atmtyp(:)*(len_type)`。`akaikkr_cpa2021v01/source/specx.f` は `*40` の直書き。3 ビルドとも 40 |
| 超えたとき | **エラーにならず、41 文字目以降を黙って捨てる** | `readin.f` の `type(i)=token`, `atmtyp(i)=token` は Fortran の文字代入なので切り詰め。type と atmtyp が同じ長さで切られるので `ty2ity.f` の照合は通る。出力（`type=` 行、`*** type-... ***` 行、反復開始行）には 40 文字だけ出る。41 文字と 45 文字で確認 |
| 落とし穴 | 40 文字目までが同じ 2 つの type は同一とみなされ、`ty2ity` は先に定義した方に結び付ける | 上と同じ理由。pyakaikkr 側で検査して弾く（`SiteComposition.type_name` と、将来 `make_inputcard` にも） |
| 使える文字 | 空白とカンマ以外 | トークンの区切りは空白または単一のカンマ（`akaikkr_common/source/xtoken.f`）。`.`、`_`、数字は可（`Rh0.5Pt0.5_1`, `B0.975Vc0.025_1` が実例） |
| 1 行の長さ | 1024 文字 | `m_param.f`: `len_token = 1024`（`xtoken.f` の行バッファ） |
| 反復開始行の表示 | 型名の後ろの `_2` 等は同じ type の原子数 | `spmain.f` 707 行: `write(cwtyp,'(a,i3)')'_',iwtyp(i)`。FeB195 の `B0.975Vc0.025_1_2` は名前ではなく「2 原子」の意味。title は 2000 文字 |

pyakaikkr の `get_type_of_site` と `get_component_moment` は空白区切りで型名を取るので 40 文字までなら影響しない。ASE calculator の `check_kkr_output_structure` は入力の型名と出力の型名を比べるので、41 文字以上を渡すと不一致になる。Cif2Kkr が付ける名前（元素 + 分率 + `_<site>`）は 8 元・3 桁分率で 40 文字を超え得る（例 `Al0.125Si0.125Sc0.125Ti0.125V0.125Cr0.125Mn0.125Fe0.125_1` は 57 文字）ので、`SiteComposition.type_name(ndigits=2)` や短縮名を使う。

scheme 側は「inputcard 辞書（`AkaikkrJob.make_inputcard` に渡せるもの）」だけを受けるので、HEA 以外は CIF（`AkaikkrJob.read_structure`）や ASE（`pyakaikkr.ase.atoms_to_kkr_param`）から作った辞書をそのまま渡せる。

### 3.2 ewidth_dos と dos の窓

- go と j は同じ ewidth_go を使う。dos の edelt は go の edelt に連動させず `edelt_dos`（既定 1e-4 = STEP1 の edelt）に固定する。**edelt を 1e-2 のように大きくするのは SCF のためであり、バンドギャップ領域の認識には適さない**（大きい edelt はギャップ領域の DOS を持ち上げる）。STEP2 で edelt を 1e-2 から始めるようにしたため（§7）、go の edelt を dos に使うと DOS が広がってギャップ区間が縮み、判定が段ごとに変わってしまう（実測: AlGeHfBi fcc で edelt 1e-2 の dos はギャップが (−2.10, −2.01) しか残らず即 `fail`）。移植元は go と同じ edelt を dos にも使っていた（`compat=True` で再現）。dos は ewidth_dos（既定 3.0、移植元と同じ）で、窓は [E_F − 0.75·3.0, E_F + 0.25·3.0] = [−2.25, 0.75] Ry、201 点、刻み 0.015 Ry。§4 の判定はこの mesh 上で行う。
- 窓の下端が −ewidth_go − eth より上にあると、ewidth_go の下側に幅 eth の区間があっても検出できない。`Gaes` は各 dos の前に `ref·ewidth_dos ≥ ewidth_go + eth + ediff` を確かめ、満たさなければ ewidth_dos を `(ewidth_go + eth + ediff)/ref` に広げる（`ewidth_dos_auto=True`、既定）。ref は `dosth` と同じくパラメータ（既定 0.75。akaikkr_cnd なら 0.5、`option={"cemesr_ref": ...}` を渡すならその値）。
- mesh 点数は固定なので ewidth_dos を広げると刻みが粗くなる。eth より刻みが十分小さいこと（刻み ≤ eth/10）を確かめ、満たさなければ警告する。
- **dos の窓は 4.5 Ry まで**（`ewidth_dos_max=4.5`）。AlGaSnPb fcc で窓を 6.0 Ry（下端 −4.5 Ry）にすると、core 準位（−6 Ry より深い）も無い −2〜−4.5 Ry に 10〜80 states/Ry の偽の構造が出た（エネルギー補間の適用範囲外）。また窓が深い core 準位を跨ぐと akaikkr ビルドは `reconf` で止まる。窓の下端に接する低 DOS 区間は候補にせず（`Decision.window_limited`）、それが理由で候補が無いときだけ Gaes が窓を 1.5 倍ずつ上限まで広げて dos を取り直す。広げた dos が specx で止まれば元の窓に戻して判定する。
- **`min_ewidth=1.0`, `max_ewidth=2.0`**（既定、2026-09-25 ユーザー指定）: 候補 ewidth の下限と上限（Ry）。Ga 3d / Pb 5d / Sn 4d のように semicore が並んで幅 eth のギャップが無い系では、候補が semicore の下の深い区間へ潜り込み、semicore をまとめて valence に切り替えてしまう。上限を与えるとそれを止めて `fail` にできる。

## 4. バンドギャップ区間の検出と ewidth の判定（`gap.py`, `ewidth.py`）

判定は DOS の値の threshold 比較だけで行う（§0.2。微分は使わない）。

### 4.1 区間検出

```python
@dataclass(frozen=True)
class GapRegion:
    e1: float; e2: float                 # Ry, E - EF
    i1: int; i2: int                     # mesh index（両端を含む）
    @property
    def width(self) -> float

def gap_regions(energy: np.ndarray, curves: Sequence[np.ndarray], dosth: float = 1e-3,
                *, labels: Sequence[str] | None = None) -> list[GapRegion]
def dos_curves_from_outputs(jobs_and_files: Sequence[tuple[AkaikkrJob, str]], *,
                            spin_sum: bool = True) -> tuple[np.ndarray, list[np.ndarray], list[str]]
def pdos_curves_from_output(job: AkaikkrJob, outfile: str, *, l_sum: bool = True,
                            spin_sum: bool = True) -> tuple[np.ndarray, list[np.ndarray], list[str]]
```

- `energy` は E − E_F（Ry）の mesh、`curves` は同じ mesh 上の DOS 曲線の列（polytyp ごとの total DOS、成分ごとの PDOS、あるいはその両方を連結したもの）。全曲線の mesh が一致しなければ `GaesError`。
- 各 mesh 点で「すべての曲線が dosth 未満」を取り（AND）、連続して真である mesh 区間を `GapRegion` にする。区間の上端 e2 は**最後に条件を満たした点の次の点**（移植元と同じ。`i2` はその index）。末尾まで続く場合は最後の点（移植元の NameError 箇所）。
- `dos_curves_from_outputs` は `get_dos_as_list` のスピン和を曲線にする。`pdos_curves_from_output` は `get_pdos_as_list(output_format="spin_separation")` を成分ごとに l 和・スピン和して曲線にする（§0.3 の PDOS 拡張の入口。いまは用意だけして `Gaes` からは使わない）。

### 4.2 ewidth の判定と候補

```python
def check_ewidth(regions: Sequence[GapRegion], ewidth: float, *, eth: float = 0.30,
                 ediff: float = 0.20) -> GapRegion | None
def choose_ewidth(regions: Sequence[GapRegion], ewidth: float, *, eth: float = 0.30,
                  ediff: float = 0.20, margin: float = 0.01) -> tuple[str, float | None, list[float]]
```

- `check_ewidth`: 区間を高エネルギー側から見て、幅 `e2 − e1 > eth` かつ `e1 < −ewidth < e2` かつ `−ewidth < e2 − ediff` を満たす最初の区間を返す（無ければ None）。「E_F − ewidth_go がバンドギャップの中にあり、上端から ediff 以上離れている」の判定そのもの。
- `choose_ewidth` は移植元の `calc_new_ewidth`:
  1. `check_ewidth` が区間を返せば `("old", ewidth, candidates)`。
  2. 無ければ幅 `> eth` の各区間から `−(e2 − ediff − margin)` を候補にして高エネルギー側から並べ、`("new", candidates[0], candidates)`。候補が無ければ `("fail", None, [])`。
- 候補は「区間の上端の少し下」なので、valence 帯の直下のギャップがまず選ばれる。semicore を valence に含めたいときは `min_ewidth` で下限を与えて候補を絞る（新規、既定 None）。

### 4.3 ewidth 判定の 2 方式: Method 1 と Method 2

§4.2 の判定（2019 年のスクリプトと同じ、threshold 1 つ）を **Method 1** と呼ぶ。2026-09-25 の Bi 系のテスト（§13.2）で、Method 1 は threshold 2e-2 では semicore ピークの裾（DOS 3〜5e-3）に ewidth を置いても `old` にし、threshold 1e-3 では区間が eth より狭くなって `fail` にすることが分かった。これを直すのが **Method 2**（threshold 2 段）で、`Gaes(method=2)` で選ぶ（既定は 2。`compat=True` は 1）。

**Method 2 の考え方**: 粗い threshold `dosth` でバンドギャップ「らしい」広い区間を認め、その中で細かい threshold `dosth2` を満たす（DOS が本当に小さい）部分に、valence 帯の底から ediff の余裕を取って ewidth を置く。

```python
def choose_ewidth2(energy, curves, ewidth, *, dosth=2e-2, dosth2=1e-3, eth=0.30, ediff=0.20,
                   margin=0.01, dosth2_relax=2.0, min_ewidth=1.0, max_ewidth=2.0) -> Decision
```

手順:

1. **粗い区間**: `gap_regions(energy, curves, dosth)` で区間を求め、幅 > eth のものだけ残す（Method 1 の条件 1 と同じ）。mesh の下端に接する区間（e1 = mesh の最小値）は候補の順位を最後にする。その区間の上端は実際の帯の縁なので候補にも `old` 判定にも使うが、位置は「積分路の下端の下に幅 eth 以上の低 DOS が確認できる」範囲（−ewidth ≥ e_min + eth）に限る。それで使える位置が無いときだけ `window_limited` を立て、呼び出し側が窓を広げる（2026-09-25 修正。それ以前は下端に接する区間を丸ごと捨てていたため、Pm-Mn-Fe-Co などで 5p 帯の下のギャップが候補にならなかった）。
2. **細かい部分区間**: 各粗い区間 [e1, e2] の内側で `gap_regions(energy, curves, dosth2)` を求める（粗い区間の内側にある区間だけ。幅の下限は課さない）。
3. **上端の余裕**: 粗い区間の上端（valence 帯の底）から `top = e2 − ediff` を取る。−ewidth は `top` より下でなければならない（Method 1 の条件 3 と同じ）。
4. **判定 `old`**: いまの −ewidth が、いずれかの粗い区間の細かい部分区間の中にあり、かつ −ewidth < top なら `old`。
5. **候補**: 各粗い区間について、各細かい部分区間 [f1, f2] に対して `hi = min(f2, top)` を取り、`hi − margin > f1` なら候補 `−(hi − margin)`。つまり「DOS < dosth2 の範囲の中で、valence 帯の底から ediff 以上離れた最も浅い位置」に積分路の下端を置く。細かい部分区間が複数あれば上（高エネルギー側）のものを先に、粗い区間は高エネルギー側から並べる（mesh 下端に接するものは最後）。候補が `min_ewidth`〜`max_ewidth` の外にあるときは、**捨てずに範囲の端へ寄せ**、寄せた位置がまだ細かい部分区間の中（かつ top より下）にあればそれを候補にする（2026-09-25 修正。それ以前は候補 1 点が範囲外なら区間ごと捨てていたため、Pr/Nd/Pm-Mn-Fe-Co の 5p 帯の上のギャップ [−1.18, −0.77]（候補 0.78）が min_ewidth=1.0 で消え、1.0〜1.17 という範囲内の位置を試さなかった）。Method 1 の候補も同じ扱い。`fail` のときは、区間ごとに候補が出なかった理由を `Decision.reasons`（`Judgement.reasons`、`KeyResult.message`、aiida の history `reasons`、`kkr-gaes check` の `note:`）に残す: 幅 eth 未満、DOS < dosth2 の部分区間無し（最小 DOS を併記）、ediff で除外（ediff がいくつ未満なら候補がいくつになるかを併記。例: In-Mn-Fe-Co の In 4d core、[−1.10, −0.94] は「ediff < 0.17 なら候補 ≈ 1.09」）、範囲外、窓の下端。**ギャップ区間かどうかの判定（手順 1〜4）は `min_ewidth` / `max_ewidth` に依らず dos の窓全体で行う**。範囲 [min_ewidth, max_ewidth] は ewidth を選ぶ（手順 5）ときだけに使う（2026-09-25 ユーザー指定）。したがって範囲の外にあるギャップ区間も判定結果（`regions` / `fine_regions`、DOS の図の緑 / 青）には残る。DOS の図には [E_F − max_ewidth, E_F − min_ewidth] の帯（斜線）も描く。
6. **細かい部分区間が無いとき**（ギャップの底が dosth2 をわずかに超える場合。例: AlSiGeBi の底 1.04e-3）: dosth2 を `dosth2 × dosth2_relax`（既定 2 倍）に緩めて 2〜5 をやり直す。それでも無ければその粗い区間は候補を出さない。緩めた事実は `Decision.relaxed=True` として記録する。
7. 候補が 1 つも無ければ `fail`、あれば `new` と先頭候補。

`Decision` は (flag, ewidth, candidates, coarse_regions, fine_regions, relaxed) を持ち、`Judgement` に粗い区間と細かい区間の両方を残す。

Method 1 との違い:

| | Method 1 | Method 2 |
|---|---|---|
| threshold | 1 つ（dosth） | 2 つ（dosth 粗、dosth2 細） |
| 幅の条件 eth | dosth の区間に課す | dosth の区間に課す（細かい区間には課さない） |
| −ewidth の置き場 | dosth の区間内で上端から ediff 下 | dosth2 の部分区間内で、かつ上端から ediff 下 |
| semicore の裾 | dosth より低ければ区間に含めてしまう | dosth2 を超える裾は除外される（下端側の余裕が自動的に入る） |
| ギャップの底が dosth2 を超える系 | – | dosth2 を relax 倍まで緩めて拾う |

**Bi 系での値**（§13.2 の ewidth_go 1.6 の DOS、dosth 2e-2、dosth2 1e-3、eth 0.30、ediff 0.20）:

| 系 | 粗い区間（dosth 2e-2） | top = e2 − ediff | dosth2（使用値） | 細かい部分区間 | Method 2 の候補 | Method 1（2e-2） | Method 1（1e-3） |
|---|---|---|---|---|---|---|---|
| AlSiRhBi fcc | [−1.688, −1.073] | −1.272 | 1e-3 | [−1.387, −1.147] | **1.2825** | old 1.6 | fail |
| AlScNiBi fcc | [−1.673, −1.042] | −1.242 | 1e-3 | [−1.343, −1.117] | **1.2525** | old 1.6 | fail |
| AlSiGeBi fcc | [−1.673, −1.028] | −1.228 | 1e-3 では無し → **2e-3**（relaxed） | [−1.508, −1.042] | **1.2375** | old 1.6 | fail |

3 系とも −1.6（Bi 5d の裾）は細かい部分区間の外なので `new` になり、候補は 1.24〜1.28 Ry で、ギャップの底（DOS 最小 8e-4〜1e-3 の位置 −1.13〜−1.25 Ry）のすぐ下に落ちる。mesh 下端に接する粗い区間 [−2.243, −1.8〜−1.9]（Ge 3d / Bi 5d より下）からも候補 2.03〜2.12 が出るが順位は最後。

実装上の注意: `gap_regions` は DOS が負の点（窓の下端付近で数値的に −1e-3 程度になることがある）も threshold 未満として区間に入れる。粗い区間が mesh の下端に接するかどうかは `GapRegion.i1 == 0` で判定する。

## 5. 「まだ収束に向かっている」判定（`convergence.py`）

```python
def is_converging(err_history: Sequence[float], moment_history: Sequence[float], *,
                  last: int = 100, r2_th: float = 0.80, mae_err_th: float = 5e-4,
                  mae_moment_th: float = 5e-5, mae_err_loose_th: float = 6e-3) -> bool
```

- 直近 `last` 点の err（log10 rms）と moment に、反復番号を説明変数にした最小二乗直線を引き、標準化した値で r² を、生の値で MAE を求める。
- 判定は移植元の生きている条件だけ: `r2_err > r2_th`、または `mae_err < mae_err_th`、または `mae_moment < mae_moment_th and mae_err < mae_err_loose_th`。タイポで死んでいた条件は入れない。
- scikit-learn は使わず `numpy.polyfit` で書く（定数列は r² = 0）。点数が 3 未満なら `False`。

## 6. 1 ディレクトリの実行（`runner.py`）

```python
class KkrRunner:
    def __init__(self, akaikkr_exe: str, directory: str, param_go: dict, *,
                 ewidth_dos: float = 3.0, files: dict | None = None, compat: bool = False)
    def run_go(self, copy_potential_from: str | None = None, force: bool = False) -> bool  # converged
    def run_dos(self, force: bool = False) -> None
    def run_j(self, force: bool = False) -> None
    def run_all(self, copy_potential_from=None, with_j: bool = True) -> bool
    def result(self) -> dict     # §9 の 1 行分
```

- `files` の既定は `{"inputcard_go": "inputcard_go", "out_go": "out_go.log", "inputcard_dos": "inputcard_dos", "out_dos": "out_dos.log", "inputcard_j": "inputcard_j", "out_j": "out_j.log", "potential": "pot.dat"}`（移植元と同名）。
- `copy_potential_from` が与えられれば `AkaikkrJob.copy_potential` で pot.dat をコピーしてから実行する。
- 再利用の判定（誤り #1 の修正）: 既存の inputcard を読んで、これから書く内容と文字列として一致し、かつ出力が正常終了（`sbtime report` あり、`***err` 無し、`check_option_error` が None）していれば実行しない。inputcard が異なるときは既存出力を `out_go.log.bak-<n>` に退避して再実行する。`force=True` なら無条件に実行する。
- 実行は `AkaikkrJob.run`。`KKRFailedExecutionError` はそのまま上げる。収束は `AkaikkrJob.get_convergence`。
- `run_all` は go の後、`compat=False`（既定）なら収束したときだけ dos と j を、`compat=True` なら常に実行する。`with_j=False` で j を省く（磁性でない系の一括処理用）。
- dos と j の inputcard は `param_go` を `deepcopy` し `go` と（dos なら）`ewidth` だけ差し替える。`option` は `normalize_option` を通す。
- `result()` は `AkaikkrJob` の `get_total_energy, get_total_moment, get_component_moment, get_lattice_constant, get_unitcell_volume, get_rms_error, get_convergence, get_ewidth, get_option, get_emesh_param, get_curie_temperature(out_j), get_jij_as_dataframe(out_j), get_dos(out_dos)` で §9 の項目を返す。

### 6.1 ディレクトリ命名（`layout.py`）

```python
@dataclass
class RunPoint:
    key: str; iew: int; ewidth: float; ied: int; edelt: float; polytyp: str; ipm: int; pmix: float

class Layout:
    def __init__(self, prefix: str = "RUN", version: int = 2, sep: str = ",")
    def dirname(self, p: RunPoint) -> str
    def parse(self, dirname: str) -> RunPoint
    def find(self, **fixed) -> list[RunPoint]
```

- `key` は系を識別する文字列（既定は `SiteComposition.key()`、一般には呼び出し側が与えるラベル。2019 年の RUN では heakey）。`polytyp` は同じ系の構造バリアント（bcc/fcc など、1 つでもよい）。
- `version=1`（互換）: `key_{key},ew_{iew:03},ed_{ied:03},polytyp_{polytyp},pm_{ipm:03}`。移植元と同じ。2019 年の RUN を読むときに使う。
- `version=2`（既定）: `key_{key},ew_{iew:03}-{ewidth:.4f},ed_{ied:03}-{edelt:.0e},polytyp_{polytyp},pm_{ipm:03}-{pmix:.0e}`。ewidth / edelt / pmix の実値を含むので、同じ `iew` で ewidth が変わっても別ディレクトリになる。
- `find` は `os.scandir` で列挙して `parse` で戻す。区切りのカンマは移植元の資産と互換のため既定にし、`sep=";"` 等も許す。

## 7. GAES（`scheme.py`）

```python
class Gaes:
    def __init__(self, akaikkr_exe: str, layout: Layout, *,  # fresh_retry=True, bzqlty_steps=("+4",) も受ける
                 ewidth_init=1.2, ewidth_dos=3.0, ewidth_dos_auto=True, ref=0.75,
                 method=2, dosth=2e-2, dosth2=1e-3, eth=0.30, ediff=0.20, margin=0.01, min_ewidth=1.0, max_ewidth=2.0,
                 edelt_init=1e-4, edelt_steps=(1e-2, 1e-3, 1e-4), pmix_steps=(1e-2, 5e-3, 1e-3, 5e-4, 1e-4),
                 pmix_init=0.005, maxitr_init=500, maxitr_2nd=200, maxitr_pm=300,
                 max_pm_iter=20, max_ew=10, with_j=True, compat=False, logger=None)
    def run(self, key: str, params: dict[str, dict]) -> KeyResult   # polytyp -> param_go 辞書
    def run_many(self, items: Iterable[tuple[str, dict[str, dict]]]) -> list[KeyResult]
```

`params` は polytyp 名 → go の inputcard 辞書（`ewidth`, `edelt`, `pmix`, `maxitr` は scheme が上書きする）。単一サイト CPA なら `{"bcc": make_single_site_param(comp, "bcc"), "fcc": make_single_site_param(comp, "fcc")}`。1 キーの処理:

1. **STEP1（ewidth の粗い決定）**: `iew=0, ewidth=ewidth_init, edelt=edelt_init, pmix=pmix_init, maxitr=maxitr_init` で各 polytyp を 1 回だけ go → dos（→ j）。全 polytyp の total DOS（スピン和）を `gap_regions` に渡し `choose_ewidth`。
   - `old` かつ全 polytyp 収束 → **finished**。
   - `old` だが未収束 → STEP2 へ。
   - `new` → `iew += 1`、ewidth を候補に置き換えて STEP1 をやり直す（誤り #1 の修正。`version=2` の layout なら別ディレクトリになる）。同じ ewidth を二度試さない。`iew >= max_ew` で **ewidth_exhausted**。
   - `fail` → **ewidth_fail**（収束扱いにしない）。go が全 polytyp 収束していれば結果は残す。
2. **STEP2（SCF の追い込み）**: `maxitr=maxitr_2nd`。`edelt_steps` を**大きい方から**（既定 1e-2 → 1e-3 → 1e-4。移植元は 1e-4 → 1e-3 → 1e-2 だったが、2026-09-25 にユーザー指示で逆順にした。大きい edelt で滑らかにして収束させてから締める）、各 edelt で `pmix_steps` を順に試す。各 pmix では前段の pot.dat をコピーして `maxitr=maxitr_pm` で最大 `max_pm_iter` 回まで継続し、収束したら次の polytyp、`is_converging` が偽なら次の pmix。全 pmix を使い切ったら次の edelt。
   - 各 edelt の後に DOS を見て `choose_ewidth`。`old` かつ全 polytyp 収束 → **finished**。`new` → 新 ewidth で STEP1 からやり直し。`fail` → **ewidth_fail**。
3. `iew` が `max_ew` に達する、または候補が尽きたら **not_converged**。

- `KeyResult` は `status`（`finished | ewidth_fail | ewidth_exhausted | not_converged | error`）、最終 `RunPoint`（polytyp ごと）、試した ewidth の列、採用したバンドギャップ区間、§9 の結果行、例外メッセージを持ち、`RUN/key_<key>.json` に書く。1 キーの例外は捕まえて `error` にし、次のキーへ進む。
- 判定に使う曲線は既定で polytyp ごとの total DOS。`curves="pdos"` にすると各 polytyp の全成分 PDOS を連結して渡す（§0.3。`dosth_pdos` を使う。実装は第 2 段階）。
- `compat=True` のときの差: `Layout(version=1)`、dos / j を常に実行、`dosth=2e-2`、STEP1 の `new` で `iew` を進めずに既存結果を再利用（移植元の挙動の再現。2019 年の RUN の検証用で、新規計算には使わない）。

- **ポテンシャルの引き継ぎ**（2026-09-25 ユーザー規則、実装済み）:
  1. ewidth を変えたら古い pot.dat は使わず、新しい pot（specx の初期化）から始める。新 ewidth の STEP1 は pot.dat をコピーしない。
  2. ewidth を変えずに収束パラメタ（edelt、pmix）を変えるときは、前の pot.dat から続ける（移植元の動作）。
  3. それでも収束しなければ、同じ設定を新しい pot で試す（`fresh_retry=True`、既定。移植元にコメントアウトで残っていた `dir_first` 案に近い）。前段の pot から始めた連鎖が「まだ収束に向かっている」まま `max_pm_iter` を使い切ったときは fresh の再試行はせず、次の pmix に進む。
  4. それでも収束しなければ bzqlty を上げて STEP2 をやり直す（`bzqlty_steps=("+4",)`、既定は 1 段階。`"+n"` は現在値に加算、数値なら置き換え。ディレクトリの `ied` は `edelt_steps` の本数ずつ進む）。
  `compat=True` では 3 と 4 を行わない。

### 7.1 GAES-Committee（`committee.py`）

```python
class GaesCommittee:
    def __init__(self, gaes: Gaes, *, ewidth_members=(0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0), ewidth_dos=None,
                 dosth_members=(1e-2, 5e-3, 1e-3, 5e-4, 1e-4),
                 quorum="majority", n_parallel=4, threads_per_member=None,
                 maxitr_member=100, reuse_member=True)
    def run(self, key: str, params: dict[str, dict]) -> KeyResult
    def vote(self, energy: np.ndarray, member_regions: Sequence[Sequence[GapRegion]]) -> list[GapRegion]
```

1 キーの処理:

1. **Committee（投機的並列実行）**: `ewidth_members` の各 ewidth_go について、各 polytyp を go → dos で 1 回ずつ走らせる。**SCF の収束は要求しない**（`maxitr_member` の既定は 100 程度の小さい値にする）。最初の ewidth はバンドギャップの外にあることが多く、そのとき AkaiKKR の SCF はむしろ収束しにくい（§13 の「未収束でも判定する」と同じ理由）。収束を待つより、多くの ewidth で多数のギャップ区間を求めて vote に回す方が重要で、メンバー数（`ewidth_members` の本数）を増やす方向に計算資源を使う。**ewidth_dos は全メンバーで固定**し、`ewidth_dos=None` なら `(max(ewidth_members) + eth + ediff) / ref` を使う（§3.2 の条件を全メンバーで満たすため。窓が届かないメンバーの DOS は投票に参加できない）。`n_parallel` 個の specx を同時に走らせ、`OMP_NUM_THREADS` を `threads_per_member`（既定は総スレッド数 / n_parallel）に分ける。ディレクトリは `Layout` の `iew` にメンバー番号を使い、ewidth の実値が名前に入る（§6.1 version 2）。
2. **Vote（投票）**: 投票者は **ewidth メンバー × DOS threshold** の組（2026-09-25 ユーザー補足）。各メンバーの DOS に対して `dosth_members`（既定 1e-2, 5e-3, 1e-3, 5e-4, 1e-4）のそれぞれで `gap_regions` を求め、mesh 点ごとに「その点がギャップ区間に入っている」票を数える（8 メンバー × 5 threshold なら 40 票）。票数が quorum（`"majority"`: 投票者数の過半、整数: その数、`"all"`: 全員）以上の連続 mesh 区間を **consensus gap** とする。threshold を変えて投票するのは、threshold 1 つの値に判定が依存しないようにするためで、大きい threshold の票は区間を広く、小さい threshold の票は狭く取るので、過半で残るのは「ほとんどの threshold で低 DOS」の芯の部分になる。polytyp は AND（GAES と同じ）。ewidth_dos と mse が固定なので E − E_F の mesh は全メンバーで同一で、EF の差は軸に吸収される。投票の単位は mesh 点（区間の重なりではない）。
3. **仕上げ**: consensus gap に `choose_ewidth` の候補生成（区間上端 − ediff − margin）を当てて初期 ewidth を決め、それを `ewidth_init` として `Gaes.run` を呼ぶ。`reuse_member=True` で、採った ewidth がメンバーの ewidth と一致すればそのディレクトリの pot.dat と結果を STEP1 の出発点として再利用する（inputcard 一致の規則 §6）。
4. consensus gap が空なら **committee_fail** を記録し、`Gaes.run` を既定の `ewidth_init` で走らせる（GAES 単独に退避）。

- `KeyResult` に `committee` を加える: メンバーごとの ewidth、収束、threshold ごとのギャップ区間、mesh 点ごとの得票数、consensus gap、採用した初期 ewidth。票の割れ（メンバー間でギャップ区間が一致しない度合い）はギャップ推定の不確かさとして残す。
- メンバーの反復数は少なく、メンバー数は多く。1 メンバーの go は `maxitr_member` で打ち切り、DOS はその未収束ポテンシャルで計算する。仕上げの GAES（STEP2 の追い込み）だけが収束を目指す。
- 使い分け: GAES は 1 本ずつ直すので specx の実行回数は少ないが逐次。GAES-Committee はメンバー数だけ余分に走らせる代わりに初期 ewidth の当たりがよく、valence 帯の中や semicore の中から始めて何度もやり直す事態を避けられる。並列資源があるときは Committee、無いときは GAES。
- PDOS 拡張（§0.3）は vote の前段の `gap_regions` に PDOS 曲線を渡すだけで、Committee 側は変わらない。

## 8. 2019 年の RUN の読み直し（互換モード）

- `collect_legacy(prefix) -> pd.DataFrame`: `Layout(version=1)` で列挙し、各ディレクトリの inputcard_go から ewidth / edelt / pmix を読み直して（ディレクトリ名の `iew` は信用しない）、`KkrRunner.result()` で §9 の行を作る。
- 残っているのは最終 ewidth の計算だけなので、「ewidth の試行履歴」は復元できない。str.out があれば `('old'|'new'|'fail', ...)` 行から試行を読む `parse_str_out(path)` を補助として用意する。
- 検証: 例 `key_13487580,ew_000,ed_000,polytyp_fcc,pm_000`（AlCdReHg fcc）で `total_energy=-21075.683771993`、`a=8.01512`、`Tc=0.0`。key 13142122 の DOS に `gap_regions(dosth=2e-2)` を当てて str.out の `[-2.2425,-2.1075],[-2.0025,-0.7275]` が再現し、`choose_ewidth(1.2)` が `("old", 1.2, [0.9375])` になること。

## 9. 出力（結果行）

`KkrRunner.result()` と `KeyResult` の 1 polytyp 分:

| 列 | 取得元 |
|---|---|
| key, polytyp, ewidth, edelt, pmix, iew, ied, ipm, directory | `RunPoint` |
| ewidth_dos, ref, dosth, gap_regions, gap_used | scheme のパラメータと §4 の結果（`gap_used` は `check_ewidth` が返した区間） |
| a_bohr, volume_bohr3 | `get_lattice_constant`, `get_unitcell_volume`（out_go） |
| converged, n_iter, last_err | `get_convergence`, `len(get_rms_error)`, `get_rms_error()[-1]` |
| total_energy_Ry, total_moment, component_moment[] | out_go |
| option, emesh | `get_option`, `get_emesh_param`（out_go） |
| Tc_K, jij_csv | out_j（`with_j` のとき） |
| dos_csv | `get_dos(out_dos, "dataframe")` を保存したパス |

集計は `Gaes.collect(prefix) -> pd.DataFrame`（全 `key_*.json` を読む）。

## 10. CLI

```
kkr-gaes run   --exe /path/to/specx --input params.json [--prefix RUN] [--dosth 1e-3]
                      [--ewidth-init 1.2] [--ewidth-dos 3.0] [--ref 0.75] [--no-j] [--compat] [--threads 24]
                      [--committee 0.8,1.0,1.2,1.5,2.0 [--quorum majority] [--n-parallel 4]]   # GAES-Committee で初期 ewidth を決める
kkr-gaes site  --exe /path/to/specx --comp Rh0.5Pt0.5 [--comp AlSiScTi ...] [--polytyp bcc fcc]
                      [--lattice expr|mjw|<bohr>] ...（run と同じオプション）   # 単一サイト CPA を組成から
kkr-gaes hea   --exe /path/to/specx --keys heakeylist0.csv [--polytyp bcc fcc] [--start N --stop M] ...
                      # 2019 年の heakey リストから（互換用）
kkr-gaes check --dos out_dos.log [--dos out_dos_fcc.log ...] --ewidth 1.2 [--dosth 1e-3]
                      # 実行せず、既存の dos 出力に対して §4 の判定だけを出す
kkr-gaes collect --prefix RUN [--legacy] -o result.csv
kkr-gaes comp  Rh0.5Pt0.5 | AlSiScTi | 13142122   # 組成・type 名・heakey の相互変換と type 名の長さ検査
```

- `params.json` は `{key: {polytyp: param_go}}`。`site` と `hea` サブコマンドは組成（または heakey リスト）からこれを作って `run` と同じ処理をする。
- `check` は人間が手でやっていた判定の代替で、scheme を回さずに使える。
- `--threads` は `OMP_NUM_THREADS` を子プロセスにだけ設定する。ログは `prefix/kkr-gaes.log`。`--committee` を付けると `GaesCommittee`、無ければ `Gaes` を使う。

## 11. テスト（`tests/gaes/`）

specx 不要:

- `test_gap.py`: 合成 DOS（ガウス 2 山 + 定数 5e-4）で区間が 1 本、`dosth` を変えると幅が変わる、末尾まで低 DOS が続くケース（旧 NameError）、mesh 不一致で例外、2 曲線の AND。
- `test_ewidth.py`: `check_ewidth` / `choose_ewidth` の `old` / `new` / `fail`、`min_ewidth`。
- `test_legacy_dos.py`: `run0/RUN/key_13487580,...,polytyp_fcc,pm_000/out_dos.log` と str.out の値の再現（§8。ディレクトリがあるときだけ）。
- `test_convergence.py`: 単調減少列で `True`、雑音のみで `False`、`last` 未満の長さ。
- `test_layout.py`: v1 / v2 の往復、旧ディレクトリ名の `parse`。
- `test_composition.py`: `from_type_name("Rh0.5Pt0.5_1")` → `{Rh: 0.5, Pt: 0.5}`、`"B0.975Vc0.025"` → Z (5, 0)、`from_elements("AlSiScTi")` が等比、`type_name()` の往復、41 文字以上・空白・カンマで `ValueError`、`heakey_to_composition("13142122")` ⇄ `("Al","Si","Sc","Ti")`、奇数桁で `ValueError`、5 元。`make_single_site_param(..., type_name="HEA")` の inputcard が例の `inputcard_go`（AlCdReHg fcc）と `#` 行と空白の正規化を除いて一致。
- `test_result_legacy.py`: 例ディレクトリの out_go / out_dos / out_j を `KkrRunner.result()` で読む。
- `test_committee_vote.py`: 合成 DOS 5 本（ギャップ区間を少しずつずらしたもの、1 本は窓が届かず区間無し）× threshold 5 個で `vote` が過半の区間を返す、threshold が大きいほど票の区間が広いこと、`quorum="all"` で共通部分だけになる、全員不一致で空、`ewidth_dos=None` の自動決定が §3.2 の条件を満たす。
- `test_gap.py` にギザギザした合成 DOS（ギャップ内に幅 1 mesh の 5e-4 の凹凸、valence 帯に 1 mesh の落ち込み）を加え、threshold 比較だけで区間が正しく出ること（微分を使えば誤る形）。

specx が必要（`AKAIKKR_PROGRAM_PATH`）:

- `test_run.py`: Cu fcc（`lattice=6.82`、`tests/akaikkr/reference/ifort.json` の Cu_go と同じ ewidth / bzqlty）で `KkrRunner.run_all(with_j=False)` が収束し te が一致。`gap_regions` が E_F − 1.0 を含む区間を返し `check_ewidth(1.0)` が None でないこと（Cu の 3d 価電子帯の下に semicore は無い）。
- `test_gaes_small.py`: 1 系（bcc のみ、`max_pm_iter=2`, `maxitr_init=50`）で `Gaes.run` が例外なく `status` を返すこと。ewidth_init を意図的に valence 帯の中（例 0.3）にして `new` が出て別ディレクトリで再計算されること（誤り #1 の修正の確認）。

## 12. 実装時の注意

- `AkaikkrJob.make_inputcard` は `atmicx` を `["0.0a","0.0b","0.0c","HEA"]` の形で受ける。移植元の `0 0 0 HEA` と等価。
- `record=2nd` で pot.dat が無いときの `***wrn in spmain...eof detected` は正常。
- specx は `a=1000000` をヘッダで `a=*********` と印字する。`get_lattice_constant` は `bravais=` の行から読むので影響しない。
- DOS の閾値判定は E − E_F の軸で行う。`get_dos` が返すエネルギーは out_dos.log の 1 列目そのもの（E_F 基準）。
- 移植元の `pot.dat.info` は読まない。
- **ewidth を微少動かしただけで SCF が収束する / しないことがある**（ユーザー注記、2026-09-25）。積分路の下端 E_F − ewidth が準位や帯の端の近くにあると、わずかな ewidth の違いで下端が準位をまたぐ（`*` の付け替え、reconf の結果が変わる）、あるいは帯の裾を切る量が変わり、収束の可否が反転する。例: Hf 4f が −1.2 Ry にある系で ewidth 1.2 は収束せず 0.9 は収束する（§13.1）、Rb 4p core 指定で 0.5875（valence 帯の底 −0.58 に接する）は発散する（§15.5.1）。したがって、収束しないことだけを理由に「その ewidth 付近に解が無い」とは言えず、ギャップ判定は DOS で行う（§0.2）。逆に、収束したことは積分路がギャップにある証拠にならない（In 4d core の 0.641 は valence 帯の中で収束した）。候補 ewidth は必ずギャップの位置から決め、収束パラメタ（edelt、pmix、bzqlty）の変更と ewidth の変更は別の段階で行う（§7）。図 `docs/data/scf_convergence_vs_ewidth_dos.png`（各面の 2 行目に、その go で specx がその軌道を core と valence（`*`）のどちらで扱ったかを書く。AlSiSnHf bcc は cpa2021v01 ビルド: 0.9 収束 / 1.2 未収束、どちらも 4f は valence `*` で −0.91 / −1.20 Ry、積分路の下端に重なる。In 4d は 0.641 で core（−1.35 Ry）、1.2 で valence `*`（−1.11 Ry）。Tl 5d は 0.386 で core（−3.21 Ry、窓の外）、1.2 で valence `*`（−0.94 Ry）。InMnFeCo In 4d core: 0.641 は valence 帯の中で収束、1.2 も収束。TlMnFeCo Tl 5d core: 0.386 は発散し 5d が −3.2 Ry へ、1.2 は収束。Rb 4p core の 1.2 / 0.5875 は `hea_XMnFeCo_fcc_orbital_rules_dos.png`）。スクリプト `tests/gaes/tools/plot_conv_pairs.py`。
- DOS 図を出すときは go の ewidth の線を入れる（`DosEXPlotter(..., go_outfile="out_go.log")`、[dos_plot_ewidth_line.md](dos_plot_ewidth_line.md)）。
- 2019 年の RUN を読むときは、誤り #1 のため `ew_000` に「ewidth=1.2 で計算した結果」しか無い。ewidth は必ず inputcard_go から読み直す。
- `test_committee_small.py`: 同じ系で `GaesCommittee(ewidth_members=(0.5, 1.0, 1.5), n_parallel=3, maxitr_member=30).run` が `committee` の票と consensus gap を返し、採用した初期 ewidth がその区間に入ること。

## 13. 実装メモ（2026-09-25、GAES 逐次方式）

実装: `library/PyAkaiKKR/src/pyakaikkr/gaes/`（`gap.py`, `ewidth.py`, `convergence.py`, `layout.py`, `runner.py`, `composition.py`, `legacy_hea.py`, `scheme.py`, `cli.py`）、`tests/gaes/`、CLI `kkr-gaes`（`setup.cfg` の console_scripts）。GAES-Committee（§7.1）は未実装。

仕様からの差:

- **SCF が未収束でも DOS でバンドギャップを判定する**（2019 年どおり、設計上の前提）。ewidth がギャップに無いこと自体が density / potential の不整合を生んで収束を妨げるので、収束を待ってから判定するのは順序が逆になる。また ewidth を変えれば DOS もギャップ区間も変わり、初期 ewidth がギャップ内でも多くの場合は SCF が進むにつれてギャップの下端が下がる。`Gaes(tighten_before_fail=True)` にすると未収束の `fail` の前に STEP2 を挟むが、既定は False。
- **dos は go が未収束でも実行する**（`KkrRunner.run_all(dos_always=True)` を scheme が使う）。判定に DOS が要るため。j は収束したときだけ（`compat=True` なら常に）。
- **dos の窓が届かないときの縮小**: `KkrRunner.run_dos` は specx が dos で止まったとき（`***err in reconf`、§13.1）、`ewidth_dos` を 0.25 ずつ `(ewidth_go + ediff)/ref` まで狭めて再試行する。使った値は `KkrRunner.ewidth_dos_used` と `Judgement.ewidth_dos` に残る。
- **ref の検出**: dos の mesh の下端から実効 ref（`-emin/ewidth_dos`、半刻み補正）を求め、`Gaes(ref=...)` と 0.05 以上ずれればログに警告する（ビルドごとの ref 既定の違い、§13.1）。
- `SiteComposition.conc` は百分率が整数になるときだけ整数を書く（等比 4 元で `25`、2019 年の inputcard と一致）。

### 13.1 AkaiKKR 2022.0721 の各ビルドと `reconf`

2019 年の HEA 系（AlGeHfBi、AlSiSnHf、AlScTiHf）で GAES を試したときに分かったこと。

| ビルド | `reconf.f` | dos の ref | `begin_option` | Hf を含む go |
|---|---|---|---|---|
| akaikkr | あり、有効（core 状態の再配置を行う） | 0.75（`cemesr_ref=` で変更可） | あり | **`***err in reconf...ecor not found` で開始直後に停止**。Hf 4f の core 準位（−0.68 Ry、絶対値）が再配置の探索窓に入り、Newton 反復（20 回）が収束しない。ewidth 1.0〜2.5、ng / mse / mxl / reltyp / sdftyp / edelt / rmt を変えても同じ。Al, Si, Sc, Ti, Ge, Sn, Bi の単体は通る |
| akaikkr_cpa2021v01 | あり、**`supprs=.true.` で再配置を抑止**（ebtm より下は core のまま、上は valence）。`emrgn` は 0.2（akaikkr は 0.7） | **0.5 固定**（`m_optn` を使わないので `cemesr_ref=` 不可） | 無し | 通る。AlCdReHg fcc は akaikkr ビルドと同じ 148 反復で収束（全エネルギー −21075.6451 と −21075.6422） |
| akaikkr_cnd | 無し | 0.5 | あり | `displc` 必須 |

- `reconf.f` は 2020 年 10 月に加わったもので 2019 年版には無い。2019 年の RUN は再配置無しの計算であり、cpa2021v01 ビルド（抑止）がそれに近い。
- したがって `tests/gaes/test_gaes_run.py` は `GAES_SPECX_CODE`（既定 `akaikkr_cpa2021v01`）で実行ファイルを選び、dos の窓は ref=0.5 に合わせて `ewidth_dos=4.5`（下端 −2.25 Ry、2019 年の ewidth_dos 3.0 × ref 0.75 と同じ）にする。akaikkr ビルドで Hf 系を走らせるには `reconf.f` の `supprs` 相当（コメントで用意されている 7 行）を有効にして作り直す必要がある。
- akaikkr ビルドでも、dos の窓が深い core 準位を跨ぐと同じ `reconf` で止まる（AlCdReHg: ewidth_dos 3.0 は停止、2.5 は通る）。上の「窓の縮小」はこのため。
- 2022.0721 では go の既定 mse / ng が 2019 年版と違う（Cu 相当で mse 43 / ng 21、2019 年は 65 / 15）ので、全エネルギーは 2019 年の値と 0.04 Ry 程度ずれる。`begin_option` の `mse=`, `ng=` で揃えられる。
- inputcard の末尾に改行が無い状態で `begin_option` を追記すると原子行に連結され `ty2ity...type not defined` になる。`make_inputcard` は末尾改行を付けないので、手で追記するときは注意。

#### 13.1.1 Hf 4f と reconf の詳細（2026-09-25、cpa2021v01 ビルド、AlSiSnHf bcc）

out_go.log から分かること:

- Hf の core 配置は 4f が up 7 / down 7（14 電子すべて core）。入力の `lmxtyp=2` では valence 基底が s, p, d までなので、`valence charge in the cell` に (f) は無く、DOS / PDOS に 4f の帯は現れない。
- `cstate.f` の `rearth=(nclr-57)*(nclr-71)*(nclr-89)*(nclr-103) .le. 0` は La〜Lu、Ac〜Lr だけを真にし、Hf（Z = 72）は偽。したがって「open f-core」（`reconf.f` 先頭の `opnc .and. rearth .and. jj .eq. 4` で占有を固定する分岐）は Hf には効かず、`noc` / `nopnc` も無関係。Hf 4f は他の semicore と同じく、準位が ebtm = E_F − ewidth より上に来ると reconf の再配置対象になる。2019 年の RUN には「open f-core」の表示も `4f)*` も無い（reconf 以前の版）。
- `total charge` は muffin-tin 球内の core 電荷 + セル内の valence 電荷で、球の外の core の裾は含まない。Hf の 70.7〜70.9（Z = 72）、Sn の 48.9〜49.1（Z = 50）は 2019 年の収束した走行でも同じで、欠損ではない。脱占有の量は `core charge in the muffin-tin sphere` の差で見る。

ewidth を振った結果（2022.0721 cpa2021v01、E_F 基準の 4f 準位は絶対値 − ef）:

| lmxtyp | ewidth | ebtm（絶対値） | 4f 準位（絶対値、up / down） | `*` | 球内 core 電荷 | valence (f) | SCF |
|---|---|---|---|---|---|---|---|
| 2 | 0.9 | −0.357 | −0.370 / −0.370 | down | 67.925 | — | 収束（152） |
| 2 | 1.15 | −0.635 | −0.785 / −0.732 | down | 67.944 | — | 未収束（rms 0.11） |
| 2 | 1.2 | −0.676 | −0.678 / −0.678 | 両方 | 67.941 | — | 未収束 |
| 2 | 1.5 | −1.013 | −1.151 / −1.058 | down | 67.955 | — | 未収束 |
| 3 | 0.9 | −0.355 | −0.368 / −0.368 | 両方 | 67.925 | 0.005 / spin | 収束（168） |
| 3 | 1.15 | −0.626 | −0.816 / −0.756 | down | 67.945 | 0.016 | 未収束（rms 0.19） |
| 3 | 1.2 | −0.679 | −0.871 / −0.806 | 両方 | 67.947 | 0.016 | 未収束（rms 0.22） |
| 3 | 1.5 | −1.022 | −1.081 / −1.004 | down | **60.956** | 0.015 | 発散（rms 1.0） |

読み方:

- 4f 準位は ewidth を変えると ebtm に追随して動く（0.9 で −0.37、1.2 で −0.68〜−0.87、1.5 で −1.0〜−1.15）。部分的に空いた 4f がポテンシャルを変え、準位が積分路の下端に張り付く自己無撞着な帰結で、2019 年版の −1.84 Ry（E_F 基準）とはまったく違う位置にある。
- 準位が ebtm のすぐ下にある間（0.9、1.15、1.2）は、reconf の幅 3 mRy のなましで脱占有は 0.02〜0.03 電子にとどまり、0.9 は収束、1.15 と 1.2 は振動する（1.15 は 2026-09-25 のユーザー指定で追加。up / down の 4f 準位が 0.05〜0.06 Ry 割れ、down 側だけ `*`）。1.5（lmxtyp=3）では down スピンの 4f が ebtm より上に出て 7 電子がまるごと core から外れ（球内 core 電荷 60.96）、valence の f には 0.015 電子しか入らないので発散する。
- **k 点を増やしても変わらない**（2026-09-25、ユーザー指定、ewidth 1.15、lmxtyp=2）: bzqlty 10 / 14 / 18（nk = 256 / 624 / 1240）で SCF の履歴は 500 反復まで 7 桁一致し（190 反復目 te = −10894.003679、rms 0.766）、最後の rms 0.106、球内 core 電荷 67.9438、4f の `*` も同じ。振動は k 点の粗さではなく、ebtm に張り付いた 4f の reconf に由来する。SCF 履歴（`neu=` 行）を見ると、6 反復ごとに neu（電子数の誤差）が −1.2〜−1.8 に跳び、全エネルギーが 13 Ry 変わる（itr 7, 13, ...）。4f の占有の再配置がその反復で反転し、SCF を毎回やり直している形で、bzqlty を変えても反復番号まで同じ。図 `docs/data/hf_bzqlty_AlSiSnHf_ew1.15_dos.png`（DOS と SCF 履歴）、出力 `tests/gaes/RUN_Hf_stopped/AlSiSnHf_bcc_lmxtyp3/lmx2_ew_1.15_bz*/`。
- `lmxtyp=3` で f を valence 基底に入れても、4f は core 配置のままなので valence の f 帯にはならず、収束の可否も変わらない。**Hf を含む系を 2022.0721 で扱うには、4f が ebtm から十分離れる小さい ewidth（この系では 0.9 が収束）を選ぶか、§15 の `Hf4f=core` 指定で max_ewidth を |E_4f − E_F| − ediff に抑えるしかない**。ただし準位が ebtm に追随するため、`core` 指定でも 1 回の go ごとに範囲を決め直す必要がある（§15.2 の「見た準位のうち最も浅い値」）。

図 `docs/data/hf_lmxtyp3_AlSiSnHf_dos.png`（左 lmxtyp=2、右 lmxtyp=3、上から ewidth 0.9 / 1.15 / 1.2 / 1.5。紫 = Hf の f-PDOS × 濃度）、出力 `tests/gaes/RUN_Hf_stopped/AlSiSnHf_bcc_lmxtyp3/ew_*/`、スクリプト `tests/gaes/tools/plot_hf_lmxtyp3.py`。

### 13.2 Bi 系のテスト（2026-09-25）

Hf 系（§13.1）の代わりに、RUN/ の走査で見つけた Bi 5d semicore の 3 系 AlSiRhBi fcc（key 13144583）、AlSiGeBi fcc（13143283、Ge 3d も）、AlScNiBi fcc（13212883）を、初期 ewidth 1.6（Bi 5d ピーク −1.75 Ry の肩）で走らせた（`tests/gaes/test_gaes_run.py`、`GAES_EWIDTH_INIT`、`GAES_DOSTH`、akaikkr ビルド、ref 0.75、ewidth_dos 3.0）。

- 3 系とも STEP1 の go は 123〜128 反復で収束（全エネルギー −13423.2853、−12080.5787、−12028.6003 Ry）。
- Method 1、dosth 2e-2: 3 系とも `old`（ギャップ [−1.69, −1.07] が 5d の肩まで届き −1.6 を含む）→ `finished` at 1.6。
- Method 1、dosth 1e-3: 3 系とも `fail`（区間幅 0.23〜0.24 < eth、AlSiGeBi は底 1.04e-3 で区間無し）。
- Method 2（§4.3、実装済み、既定）で ewidth 1.6 から走らせた結果: 3 系とも STEP1 で `new` → 新 ewidth の STEP1 で `old`、合計 2 回の go で `finished`。

  | 系 | ewidth の経過 | dosth（粗） | dosth2（使用値） | 最終ギャップ（細）| 収束 | 全エネルギー (Ry) |
  |---|---|---|---|---|---|---|
  | AlSiRhBi fcc | 1.6 → 1.2825 | 2e-2 | 1e-3 | [−1.373, −1.133] | 125 反復、err −6.16 | −13423.2973 |
  | AlSiGeBi fcc | 1.6 → 1.2375 | 2e-2 | 2e-3（1e-3 では細かい区間が無く 2 倍に緩和）| [−1.478, −1.042] | 120 反復、err −6.04 | −12080.6142 |
  | AlScNiBi fcc | 1.6 → 1.2525 | 2e-2 | 1e-3 | [−1.312, −1.103] | 136 反復、err −6.11 | −12028.6194 |

  ewidth 1.6 の解と比べて全エネルギーは 0.012〜0.036 Ry 低い（Bi 5d の裾を積分路が跨がなくなった分）。新 ewidth の DOS では細かい区間の位置が 0.01〜0.03 Ry ずれるだけで、候補は変わらない（`old`）。

### 13.3 X-Mn-Fe-Co fcc の 18 系（2026-09-25、Method 2、min_ewidth 1.0 / max_ewidth 2.0）

§14.1.1 の表から希ガスと reconf で止まる元素（Hf など）を除いた X = Ga, As, Se, Rb, In, Sb, Te, Ba, La, Ce, Pr, Nd, Pm, Sm, Yb, Lu, Tl, Bi について、X-Mn-Fe-Co 等比 fcc（単一サイト CPA、akaikkr ビルド、ref 0.75、ewidth 初期値 1.2）を走らせた。スクリプト `tests/gaes/tools/run_hea_XMnFeCo.py`、図 `tests/gaes/tools/plot_hea_XMnFeCo.py`、結果 `docs/data/hea_XMnFeCo_fcc_gaes_summary.json`、図 `docs/data/hea_XMnFeCo_fcc_gaes_dos.png`（各成分の core 準位を重ね描き、斜線は [E_F − max_ewidth, E_F − min_ewidth]）。

§4.3 の 2 つの修正（範囲外の候補を端へ寄せる、窓の下端に接する区間も使う）の後、18 系すべてが `finished`、SCF 収束。

| X | ewidth | 経過 | 使ったギャップ (Ry) | X の浅い準位 (E − E_F) |
|---|---|---|---|---|
| Ga | 1.416 | 1.2 → 1.416 | [−3.16, −1.32] | 3d* −0.97 |
| As | 1.324 | 1.2 → 1.324 | [−2.35, −1.14] | 4s* −0.89 |
| Se | 1.414 | 1.2 → 1.414 | [−3.03, −1.25] | 4s* −1.04 |
| Rb | 1.2 | old | [−1.93, −1.18] | 4p* −0.94 |
| In | 1.506 | 1.2 → 1.439 → 1.506 | [−3.16, −1.45] | 4d* −1.13 |
| Sb | 1.283 | 1.2 → 1.283 | [−1.90, −1.10] | 5s* −0.84 |
| Te | 1.403 | 1.2 → 1.403 | [−2.24, −1.22] | 5s* −0.99 |
| Ba | 1.388 | 1.2 → 1.373 → 1.388 | [−1.99, −1.34] | 5p* −1.07 |
| La | 1.673 | 1.2 → 1.673 | [−2.24, −1.58] | 5p* −1.32 |
| Ce | 1.2 | old | [−1.25, −0.79] | 5p −1.56 |
| Pr | 1.0 | 1.2 → 1.0（下限に寄せた） | [−1.09, −0.79] | 5p −1.38 |
| Nd | 1.0 | 1.2 → 1.0（下限に寄せた） | [−1.13, −0.77] | 5p −1.41 |
| Pm | 1.0 | 1.2 → 1.0（下限に寄せた） | [−1.18, −0.77] | 5p −1.44 |
| Sm | 1.2 | old | [−1.22, −0.76] | 5p −1.47 |
| Yb | 1.2 | old | [−1.36, −0.74] | 5p −1.63 |
| Lu | 1.2 | old | [−1.70, −0.74] | 5p −1.95 |
| Tl | 1.304 | 1.2 → 1.304 | [−3.16, −1.27] | 5d* −0.91 |
| Bi | 1.373 | 1.2 → 1.373 | [−1.42, −1.22] | 5d −1.80、6s* −0.95 |

- 図の緑（粗い区間）と青（細かい部分区間）は、保存した判定ではなく最終 dos を現在の規則で判定し直して描いている。Ce / Sm / Yb / Lu の RUN は窓の下端に接する区間を調べない旧版で走ったため、保存した判定には 5p 帯の下の細かい部分区間（Ce [−2.24, −1.82]、Sm [−2.24, −1.81]、Yb [−2.24, −1.87]）が無い。判定し直すと候補に 1.91 / 1.91 / 1.90 が加わるが、1.2 が上のギャップに入っているので `old` のまま。
- Pr / Nd / Pm は 5p 帯の上のギャップ [−1.1〜−1.2, −0.77] の候補 0.78〜0.79 が `min_ewidth` 1.0 に寄せられた解。5p 帯の下（1.78〜1.87）に置きたければ `min_ewidth` を 1.2 にする。
- Te / In / As / Bi の −0.7 Ry 付近（X の浅い s / d 帯と valence 帯の底の間）の低 DOS 区間は、DOS < 2e-2 の幅が 0.16〜0.29 Ry で eth = 0.30 に届かず、DOS < 1e-3 の点も無いので、ギャップ区間にならない。`min_ewidth` 0.5 で再実行しても同じ解になる（Te 1.4025、In 1.5025、As 1.3125、Bi 1.3725）。In は幅 0.29 で境界にあり、eth 0.25 なら In 4d を core 側に置く解（≈ 0.8）が候補になる。

## 14. 実行中に原子準位から ewidth を予測する方法（提案、2026-09-25）

### 14.1 使える事実

- out_go.log の成分ごとの core 準位（`core level` 行）は、その準位が dos の窓に入っていれば total DOS の semicore ピークと数 mRy 以内で一致する（[core_levels_vs_dos.md](core_levels_vs_dos.md) §2、fcc 7,505 系）。つまり **semicore の位置は DOS を計算しなくても go の出力から分かる**。
- core 準位は SCF の各反復で更新されるが、出力には最後の 1 回しか出ない。E_F（`ef=`）も同じ。ただし反復ごとの `neu`（電荷中性のずれ）と `te` は出るので、SCF が「ある程度」進んだかは判定できる。
- `*` 付きの準位（ebtm より上で valence に切替）は ewidth に応じて 0.05〜0.1 Ry 動く（In 4d: ewidth 1.2 → 1.46 → 1.53 で −1.10 → −1.17 → −1.13 Ry）。core のままの準位（Ge 3d、Sn 4d、Bi 5d …）は系や ewidth によらず ±0.03 Ry しか動かない。
- valence 帯の底は DOS でしか分からない。ただし 4 元 HEA では −0.75〜−1.0 Ry の範囲に収まっていた（valence 帯の幅は成分の s, p 帯で決まる）。

### 14.1.1 元素ごとの準位表（段階 0: 計算前の予測に使う）

ユーザー提案（2026-09-25）: 準位表を §14.1 に加え、まず E_F を仮定して ewidth の予測値を得る。表は 2 種類ある。

**A. 収束済みの表（主表）**: 2019 年の RUN（fcc 7,505 系）の out_go.log から取った、成分ごとの core 準位の **E − E_F**（[core_levels_vs_dos.md](core_levels_vs_dos.md) §2、`docs/data/converged_core_levels_2019.csv`）。既に E_F 基準なので E_F の仮定は要らない。系による散らばりは ±0.05 Ry。dos の窓（−2.3 Ry）に入る準位だけを載せる。

| 元素 | 準位 | E − E_F 中央値 (Ry) | 範囲 | 扱い |
|---|---|---|---|---|
| Tl | 5d | −0.90 | −0.95〜−0.79 | valence（`*`） |
| Bi | 6s | −0.88 | −0.94〜−0.80 | valence（`*`） |
| Ga | 3d | −1.10 | −1.17〜−1.01 | valence（`*`） |
| In | 4d | −1.10 | −1.17〜−1.04 | valence（`*`） |
| Pb | 5d | −1.29 | −1.35〜−1.23 | core |
| Sn | 4d | −1.62 | −1.67〜−1.56 | core |
| Y | 4p | −1.68 | −1.72〜−1.64 | core |
| Bi | 5d | −1.73 | −1.78〜−1.67 | core |
| Ge | 3d | −1.81 | −1.84〜−1.79 | core |
| Hf | 4f | −1.84 | −1.89〜−1.80 | core（2019 年版。2022 版では −1.20） |
| Zr | 4p | −1.98 | −2.02〜−1.96 | core |
| Sc | 3p | −2.06 | −2.12〜−2.03 | core |
| Hf | 5p | −2.26 | −2.29〜−2.23 | core |
| Nb | 4p | −2.28 | −2.30〜−2.26 | core |

2019 年の 37 元素のうちここに無い元素（Al, Si, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Mo, Tc, Ru, Rh, Pd, Ag, Cd, Ta, W, Re, Os, Ir, Pt, Au, Hg）は −2.3 Ry より浅い core 準位を持たない（Cd 4d、Hg 5d は valence 配置）。

**B. 初期原子ポテンシャルの表（補助）**: H〜Bi の単体 fcc に `go=dsp` を新規で掛けた絶対値の準位（`docs/data/atomic_core_levels_dsp.csv`、図 `atomic_core_levels_dsp.png`）。A に無い元素（アルカリ、アルカリ土類、La / Ce、希ガス、非金属）に使う。絶対値なので **E_F = 0.5〜0.7 Ry** を仮定して E − E_F に直す（テスト物質の E_F は 0.42〜0.87 Ry、HEA は 0.44〜0.60 Ry）。初期値は収束値より深く出る（Pb 5d −1.86 → −1.29、Tl 5d −1.37 → −0.90、Bi 6s −1.23 → −0.88、Sn 4d −1.89 → −1.62、In 4d −1.34 → −1.10、Ga 3d −1.24 → −1.10、Hf 4f −1.33 → −1.20）ので、禁止区間の上側を 0.3〜0.5 Ry 広げる。

**段階 0 の規則**（dos も go も無しに ewidth の初期値を決める）:

1. 組成の各元素について、表 A（無ければ B）の準位 E_c（E − E_F）を取り、禁止区間 [E_c,min − w, E_c,max + w]（w = 0.15 Ry。表 B なら E_F の範囲 0.5〜0.7 を含めてさらに上側 +0.3）を作る。
2. ewidth の候補範囲 [min_ewidth, max_ewidth] = [1.2, 2.0]（下限は valence 帯の底 −0.75〜−1.0 Ry に ediff 0.2 を取った値）を 0.01 刻みで走査し、−ewidth がどの禁止区間にも入らない連続区間を **許容区間** とする。
3. 許容区間のうち最も浅いものの下端を初期 ewidth にする。許容区間が無ければ「semicore の階段」（AlGaSnPb 型）として予告し、Method 2 で `fail` になる可能性が高いと報告する。

**検証（2026-09-25 の GAES 結果との比較、表 A、w = 0.15）**:

| 系 | 禁止区間 (E − E_F) | 許容 ewidth | 段階 0 の予測 | GAES（Method 2）の結果 |
|---|---|---|---|---|
| AlSiRhBi | Bi 5d [−1.93, −1.52], Bi 6s [−1.09, −0.65] | [1.2, 1.52], [1.94, 2.0] | 1.2 | 1.2825 |
| AlSiGeBi | Ge 3d [−1.99, −1.64], Bi 5d, Bi 6s | [1.2, 1.52], [1.99, 2.0] | 1.2 | 1.2375 |
| AlScNiBi | Sc 3p [−2.27, −1.88], Bi 5d, Bi 6s | [1.2, 1.52] | 1.2 | 1.2525 |
| AlSiScPb | Sc 3p, Pb 5d [−1.50, −1.08] | [1.51, 1.87] | 1.51 | 1.6525 |
| AlMnFeIn | In 4d [−1.31, −0.89] | [1.32, 2.0] | 1.32 | 1.5288 |
| AlCoNiSn | Sn 4d [−1.82, −1.41] | [1.2, 1.40], [1.83, 2.0] | 1.2 | 1.2 |
| AlSiGeTl | Ge 3d, Tl 5d [−1.10, −0.64] | [1.2, 1.63], [1.99, 2.0] | 1.2 | 1.3525 |
| AlGaSnPb | Ga 3d [−1.31, −0.86], Sn 4d, Pb 5d | [1.83, 2.0] | 1.83（階段の予告） | fail |

8 系すべてで GAES の最終 ewidth は段階 0 の許容区間の中にあり、AlGaSnPb は許容区間が深い側にしか無いことから `fail` を予告できる。段階 0 の予測値（許容区間の下端）は GAES の値より 0.03〜0.2 Ry 浅いが、これは valence 帯の底との ediff を DOS で合わせる分で、§14.3 の B（dos 1 回）で埋まる。

段階 0 → A（短い go で core 準位を実測し表を更新）→ B（dos 1 回で valence 帯の底に合わせる）の順に精度を上げる。

### 14.2 「DOS がある程度収束したら E_F と semicore 位置が分かるか」

分かる。手順は次のとおり。

1. go を `maxitr` を小さく（50〜100）して回す。収束していなくても、`ef=` と core 準位は出る。`neu` の絶対値が 0.1 未満まで下がっていれば、E_F は最終値から ±0.05 Ry 以内、core 準位（core のもの）は ±0.03 Ry 以内にある（今日の In / Sn / Bi の実行での最終値との比較から。数値は実装時に検証する）。
2. 成分ごとの core 準位 E_c − E_F を集め、dos の窓に入る範囲（−2.3 Ry 以上）にあるものを **semicore の候補位置** とする。`*` の有無も記録する。
3. semicore の位置が分かれば、ギャップの**下端**はその準位から上に 0.1〜0.2 Ry（ピークの裾、DOS < 1e-3 になる位置。d 状態は幅 0.1、p 状態は 0.2 程度）で予測できる。ギャップの**上端**（valence 帯の底）は DOS が要る。

### 14.3 提案: core 準位を使った 3 段の予測

**A. 予測（dos 無し、go の短い実行 1 回）**

- `maxitr_probe`（100）で go を 1 回回し、`neu` が閾値を切ったら止める（AkaiKKR には反復途中で止める入力が無いので、`maxitr` を小さくして走らせ、`neu` が大きければもう 1 回続ける）。
- core 準位から semicore 候補 {E_c} を作る。**禁止区間**を E_c ± w（w = 0.15 Ry、`*` 付きは 0.25 Ry）とする。
- ewidth の候補は、E_F から下へ見て、禁止区間に入らず、かつ上の禁止区間（または E_F − ewidth_min）から ediff 以上離れた位置。具体的には、E_F 基準で最も浅い禁止区間の上端 e_top を取り、`ewidth_pred = −(e_top + ediff + margin)` … ただし e_top より上に valence 帯の底があることは DOS で確かめないと分からないので、A の段階では「semicore に当たらない範囲」だけを決める。
- 禁止区間が無い（−2.3 Ry より上に core 準位が無い: Al と 3d 遷移金属だけの系、Pt / Au / Hg など）なら、ewidth は valence 帯の底だけで決まる。この場合 A は「制約無し」を返す。

**B. 確認（dos 1 回）**

- A の候補で go（通常の maxitr）→ dos を回し、Method 2 で判定する。A で禁止区間を避けているので、Bi 系の「5d の裾に乗る」や In 系の「4d を切る」ような初期値の失敗が無くなり、`new` は valence 帯の底の位置調整だけになる。
- dos で得た valence 帯の底 e_v と、A の禁止区間の上端 e_top の間がギャップ。幅 e_v − e_top < eth なら `fail`（AlGaSnPb のような semicore の階段）を DOS を待たずに A の段階で予告できる（e_v は −0.75〜−1.0 Ry と仮定して警告する）。

**C. Committee への組み込み**

- GAES-Committee（§7.1）のメンバー ewidth は、A の禁止区間を避けた値だけにする。メンバー数を減らせる。
- vote には DOS 由来のギャップ区間に加え、core 準位由来の禁止区間の補集合を 1 票として入れる（threshold に依存しない票）。

### 14.4 実装の場所

- `pyakaikkr`: `AkaikkrJob.get_core_levels_by_component(outfile)`（成分・軌道・値・`*`）、`gaes.corelevels.forbidden_regions(job, outfile, w, w_star)`、`gaes.ewidth.decide` に `forbidden` を渡して候補と `old` 判定から除外する。
- `Gaes` / `AkaikkrGaesWorkChain`: STEP1 の前に probe go（`maxitr_probe`）を入れ、`ewidth_init` を A で置き換える（`predict_from_core_levels=True`）。
- 検証: 今日の 4 系（AlSiRhBi、AlSiScPb、AlMnFeIn、AlCoNiSn）と Bi / Pb / Tl の 2019 年 RUN で、A の予測が Method 2 の最終 ewidth と一致するか、go の回数が減るかを見る。

### 14.5 限界

- valence 帯の底は core 準位からは出ないので、DOS を完全に省くことはできない。A は「どこに置いてはいけないか」を決める段。
- `*` 付きの準位は ewidth で動くので、A の予測は ±0.1 Ry の不確かさを持つ。禁止区間の幅で吸収する。
- 2022.0721 の akaikkr ビルドでは Hf 4f のように準位そのものが 2019 年版と 0.6 Ry 違う場合がある。予測は同じビルドの go 出力から行う。

## 15. 軌道の valence / core 指定から ewidth の範囲を決める（仕様、2026-09-25、branch `ewidth_select_orbital`）

ユーザー要望: 「ある元素のある軌道を occupied / unoccupied にするように指定できるか。軌道準位が分かるので min_ewidth が決まるはず。例: Rb 4p、Se 4s、Bi 6s」。

### 15.1 用語と意味

積分路 [E_F − ewidth, E_F] に入る状態は valence として自己無撞着に扱われ、積分路より下（ebtm より下）の core 配置の状態は core として固定される（§0、[core_levels_vs_dos.md](core_levels_vs_dos.md) §1 の `*`）。ユーザーの語をこれに対応させる。

| ユーザーの語 | 正準名 | 意味 | specx 出力での確認 |
|---|---|---|---|
| occupied | `valence` | その軌道を積分路に入れる（E_F − ewidth を準位より **下** に置く） | out_go.log の core level 表でその準位に `*` が付く（core 配置に無い軌道は表に出ず、常に valence） |
| unoccupied | `core` | その軌道を積分路から外す（E_F − ewidth を準位より **上** に置く） | 同じ表で `*` が付かない |

「occupied / unoccupied」は「積分路が占める / 占めない」の意味であり、電子の占有数のことではない（core も占有されている）。API は `valence` / `core` を正準とし、`occupied` / `unoccupied` を別名として受け付ける。

指定できる軌道は、その元素の **core 配置にある軌道**（out_go.log の `core configuration` 表、`docs/data/atomic_core_levels_dsp.csv` に載る軌道）だけ。core 配置に無い軌道（Cd 4d、遷移金属の 3d など）は常に valence なので、`valence` 指定は無条件に満たされ（制約無し）、`core` 指定は満たせない（`GaesError`: specx には ewidth 以外に軌道を core に固定する手段が無い）。

### 15.2 規則: 準位から範囲を決める

軌道 (元素 X, 軌道 nl) の準位を ε = E − E_F（< 0）とする。

- `valence`: ewidth ≥ |ε| + ediff。準位のすぐ下に積分路の下端を置くと、準位の裾（`*` 準位は幅 0.05〜0.1 Ry のピーク）を切るので、ediff（既定 0.2 Ry）の余裕を取る。→ **min_ewidth_orb = |ε| + ediff**。
- `core`: ewidth ≤ |ε| − ediff。→ **max_ewidth_orb = |ε| − ediff**。
- 複数の指定は AND: min_ewidth = max(min_ewidth_orb)、max_ewidth = min(max_ewidth_orb)。
- 軌道指定から出た範囲は **既定の [1.0, 2.0] を置き換える**（既定値は「指定が無いときの規定」）。ユーザーが `min_ewidth` / `max_ewidth` を明示して同時に与えた場合は両方を満たす範囲（共通部分）を使い、空なら `GaesError` で両方の値を示す。
- min_ewidth > max_ewidth になったら（例: Rb 4p を core、Se 4s を valence）`GaesError`。矛盾した組を示す。
- 範囲を決めた後の動作は §4.3 と同じ: ギャップ区間の判定は範囲に依らず窓全体で行い、ewidth の **選択** だけを範囲で絞る。範囲内にギャップの候補が無ければ `fail`（message に範囲と、範囲内の DOS の最小値を書く）。§4.3 の `old` 判定に **−ewidth が範囲内にあること** を条件として加える（現状は範囲を見ずに `old` を返すので、例えば Rb 4p を core と指定しても現在の 1.2 がギャップに入っていれば `old` になってしまう。§15.5 の表の「現状」列）。

準位は次の順で取る。上のものが無ければ下へ。

1. **実行中**: 直前の go の out_go.log の成分ブロックから、その元素・軌道の準位（E − E_F = 準位 − ef、`*`）。判定のたびに読み直す（`*` 準位は ewidth で 0.05〜0.1 Ry 動く、§14.5）。SCF が収束していなくても使う（core 準位は数十反復で mRy まで落ち着く、§14.2）。
2. **段階 0（最初の go の前）**: `docs/data/converged_core_levels_2019.csv`（E − E_F の中央値、§14.1.1 表 A）。表 A に無い元素・軌道は `docs/data/atomic_core_levels_dsp.csv`（表 B、muffin-tin ゼロ基準の絶対値）から ε ≈ energy_Ry − ef_assumed（既定 0.6 Ry）で見積もる。表 B は原子の初期ポテンシャルの値で、収束値より 0.2〜0.5 Ry 深い（Rb 4p: 表 B −0.862 Ry、18 系の E_F 0.245 で ε ≈ −1.11 に対し収束値 −0.94）。したがって表 B から決めた範囲は **初期 ewidth の選択にだけ** 使い、最初の go の後は必ず 1. の値で決め直す。表 A の値も系による散らばり ±0.05 Ry があるので同様に決め直す。

初期 ewidth（`ewidth_init`）は範囲の中に無ければ範囲の中央に置き換える。

### 15.3 確認と再判定

各 go の後に、指定した軌道の `*` を読み、指定と一致するか確かめる。

- `valence` 指定なのに `*` が無い（積分路が準位より上に来た）、または `core` 指定なのに `*` が付く → 準位が動いて範囲の外に出た。範囲を 1. の準位で決め直し、`new` として次の ewidth へ（`max_ew` まで）。history に `orbital_mismatch` を記録する。
- 一致していれば通常の判定。

`finished` の条件に「指定した全軌道の `*` が指定と一致」を加える。

### 15.4 API

pyakaikkr（`pyakaikkr.gaes.orbital`、新規）:

```python
@dataclass
class OrbitalRule:
    element: str        # "Rb"
    orbital: str        # "4p"
    role: str           # "valence" | "core"

parse_orbital_rules(["Rb4p=occupied", "Se4s:core", "Bi 6s valence"]) -> [OrbitalRule, ...]
    # 区切りは "=", ":", 空白のどれでもよい。元素記号 + n + l（s/p/d/f）。role は
    # valence|occupied|core|unoccupied（大文字小文字を区別しない）。
levels_from_go(outfile) -> {(element, orbital): (E_minus_EF, star)}     # AkaikkrJob.get_core_levels_by_component
levels_from_tables(elements, ef_assumed=0.6) -> {(element, orbital): E_minus_EF}   # 表 A、無ければ表 B
bounds_from_rules(rules, levels, ediff=0.2, user_min=None, user_max=None)
    -> Bounds(min_ewidth, max_ewidth, details=[(rule, level, bound), ...])   # 矛盾は GaesError
check_rules(rules, levels) -> [(rule, star)]   # 不一致の一覧（§15.3）
```

`Gaes(..., orbitals=["Rb4p=valence"])`、`kkr-gaes run/site/hea --orbital Rb4p=occupied --orbital Bi6s=core`（複数可）。`Judgement` に `orbital_bounds: [min, max]`、`orbital_levels: {"Rb4p": [E, star]}`、`orbital_mismatch: [..]` を追加。`KeyResult.parameters` に `orbitals` を残す。

aiida-akaikkr: `gaes` Dict に `orbitals: ["Rb4p=valence"]`。WorkChain は各 go の後に `levels_from_go` を parser 出力（`results` に成分ごとの core level を加える: 新規 parser 項目 `core_levels`）から取り、history に同じ項目を書く。CLI `submit-gaes --orbital`、MCP `kkr_submit_gaes(orbitals=[...])`。

### 15.5 例（X-Mn-Fe-Co fcc、§13.3 の最終 DOS で判定し直した値、ediff 0.2）

| 指定 | 準位 ε (Ry) | 範囲 | 現在の ewidth | 結果（範囲を選択に使うだけ） | 備考 |
|---|---|---|---|---|---|
| Rb 4p valence | −0.936 `*` | [1.14, —] | 1.2 | `old`（ギャップ [−1.93, −1.18]） | 現状と同じ |
| Rb 4p core | −0.936 `*` | [—, 0.74] | 1.2 | 候補 0.59（ギャップ [−0.65, −0.5]）→ 0.59 で go | 現状は `old` を返す（§15.2 の条件追加で `new`） |
| Se 4s valence | −1.042 `*` | [1.24, —] | 1.414 | `old` | 現状と同じ |
| Se 4s core | −1.042 `*` | [—, 0.84] | 1.414 | 候補 0.76 | 4s の上の低 DOS 区間 |
| Bi 6s valence | −0.954 `*` | [1.15, —] | 1.373 | `old` | 現状と同じ |
| Bi 6s core | −0.954 `*` | [—, 0.75] | 1.373 | 候補無し（eth 0.3）→ `fail`；eth 0.2 なら 0.73 | 6s と valence 帯の間の区間は幅 0.27 |
| Te 5s core | −0.987 `*` | [—, 0.79] | 1.403 | 候補無し → `fail` | 区間幅 0.24 |
| In 4d core | −1.128 `*` | [—, 0.93] | 1.506 | 候補無し → `fail` | 区間幅 0.27〜0.29 |

`core` 指定では、その軌道と valence 帯の底の間に幅 eth のギャップが要る。Bi 6s、Te 5s、In 4d のように幅が 0.24〜0.29 Ry しか無い系では eth を 0.2 に下げないと `fail` になる。この場合の `fail` は仕様どおり（指定を満たす積分路が無い）で、message に範囲と区間幅を出す。

### 15.5.1 実装と実走の記録（2026-09-25、branch `ewidth_select_orbital`）

実装: `pyakaikkr.gaes.orbital`（`OrbitalRule`, `parse_orbital_rules`, `levels_from_go` / `levels_from_records`, `levels_from_tables`, `levels_for_step0`, `bounds_from_rules`, `check_rules`, `initial_ewidth`）、`AkaikkrJob.get_core_levels_by_component`、`Gaes(orbitals=..., ef_assumed=...)`、`kkr-gaes ... --orbital`、`kkr-gaes check --go --orbital`、aiida-akaikkr（parser `results["core_levels"]`、`AkaikkrGaesWorkChain` の `gaes.orbitals`、`submit-gaes --orbital`、`kkr_submit_gaes(orbital=...)`、終了コード 423、`plot --gaes-pk` の判定ごとの範囲と準位）。§15.2 の `old` 判定の範囲条件も入れた（`within_bounds`、Method 1 / 2 とも）。テスト `tests/gaes/test_orbital.py`（fixture `tests/gaes/data/out_go_RbMnFeCo_fcc_ew1.2.log`）。

実走で分かったこと（X-Mn-Fe-Co fcc、akaikkr ビルド、Method 2、ediff 0.2）:

- **core 扱いと valence 扱いで準位が 0.7〜0.9 Ry 違う**。Bi 6s は valence（`*`）なら −0.96 Ry、core なら −1.83 Ry。Rb 4p は −0.94 / −1.68 Ry、Se 4s は −0.97 / −1.05 Ry。ある go の準位だけから範囲を決めると、core 扱いの go の後は範囲が緩み（Bi: max 1.63）、次の go で準位が valence に戻って振動する。そこで **その key で見た準位をすべて記録し、core 指定には最も浅い値、valence 指定には最も深い値** を使う（`Gaes._levels_seen`、WorkChain の `ctx.levels_seen`）。
- 新規ポテンシャルで始める go では、reconf が原子の初期準位（表 B に相当）と ebtm を比べて `*` を決める。表 B が 0.2〜0.5 Ry 深いため、段階 0 を表 B で決めると valence 指定の初期 ewidth が深すぎる（Rb 4p valence: 1.2 → 1.86 で `old`、範囲内なので受理されてしまう）。段階 0 は **表 A（2019 年の収束値）だけ** で決め、表に無い元素は `ewidth_init` のまま始めて最初の判定で直す（`levels_for_step0`）。表に無い元素の core 指定も段階 0 では誤りにしない（`strict=False`）。
- **段階 0 は ewidth_init を深くする方向にしか動かさない**（`initial_ewidth`）。core 指定の上限に合わせて浅くすると（Tl 5d core: 表 A から 0.39、Ga 3d 0.61、In 4d 0.64）積分路が valence 帯の中に入り、SCF が発散する（Tl: rms 1.0、5d が −3.2 Ry へ）か、ギャップの無い DOS になる。ewidth_init（1.2）で 1 回 go を回せば、その DOS から範囲内の候補が取れる（Bi 6s core, eth 0.2: 0.4 → 1.40 → 0.73 の 3 回が、1.2 → 0.7275 の 2 回になった）。
- valence 指定の下限が今の積分路より深いとき、dos の窓を E_F − min_ewidth − eth − ediff まで広げてから判定する（Rb 4s valence: min 2.3 → ewidth_dos 3.73）。
- 結果:

  | 指定 | ewidth の経過 | 範囲の経過 | 結果 |
  |---|---|---|---|
  | Rb 4p valence | 1.2 | [1.14, —] | `finished` 1.2（Rb 4p は `*`） |
  | Bi 6s core | 0.4（表 A から）→ 1.4025 → | [—, 1.63] → [—, 0.754] | `ewidth_fail`（eth 0.3 では 6s と valence 帯の間に幅 0.3 の区間が無い） |
  | Bi 6s core, eth 0.2 | 0.4 → 1.4025 → 0.7275（修正版: 1.2 → 0.7275） | [—, 1.63] → [—, 0.754] → [—, 0.754] | `finished` 0.7275（Bi 6s は core、−0.96 Ry、SCF 収束） |
  | Rb 4p core | 1.2 → 0.5875 | [—, 0.736] → [—, 0.736] | `ewidth_fail`（0.5875 で SCF 発散、rms 0.95 / 500 反復、低 DOS 区間無し。積分路の下端 −0.59 が valence 帯の底 −0.58 に接する。Rb 4p を core にした解はこの系に無い） |
  | Se 4s core | 1.2 → 0.7725 | [—, 0.844] → [—, 0.844] | `finished` 0.7725（Se 4s は core、−1.05 Ry、SCF 収束、ギャップ [−0.89, −0.67]） |

  他の元素（同じ X-Mn-Fe-Co fcc、ewidth_init 1.2、修正版）:

  | 指定 | 準位（1.2 の go、E − E_F） | 範囲 | ewidth の経過 | 結果 |
  |---|---|---|---|---|
  | Pr 5p valence | −1.41（core） | [1.65, —] | 1.2 → 1.7775 | `finished`（5p は `*`、−1.40 Ry、5p 帯の下のギャップ [−2.27, −1.66]） |
  | Ce 5p valence | −1.56（core） | [1.76, —] | 1.2 → 1.9125 | `finished`（5p は `*`、−1.55 Ry） |
  | La 5p core | −1.32 `*` | [—, 1.11] | 1.2 → 0.8125 | `finished`（5p は core、−1.33 Ry、ギャップ [−0.98, −0.80]） |
  | Ba 5p core | −1.06 `*` | [—, 0.86] | 1.2 → 0.6375 | `finished`（5p は core、−1.08 Ry。STEP1 未収束、STEP2 edelt 1e-2 で収束） |
  | In 4d core | −1.11 `*` | [—, 0.91] | 1.2 | `ewidth_fail`（4d と valence 帯の間の区間は幅 0.27〜0.29 < eth） |
  | Tl 5d core | −0.94 `*` | [—, 0.74] | 1.2 | `ewidth_fail`（同上） |
  | Ga 3d core | −1.08 `*` | [—, 0.88] | 1.2 | `ewidth_fail`（同上） |
  | Rb 4s valence | −2.10（core） | [2.30, —] | 1.2 | `ewidth_fail`（窓を 3.73 に広げても [−2.8, −2.3] に DOS < 2e-2 の幅 0.3 の区間が無い） |
  | In 4d valence | 表 A −1.17 → go −1.07 `*` | [1.27, —] | 1.565（表 A から 1.365 + 0.2） | `finished`（4d 帯の下 [−2.24, −1.40]） |
  | Tl 5d valence | −0.94 `*` | [1.15, —] | 1.2 → 1.2925 | `finished`（5d 帯の下 [−2.24, −1.25]） |
  | Ga 3d valence | 表 A −1.17 → go −1.04 `*` | [1.24, —] | 1.565（表 A から） | `finished`（3d 帯の下 [−2.24, −1.37]） |
  | Rb 4s core | −2.10（core） | [—, 1.90] | 1.2 | `finished`（1.2 は範囲内、ギャップ [−1.93, −1.18]） |
  | Se 4s valence | −1.04 `*` | [1.25, —] | 1.2 → 1.4175 | `finished`（4s 帯の下 [−2.24, −1.25]） |
  | Bi 6s valence | −0.96 `*` | [1.16, —] | 1.2 → 1.3725 | `finished`（5d と 6s の間 [−1.42, −1.22]。指定無しの結果と同じ） |

  同じ軌道を valence / core にした対の図: `docs/data/hea_XMnFeCo_fcc_orbital_pairs_dos.png`（Rb 4p、Se 4s、Bi 6s の 3 組。左 valence、右 core。Rb 4p core は SCF 発散で `ewidth_fail`、Bi 6s core は eth 0.2 の走行）。

  図 `docs/data/hea_XMnFeCo_fcc_orbital_rules_dos.png`（17 面、最終判定の DOS、青線 = 指定した準位（破線 `*`）、斜線 = その判定の範囲、赤 = −ewidth）、結果 `docs/data/hea_XMnFeCo_fcc_orbital_rules_summary.json`、スクリプト `tests/gaes/tools/plot_orbital_rules.py`。Rb 4s valence の面（窓 3.73）では −2.3 Ry より下の DOS に窓を広げたことによる偽のピークが並ぶ（§13、ewidth_dos_max の理由）。

  `core` 指定が通るのは、その準位と valence 帯の間に幅 eth のギャップがある系（La、Ba、Se、Bi は eth 0.2）だけで、In / Tl / Ga のように準位が valence 帯に近い系では `fail` になる。`valence` 指定は準位の下にギャップがあれば通る（Pr、Ce）。

  aiida-akaikkr の `AkaikkrGaesWorkChain`（SLURM、specx-akaikkr@mygardenx2-slurm）でも同じ: Se 4s core は pk 3565 で 1.2 → 0.7725 `finished`、Rb 4p core は pk 3539 で 1.2 → 0.5875 `ewidth_fail`（終了コード 420）。図 `~/aiida_work/figures/gaes_orbital/pk3565_gaes01_dos.png`（斜線帯 [—, 0.844]、Se 4s の準位線）。

  最初の版（各 go の準位だけ、表 B も段階 0 に使用）では、Rb 4p core の 2 回目の判定で範囲が [—, 1.48] に緩み（core 扱いの 4p が −1.68 Ry）、Se 4s core は表 B から 0.955 で始めて 0.955 → 0.7275 → 0.7725 と 3 回の go を要した。

### 15.6 実装しないこと

- 軌道の準位を DOS の判定（§4）に混ぜること。ギャップ判定は DOS の値だけで行う（§0.2）。準位は範囲を決めるためだけに使う。
- 準位そのものを specx に指示すること（specx に軌道を core に固定する入力は無い）。
- Committee（§5）への対応は範囲を members の ewidth 列に掛けるだけで済むので、逐次方式の実装後に行う。

### 15.7 実装の順序

1. `AkaikkrJob.get_core_levels_by_component(outfile)`（[core_levels_vs_dos.md](core_levels_vs_dos.md) §5 の未実装項目）。
2. `pyakaikkr.gaes.orbital`（parse / levels / bounds / check）とオフラインテスト（§13.3 の out_go.log を fixture に使う）。
3. `choose_ewidth` / `choose_ewidth2` の `old` 判定に範囲の条件を追加（既存テストの期待値を確認）。
4. `Gaes` と CLI に `orbitals` を通し、Judgement / KeyResult に記録。Rb 4p core（→ 0.59）、Se 4s core（→ 0.76）、Bi 6s core（→ `fail`、eth 0.2 で 0.73）で実走テスト。
5. aiida-akaikkr（parser の `core_levels`、WorkChain、CLI、MCP、docs/gaes_workchain.md）。
