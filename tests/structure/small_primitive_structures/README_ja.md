# 少原子数のバルク原始胞：CIF / POSCAR

7種類のバルク結晶。各フォルダに POSCAR、原始胞CIF、SOURCE.json を収録。
出典に記載された格子定数と座標・構造型から再構成したファイルであり、配布元のCIFをそのまま転載したものではありません。DFT構造最適化は未実施です。

## 採用パラメータと出典

長さの単位は Å。立方晶の a は慣用立方胞の辺長であり、保存した原始胞ベクトルの長さとは異なります。

|物質|原始胞原子数|空間群|採用値|出典|
|---|---:|---|---|---|
|RbCl|2|Fm-3m (225)|a_cubic=6.5810|〔AaltoNDrocksalt〕|
|ZnSe|2|F-43m (216)|a_cubic=5.667|〔AaltoNDzincblende〕、〔AFLOWND216〕|
|InP|2|F-43m (216)|a_cubic=5.8686, 300 K|〔ANLNDlattice〕、〔AFLOWND216〕|
|InAs|2|F-43m (216)|a_cubic=6.0584, 300 K|〔ANLNDlattice〕、〔AFLOWND216〕|
|InSb|2|F-43m (216)|a_cubic=6.4794, 300 K|〔ANLNDlattice〕、〔AFLOWND216〕|
|SnSe2|3|P-3m1 (164)|a=3.811, c=6.13, z_Se=0.25|〔Borges2018SnSe2〕Table I、〔AFLOWND164〕|
|Bi2Se3|5|R-3m (166)|a_hex=4.138, c_hex=28.64; z_Bi=0.399, z_Se=0.206|〔QuantumATK2016Bi2Se3〕|

## 出典の性格

- RbCl、ZnSe：Aalto Universityの構造資料に掲載された値。個別試料の実験CIFではありません。温度を独自に仮定していません。
- InP、InAs、InSb：Argonne National Laboratoryの300 K格子定数表。前の回答の丸めた値より桁数の多い、同表の値を採用しました。
- SnSe2：Borgesらの論文Table Iの0.0(1) GPa行。a=3.811(1), c=6.13(7), z=0.25の中心値を採用。機械的合金化で作製されたナノ構造試料の回折データに基づく周期的結晶モデルで、粒界や有限サイズは含みません。c=6.137など別出典の値を混ぜていません。
- Bi2Se3：QuantumATK公式チュートリアルの構築パラメータを採用。六方晶設定の3a/6cサイトから15原子胞を生成し、5原子原始胞へ変換。チュートリアルの値であり、今回のDFT緩和値ではありません。

## 原子座標と原始胞への変換

- 岩塩型、閃亜鉛鉱型：FCC原始格子ベクトル a/2(0,1,1), a/2(1,0,1), a/2(1,1,0)。原始胞分率座標は第1元素(0,0,0)、第2元素は岩塩型(1/2,1/2,1/2)、閃亜鉛鉱型(1/4,1/4,1/4)。これらの座標は構造型で決まり、自由な内部パラメータはありません。
- SnSe2：Sn=(0,0,0)、Se=(1/3,2/3,z), (2/3,1/3,-z)、z=0.25。三方晶の六方軸表示で、この3原子胞自体が原始胞です。
- Bi2Se3：六方晶設定のSe 3a=(0,0,0)、Bi 6c=(0,0,0.399)、Se 6c=(0,0,0.206)を空間群166で展開して原始胞化。出力された原始胞での座標はSOURCE.jsonに明記。

## CIFの対称性表記

原始胞の形状と全原子座標を確実に保持するため、CIFはP1表記で全サイトを明示しています。これは結晶の対称性がP1であるという意味ではありません。実際の空間群は上表、CIF冒頭コメント、SOURCE.jsonに記載。spglibで座標から再検出して一致を確認しました。

## 検証・利用

CIF / POSCAR双方を再読込し、原子数、検出空間群、格子の計量テンソル、最小像原子間距離を照合。spglibでこれ以上小さい原始胞に縮約されないことも確認。許容差は1e-5 Å。バージョンはmanifest.json参照。

POSCARはそのまま各計算フォルダに配置できます。POTCARはPOSCARの元素順に合わせて用意してください。構造ファイルだけを含み、INCAR、KPOINTS、POTCARは含みません。格子振動等に用いる場合は、採用する計算条件で必要な緩和を実施してください。

build_structures.pyで再作成可能（Python、numpy、ase、spglibが必要）。

## 引用文献・データ出典

NDは公開年不詳を意味します。論文以外のWeb資料にはDOIが付与されていないものがあります。

[AaltoNDrocksalt] Aalto University, Solid State Chemistry, NaCl (rocksalt), Web資料（公開年不詳). DOI: なし／当該ページに記載なし. URL: https://wiki.aalto.fi/spaces/SSC/pages/165132721/NaCl%2Brocksalt

[AaltoNDzincblende] Aalto University, Solid State Chemistry, ZnS (zinc blende), Web資料（公開年不詳). DOI: なし／当該ページに記載なし. URL: https://wiki.aalto.fi/spaces/SSC/pages/165133723/ZnS%2Bzinc%2Bblende

[ANLNDlattice] Argonne National Laboratory, APS 7-ID, Lattice Constants and Crystal Structures of some Semiconductors and Other Materials, Web資料（公開年不詳). DOI: なし／当該ページに記載なし. URL: https://7id.xray.aps.anl.gov/calculators/crystal_lattice_parameters.html

[Borges2018SnSe2] Borges, Z. V., Poffo, C. M., de Lima, J. C., Souza, S. M., Triches, D. M., & de Biasi, R. S., High-pressure angle-dispersive X-ray diffraction study of mechanically alloyed SnSe2, Journal of Applied Physics, 124, 215901 (2018). DOI: 10.1063/1.5053220. URL: https://doi.org/10.1063/1.5053220

[QuantumATK2016Bi2Se3] Synopsys QuantumATK, Bi2Se3 topological insulator, Web資料（version 2016.2). DOI: なし／当該ページに記載なし. URL: https://docs.quantumatk.com/tutorials/topological_insulator_bi2se3/topological_insulator_bi2se3.html

[AFLOWND164] AFLOW Encyclopedia of Crystallographic Prototypes, CdI2 (C6) structure, AB2_hP3_164_a_d-001, Web資料（公開年不詳). DOI: なし／当該ページに記載なし. URL: https://aflow.org/p/AB2_hP3_164_a_d-001/

[AFLOWND216] AFLOW Encyclopedia of Crystallographic Prototypes, Zincblende structure, AB_cF8_216_a_c-001, Web資料（公開年不詳). DOI: なし／当該ページに記載なし. URL: https://aflow.org/p/AB_cF8_216_a_c-001/

