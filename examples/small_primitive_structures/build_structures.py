from pathlib import Path
import json, zipfile, shutil
import numpy as np
import ase, spglib
from ase import Atoms
from ase.build import bulk
from ase.spacegroup import crystal
from ase.io import read, write

root = Path(__file__).resolve().parent
out = root / 'small_primitive_structures'
out.mkdir(exist_ok=True)
sources = {
 'AaltoNDrocksalt': {'title':'NaCl (rocksalt)', 'author':'Aalto University, Solid State Chemistry', 'url':'https://wiki.aalto.fi/spaces/SSC/pages/165132721/NaCl%2Brocksalt', 'detail':'RbCl: cubic a = 6.5810 angstrom; rocksalt structure.'},
 'AaltoNDzincblende': {'title':'ZnS (zinc blende)', 'author':'Aalto University, Solid State Chemistry', 'url':'https://wiki.aalto.fi/spaces/SSC/pages/165133723/ZnS%2Bzinc%2Bblende', 'detail':'ZnSe: cubic a = 5.667 angstrom; zincblende structure.'},
 'ANLNDlattice': {'title':'Lattice Constants and Crystal Structures of some Semiconductors and Other Materials', 'author':'Argonne National Laboratory, APS 7-ID', 'url':'https://7id.xray.aps.anl.gov/calculators/crystal_lattice_parameters.html', 'detail':'300 K cubic constants: InP 5.8686, InAs 6.0584, InSb 6.4794 angstrom.'},
 'Borges2018SnSe2': {'title':'High-pressure angle-dispersive X-ray diffraction study of mechanically alloyed SnSe2', 'author':'Borges, Z. V., Poffo, C. M., de Lima, J. C., Souza, S. M., Triches, D. M., & de Biasi, R. S.', 'url':'https://doi.org/10.1063/1.5053220', 'doi':'10.1063/1.5053220', 'journal':'Journal of Applied Physics, 124, 215901 (2018)', 'detail':'Table I, ambient-pressure row 0.0(1) GPa: a=3.811(1), c=6.13(7) angstrom, z(Se)=0.25. Central values used. Mechanically alloyed nanostructured sample.', 'access_url':'https://www.researchgate.net/publication/329375244_High-pressure_angle-dispersive_X-ray_diffraction_study_of_mechanically_alloyed_SnSe2'},
 'QuantumATK2016Bi2Se3': {'title':'Bi2Se3 topological insulator', 'author':'Synopsys QuantumATK', 'url':'https://docs.quantumatk.com/tutorials/topological_insulator_bi2se3/topological_insulator_bi2se3.html', 'detail':'Tutorial version 2016.2, Build the Bi2Se3 crystal: hexagonal a=4.138, c=28.64 angstrom; Se 3a (0,0,0), Bi 6c z=0.399, Se 6c z=0.206. Reconstructed hexagonal cell then converted to primitive cell.'},
 'AFLOWND164': {'title':'CdI2 (C6) structure, AB2_hP3_164_a_d-001', 'author':'AFLOW Encyclopedia of Crystallographic Prototypes', 'url':'https://aflow.org/p/AB2_hP3_164_a_d-001/', 'detail':'Sn at 1a (0,0,0), Se at 2d (1/3,2/3,z), (2/3,1/3,-z).'},
 'AFLOWND216': {'title':'Zincblende structure, AB_cF8_216_a_c-001', 'author':'AFLOW Encyclopedia of Crystallographic Prototypes', 'url':'https://aflow.org/p/AB_cF8_216_a_c-001/', 'detail':'Zincblende conventional-cell basis at (0,0,0), (1/4,1/4,1/4), with F centering.'}
}
specs = []
for name, a, sg, sid in [('RbCl',6.5810,225,'AaltoNDrocksalt'),('ZnSe',5.667,216,'AaltoNDzincblende'),('InP',5.8686,216,'ANLNDlattice'),('InAs',6.0584,216,'ANLNDlattice'),('InSb',6.4794,216,'ANLNDlattice')]:
    phase = 'rocksalt' if sg==225 else 'zincblende'
    at = bulk(name, phase, a=a)
    specs.append((name, at, sg, [sid]+(['AFLOWND216'] if sg==216 else []), {'a_cubic_angstrom':a},2))
at = crystal(['Sn','Se'],basis=[(0,0,0),(1/3,2/3,0.25)],spacegroup=164,cellpar=[3.811,3.811,6.13,90,90,120],primitive_cell=True)
specs.append(('SnSe2',at,164,['Borges2018SnSe2','AFLOWND164'],{'a_hex_angstrom':3.811,'c_hex_angstrom':6.13,'z_Se':0.25},3))
at = crystal(['Bi','Se','Se'],basis=[(0,0,0.399),(0,0,0),(0,0,0.206)],spacegroup=166,cellpar=[4.138,4.138,28.64,90,90,120],primitive_cell=True)
specs.append(('Bi2Se3',at,166,['QuantumATK2016Bi2Se3'],{'a_hex_angstrom':4.138,'c_hex_angstrom':28.64,'z_Bi_hex':0.399,'z_Se_hex':0.206},5))
records=[]
def dataset(at):
    return spglib.get_symmetry_dataset((at.cell.array,at.get_scaled_positions(),at.numbers),symprec=1e-5)
for name, at, sg, ids, params, n in specs:
    # Group species in formula order for POSCAR / POTCAR consistency.
    order=list(dict.fromkeys(at.get_chemical_symbols()))
    at=at[sorted(range(len(at)),key=lambda i:order.index(at[i].symbol))]
    at.info.clear()
    at.wrap()
    ds=dataset(at)
    assert len(at)==n and ds.number==sg, (name,len(at),ds.number)
    assert len(spglib.find_primitive((at.cell.array,at.get_scaled_positions(),at.numbers),symprec=1e-5)[2])==n
    folder=out/name
    folder.mkdir(exist_ok=True)
    pos=folder/'POSCAR'
    cif=folder/f'{name}_primitive.cif'
    write(pos,at,format='vasp',direct=True,vasp5=True,sort=False)
    lines=pos.read_text().splitlines()
    lines[0]=f'{name} primitive; source={ids[0]}; see SOURCE.json'
    pos.write_text('\n'.join(lines)+'\n')
    write(cif,at,format='cif')
    comments=[f'# {name} primitive; actual space group {ds.international} ({sg})', '# Explicit full-cell sites; CIF is serialized in P1; no symmetry reduction of sites.', '# Reconstructed from cited parameters; not DFT-relaxed.']
    for sid in ids:
        comments += [f'# Source [{sid}]: {sources[sid]["url"]}',f'# {sources[sid]["detail"]}']
    cif.write_text('\n'.join(comments)+'\n'+cif.read_text())
    for path, fmt in [(pos,'vasp'),(cif,'cif')]:
        back=read(path,format=fmt)
        assert len(back)==n and dataset(back).number==sg
        assert np.allclose(back.cell.array@back.cell.array.T,at.cell.array@at.cell.array.T,atol=1e-8)
        assert np.allclose(back.get_all_distances(mic=True),at.get_all_distances(mic=True),atol=1e-7)
    record={'material':name,'n_atoms':n,'spacegroup_number':sg,'spacegroup_symbol':ds.international,'source_ids':ids,'source_parameters':params,'source_records':[sources[i] for i in ids],'construction':'Reconstructed from cited parameters; no energy/force optimization performed.','cell_vectors_angstrom':at.cell.array.tolist(),'species':at.get_chemical_symbols(),'fractional_positions':at.get_scaled_positions().tolist(),'verification':'CIF and POSCAR reread: atom count, symmetry, metric tensor, MIC distances passed; spglib primitive count verified.'}
    (folder/'SOURCE.json').write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n')
    records.append(record)
    print(name,n,ds.international,sg)
(out/'manifest.json').write_text(json.dumps({'ase_version':ase.__version__,'spglib_version':spglib.__version__,'symprec_angstrom':1e-5,'structures':records},indent=2,ensure_ascii=False)+'\n')
readme='''# 少原子数のバルク原始胞：CIF / POSCAR

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

'''
for sid,s in sources.items():
    if 'doi' in s:
        readme+=f'[{sid}] {s["author"]}, {s["title"]}, {s["journal"]}. DOI: {s["doi"]}. URL: {s["url"]}\n\n'
    else:
        readme+=f'[{sid}] {s["author"]}, {s["title"]}, Web資料（'+('version 2016.2' if sid=='QuantumATK2016Bi2Se3' else '公開年不詳')+f'). DOI: なし／当該ページに記載なし. URL: {s["url"]}\n\n'
(out/'README_ja.md').write_text(readme)
shutil.copy2(__file__,out/'build_structures.py')
with zipfile.ZipFile(root/'small_primitive_structures.zip','w',zipfile.ZIP_DEFLATED) as z:
    for p in sorted(out.rglob('*')):
        if p.is_file():z.write(p,p.relative_to(root))
print(root/'small_primitive_structures.zip')
