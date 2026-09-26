# DOS / PDOS 図の ewidth 線

`pyakaikkr` の DOS / PDOS 図（`DosPlotter`, `DosEXPlotter`, `PDosPlotter`, `PDosEXPlotter`）は E − E_F = −|ewidth_go| に一点鎖線（赤）を引きます。ewidth_go は **go の ewidth** で、dos の ewidth ではありません。

## なぜ go の ewidth か

- go は [E_F − ewidth_go, E_F] を積分して電荷を自己無撞着にします（cemesh.f）。この線より下の状態は SCF に入っていないので、価電子帯の底が線より右にあることを図で確かめるための線です。
- dos の ewidth はメッシュ [E_F − ref·ewidth_dos, E_F + (1 − ref)·ewidth_dos] を決めるだけです（cemesr.f、ref は akaikkr ビルドで 0.75、akaikkr_cnd ビルドで 0.5）。−ewidth_dos は必ずメッシュの外に落ちるので、これで線を引いても意味がありません。fallback にも使いません。
- 線がメッシュより左（akaikkr_cnd で ewidth_go > 1.0 のとき等）なら、横軸を線まで広げます。

## ewidth_go の決め方（優先順）

1. `make(..., ewidth_go=1.0)` の引数。
2. `DosEXPlotter(directory, outfile, output_directory, go_outfile="out_go.log", read_go_outfile=True)` の `go_outfile` を `directory` から `AkaikkrJob.get_ewidth` で読む（既定で読む）。
3. どちらも無ければ `UserWarning` を出して線を描かない。`read_go_outfile=False` なら警告なしで描かない。

`DosPlotter.make(energy, dos_block, ewidth_go=...)` と `PDosPlotter.make(energy, pdos_block, typeofsites, ewidth_go=...)` は引数だけで、ファイルは読みません。`pyakaikkr.DosPlotter.resolve_ewidth_go(directory, ewidth_go, go_outfile, read_go_outfile)` が上の優先順を実装しています。

## 呼び出し側

- `GoGo.GoDos.postscript` は同じディレクトリの `out_go.log` を `go_outfile` に渡します。テストスクリプト（`tests/akaikkr*/testrun.py`）の dos.png / pdos_*.png はこれで線が付きます。
- `kkr-cmd dos <dir>/out_dos.log -o output` は `<dir>/out_go.log` を読みます（以前は存在しない `make_dos` を呼んでいて動きませんでした）。
- aiida-akaikkr の `aiida_akaikkr.plot` は provenance（`inputs.potential.creator`）から go の ewidth を取る同じ規則です。

テスト: `tests/plot/test_dos_ewidth_line.py`（`tests/testrun/akaikkr/Cu/` の出力があれば go=1.0 / dos=2.0 の取り違えが無いことも確認）。

## 2026-09-26: 描画の共通化

線を引く関数は `pyakaikkr.plot.mark_ewidth_go` に移した（`DosPlotter._mark_ewidth_go` はその薄い皮）。`DosPlotter` / `PDosPlotter` は配列を `pyakaikkr.plot.plot_dos` / `plot_pdos` に渡して描き、aiida-akaikkr の `plot_dos` / `plot_pdos` も同じ関数を使う（`contour_bottom` で go の ewidth を求めて `ewidth_go` に渡す）。規則（go の ewidth、fallback 無し）は変わらない。
