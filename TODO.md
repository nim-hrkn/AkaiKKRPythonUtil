# TODO（AkaiKKRPythonUtil / pyakaikkr）

- [ ] **GAES の ewidth 決定に PDOS を読む**（2026-09-26、ユーザー指示）。いまは total DOS（原子あたり）だけで判定する。希薄成分の semicore は total DOS では濃度分しか現れず、閾値 1e-3 と同程度なので検知できない。`gap.py` の区間検出を「energy mesh + 任意本数の DOS 曲線（total DOS でも成分ごとの PDOS でも）」を取り全曲線の AND でギャップを決める形に一般化し、`Gaes` が out_dos.log の PDOS（`get_pdos`、濃度は掛かっていない: `spmain.f` の `DOS of component i`、希薄系で要確認）を判定に加えられるようにする。閾値は PDOS ではそのまま使える見込み。軌道指定（§15）とも組み合わせる（指定した成分の PDOS で準位の上下のギャップを見る）。設計方針は docs/ewidth_tuning_scheme.md §0.3。
