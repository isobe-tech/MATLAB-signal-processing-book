# 『MATLABによる信号処理入門（ツールボックス無し・基礎徹底）』執筆用

## MATLAB（配布コード）
- `chapters/`：章ごとの実験スクリプト（`.m`）
- `src/`：再利用する関数（`sp_` 接頭辞）
- `tests/`：簡易テスト（`assert`）
- `data/`：サンプル波形など
- `figs/`：図の出力先（自動生成）

実行例（MATLAB上で）：
- 章スクリプト：`run(fullfile("chapters","ch01_environment.m"))`
- テスト：`run(fullfile("tests","test_ch01_environment.m"))`

## 本文（LaTeX/PDF）
本リポジトリは配布コード用で、本文原稿（LaTeX/PDF）は含めない。

## project/
- `project/`：総合プロジェクト（最終プロジェクト）
