# MATLABによる信号処理入門（仮） — サンプルコード

<p align="center">
  <img src="figs/Title.png" width="320" alt="cover (front)" />
</p>

このリポジトリは、書籍『MATLABによる信号処理入門』（仮題）の**配布コード**をまとめたものです。  
本書では「式 → 実装（コード） → 図で確認 → 数値（誤差）で確認」という往復で、信号処理の基礎を腹落ちさせることを狙っています。

> 注意：本文原稿（LaTeX/PDF）はこのリポジトリに含めません（配布コード用）。

## 特徴
- **Signal Processing Toolbox なし**で動くことを基本方針にしています（`fft`/`ifft`/`filter` などMATLAB標準関数中心）。
- 各章のスクリプトは原則として **Run All で図まで再現**できるようにしています。
- 再利用コードは `src/sp_*.m` として整理し、章スクリプトから呼び出します。
- 図の出力は `exportgraphics` を基本にし、`figs/` に保存します（生成図は基本的にGit管理しません）。

## 使い方（最短）
1. このリポジトリをクローン
   - HTTPS: `git clone https://github.com/isobe-tech/MATLAB-signal-processing-book.git`
   - SSH: `git clone git@github.com:isobe-tech/MATLAB-signal-processing-book.git`
2. MATLABでリポジトリ直下を Current Folder に設定
3. 章スクリプトを実行（例：第1章）
   - `run(fullfile("chapters","ch01_environment.m"))`

各章スクリプトは内部で `src/` を `addpath` し、図の出力先 `figs/` が無ければ作成します。

## ディレクトリ構成
- `chapters/`：章ごとの実験スクリプト（`.m`、図の生成まで含む）
- `src/`：共通関数（`sp_` 接頭辞）
- `project/`：総合プロジェクト
  - `project/ch10_final_project.m`：最終プロジェクト（章の総まとめ）
- `figs/`：表紙画像＋（ローカルで生成される）図の保存先
  - GitHubに含めているのは `figs/Title.png` と `figs/Back.png` のみです

## 推奨MATLABバージョン
- **R2017a以降**（スクリプト内で string（`"..."`）を使用しているため）

## 表紙
| 表紙 | 裏表紙 |
| --- | --- |
| <img src="figs/Title.png" width="240" alt="cover front" /> | <img src="figs/Back.png" width="240" alt="cover back" /> |
