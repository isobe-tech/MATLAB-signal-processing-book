# chapters/

章ごとの実験スクリプト（`.m`）を置く。
本文の図は，ここから再現できる形にする。

テンプレ：`chapters/_template.m`

## 図の再生成（MATLABで統一したい場合）

本書の本文で `\includegraphics{figs/...}` により埋め込まれるPNGは，**章スクリプト（`chapters/*.m`）を実行すると `figs/` に書き出される**。
（`tools/gen_ch08_figs.py` 等のPython版は「MATLABが無い環境でもPNGを作れる」ための補助ルートで，MATLAB生成に統一する場合は不要。）

### 0) 準備
- MATLABを起動し，Current Folder を `matlab-signal-processing-book/` にする
- 以降はコマンドウィンドウで `run(...)` するか，各 `.m` を開いて「Run」を押す
- 各章スクリプトは冒頭で `src/` を `addpath` し，出力先 `figDir=figs/` を自動で作る

### 1) 本文で使う図（第8〜10章）
以下の3つを実行すれば，本文で参照しているPNG（第8〜10章）がすべて `figs/` に生成される。

```matlab
run(fullfile("chapters","ch08_stft.m"));
run(fullfile("chapters","ch09_z_iir.m"));
run(fullfile("chapters","ch10_fir_design.m"));
```

#### 生成されるファイル名（抜粋）
- 第8章（`tex/08_stft.tex`）：
  - `figs/ch08_signals.png`
  - `figs/ch08_chirp_spectrogram.png`
  - `figs/ch08_hop_spectrogram.png`
  - `figs/ch08_tradeoff_window.png`
  - `figs/ch08_hop_effect.png`
  - `figs/ch08_stft_avg_vs_psd.png`
  - `figs/ch08_cola_weight.png`
- 第9章（`tex/09_z_iir.tex`）：
  - `figs/ch09_iir1_impulse_stability.png`
  - `figs/ch09_freqresp_direct_vs_fft.png`
  - `figs/ch09_pole_radius_ringdown.png`
  - `figs/ch09_resonator_pz_mag.png`
  - `figs/ch09_notch_pz_mag.png`
  - `figs/ch09_allpass_gd.png`
  - `figs/ch09_structure_numeric_error.png`
  - `figs/ch09_coeff_quantization_effect.png`
- 第10章（`tex/10_fir_design.tex`）：
  - `figs/ch10_ideal_lpf_sinc.png`
  - `figs/ch10_coeffs_windows.png`
  - `figs/ch10_step_response_windows.png`
  - `figs/ch10_window_compare_mag.png`
  - `figs/ch10_tap_count_effect.png`
  - `figs/ch10_linear_phase_gd.png`
  - `figs/ch10_fir_direct_vs_fft.png`
  - `figs/ch10_specs_atten_compare.png`
  - `figs/ch10_window_vs_ls_mag.png`
  - `figs/ch10_ls_weight_tradeoff.png`
  - `figs/ch10_band_filters_mag.png`
  - `figs/ch10_project_time_before_after.png`
  - `figs/ch10_project_psd_before_after.png`
  - `figs/ch10_project_spectrogram_before_after.png`

### 2) LaTeXへの反映
PNGを更新したら，ターミナルで以下を実行してPDFへ反映する。

```sh
make
```

### 3) 追加で全章を回したい場合
章1〜7のスクリプトも同様に `run(fullfile("chapters","chXX_....m"))` で実行できる（図が `figs/` に出るものもある）。
