# Understanding Signal Processing with MATLAB
# MATLABで理解する信号処理 ― 周波数解析とフィルタのはじめ方

Companion code for the book *Understanding Signal Processing with MATLAB --- Introduction to Frequency Analysis and Filters*.

本書籍『MATLABで理解する信号処理 ― 周波数解析とフィルタのはじめ方』の配布コードです。

The book builds understanding through a tight loop of **equation → implementation → visualization → numerical verification**.

> **Note:** The book manuscript (LaTeX/PDF) is not included in this repository.

## Features

- **No Signal Processing Toolbox required** --- uses only standard MATLAB functions (`fft`, `ifft`, `filter`, `conv`, etc.).
- Each chapter script reproduces **all figures via Run All**.
- Reusable functions are organized as `src/sp_*.m` and called from chapter scripts.
- Figures are exported with `exportgraphics` to `figs/`.

## Quick Start

1. Clone this repository:
   ```
   git clone https://github.com/isobe-tech/MATLAB-signal-processing-book.git
   ```
2. Set the repository root as MATLAB's **Current Folder**.
3. Run a chapter script (e.g., Chapter 1):
   ```matlab
   run(fullfile("chapters", "ch01_environment.m"))
   ```

Each script automatically adds `src/` to the path and creates `figs/` if needed.

## Directory Structure

| Path | Contents |
|------|----------|
| `chapters/` | Per-chapter experiment scripts |
| `src/` | Shared utility functions (`sp_` prefix, 30 files) |
| `project/` | Capstone project (`ch10_final_project.m`) |
| `figs/` | Figure output (only cover images are tracked in Git) |

## Chapter Guide

| Script | Topic |
|--------|-------|
| `ch01_environment.m` | MATLAB environment, vectors, complex numbers |
| `ch02_sampling_quantization.m` | Sampling theorem, quantization, dB |
| `ch03_complex_phase.m` | Complex exponentials, phase, inner products |
| `ch04_lti_convolution.m` | LTI systems, convolution, Toeplitz matrix |
| `ch05_freq_response_dtft.m` | Frequency response, DTFT |
| `ch06_dft_fft.m` | DFT, FFT, Parseval's theorem |
| `ch07_window_psd.m` | Window functions, PSD estimation |
| `ch08_stft.m` | Short-time Fourier transform, spectrogram |
| `ch09_z_iir.m` | Z-transform, poles/zeros, IIR filters |
| `ch10_fir_design.m` | FIR filter design (window method), capstone project |

## Shared Functions (`src/`)

| Function | Purpose |
|----------|---------|
| `sp_dft.m` / `sp_dft_matrix.m` | DFT from definition; DFT matrix |
| `sp_convolve.m` / `sp_convmtx_toeplitz.m` | Convolution; Toeplitz convolution matrix |
| `sp_fft_axis.m` / `sp_time_axis.m` | Frequency / time axis generation |
| `sp_freq_response.m` / `sp_group_delay.m` | Frequency response; group delay |
| `sp_stft.m` / `sp_psd_avg.m` | STFT; averaged PSD (Welch-style) |
| `sp_fir_window_design.m` / `sp_fir_ls_lpf.m` | FIR design: window method / least squares |
| `sp_window_rect.m` / `sp_window_hann.m` / `sp_window_hamming.m` | Window functions |
| `sp_iir1.m` / `sp_iir2.m` | First / second-order IIR filters |
| `sp_zplane_plot.m` | Pole-zero plot |
| `sp_book_style.m` | Consistent figure styling |

## Requirements

**MATLAB R2017a or later** (uses double-quote string syntax `"..."`).

