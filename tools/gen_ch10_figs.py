#!/usr/bin/env python3
"""
Generate Chapter 10 figures (ch10_*.png) using NumPy + Matplotlib.

This is optional; the main workflow uses MATLAB chapter scripts.
It is useful in environments where MATLAB is unavailable.

Outputs:
  matlab-signal-processing-book/figs/ch10_*.png
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np


def hann(N: int) -> np.ndarray:
    if N <= 0:
        raise ValueError("N must be positive")
    if N == 1:
        return np.ones(1)
    n = np.arange(N)
    return 0.5 - 0.5 * np.cos(2 * np.pi * n / (N - 1))


def hamming(N: int) -> np.ndarray:
    if N <= 0:
        raise ValueError("N must be positive")
    if N == 1:
        return np.ones(1)
    n = np.arange(N)
    return 0.54 - 0.46 * np.cos(2 * np.pi * n / (N - 1))


def ideal_lpf(fs: float, fc: float, N: int) -> np.ndarray:
    fs = float(fs)
    fc = float(fc)
    if not (0.0 < fc < fs / 2.0):
        raise ValueError("fc must satisfy 0 < fc < fs/2")
    n = np.arange(N)
    alpha = (N - 1) / 2.0
    m = n - alpha
    OmegaC = 2.0 * np.pi * fc / fs
    h = np.zeros(N, dtype=np.float64)
    idx0 = np.isclose(m, 0.0)
    h[idx0] = OmegaC / np.pi
    h[~idx0] = np.sin(OmegaC * m[~idx0]) / (np.pi * m[~idx0])
    return h


def fir_window_design(fs: float, fc: float, N: int, window: str) -> np.ndarray:
    h = ideal_lpf(fs, fc, N)
    if window == "rect":
        w = np.ones(N)
    elif window == "hann":
        w = hann(N)
    elif window == "hamming":
        w = hamming(N)
    else:
        raise ValueError("unknown window")
    h = h * w
    # DC normalization
    h = h / np.sum(h)
    return h


def fir_ls_lpf_design(fs: float, Fp: float, Fs: float, N: int, *, grid_n: int = 4096, Wp: float = 1.0, Ws: float = 80.0) -> np.ndarray:
    """Type-I symmetric FIR LPF by weighted least squares (transition ignored)."""
    fs = float(fs)
    Fp = float(Fp)
    Fs = float(Fs)
    N = int(N)
    if N % 2 != 1:
        raise ValueError("N must be odd (Type-I).")
    if not (0.0 < Fp < Fs < fs / 2.0):
        raise ValueError("Require 0 < Fp < Fs < fs/2.")

    f_grid = np.linspace(0.0, fs / 2.0, int(grid_n))
    omega = 2.0 * np.pi * f_grid / fs
    mask_pb = f_grid <= Fp
    mask_sb = f_grid >= Fs
    mask = mask_pb | mask_sb

    omega_fit = omega[mask]
    d = mask_pb[mask].astype(np.float64)  # 1 in PB, 0 in SB
    w = np.where(d > 0.5, float(Wp), float(Ws))

    M = (N - 1) // 2
    m = np.arange(M + 1)
    A = np.cos(omega_fit[:, None] * m[None, :])
    if A.shape[1] > 1:
        A[:, 1:] *= 2.0

    Aw = A * w[:, None]
    dw = d * w
    c, *_ = np.linalg.lstsq(Aw, dw, rcond=None)

    h = np.zeros(N, dtype=np.float64)
    alpha = M
    h[alpha] = c[0]
    for mm in range(1, M + 1):
        h[alpha + mm] = c[mm]
        h[alpha - mm] = c[mm]

    # DC normalization
    h = h / np.sum(h)
    return h


def freq_response(fs: float, h: np.ndarray, f: np.ndarray) -> np.ndarray:
    """DTFT evaluation on unit circle for FIR: H(e^{jΩ}) = sum h[n] e^{-jΩn}."""
    fs = float(fs)
    h = np.asarray(h, dtype=np.float64).reshape(-1)
    f = np.asarray(f, dtype=np.float64).reshape(-1)
    omega = 2.0 * np.pi * f / fs
    n = np.arange(h.size)
    E = np.exp(-1j * omega[:, None] * n[None, :])
    return E @ h


def group_delay_from_H(fs: float, f: np.ndarray, H: np.ndarray) -> np.ndarray:
    fs = float(fs)
    f = np.asarray(f, dtype=np.float64).reshape(-1)
    omega = 2.0 * np.pi * f / fs
    phi = np.unwrap(np.angle(H))
    return -np.gradient(phi, omega)


def psd_welch_like(x: np.ndarray, fs: float, L: int, hop: int, nfft: int) -> tuple[np.ndarray, np.ndarray]:
    """Very small Welch-like PSD (Hann window, one-sided)."""
    x = np.asarray(x, dtype=np.float64).reshape(-1)
    N = x.size
    if N < L:
        raise ValueError("x shorter than one frame")
    M = (N - L) // hop + 1
    start = np.arange(M) * hop
    idx = start[None, :] + np.arange(L)[:, None]
    X = x[idx]
    w = hann(L)
    Xw = X * w[:, None]
    Xf = np.fft.rfft(Xw, n=nfft, axis=0)
    win_pow = np.sum(w**2)
    P2 = (np.abs(Xf) ** 2) / (fs * win_pow)  # two-sided (rfft gives 0..Nyq)
    P = P2.copy()
    # one-sided correction for power
    if nfft > 1:
        if nfft % 2 == 0:
            if P.shape[0] > 2:
                P[1:-1, :] *= 2
        else:
            if P.shape[0] > 1:
                P[1:, :] *= 2
    Pxx = np.mean(P, axis=1)
    f = np.arange(Pxx.size) * (fs / nfft)
    return f, Pxx


def stft_psd(x: np.ndarray, fs: float, L: int, hop: int, nfft: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=np.float64).reshape(-1)
    N = x.size
    M = (N - L) // hop + 1
    start = np.arange(M) * hop
    idx = start[None, :] + np.arange(L)[:, None]
    X = x[idx]
    w = hann(L)
    # center time
    t = (start + (L - 1) / 2.0) / fs
    Xw = X * w[:, None]
    Xf = np.fft.rfft(Xw, n=nfft, axis=0)
    win_pow = np.sum(w**2)
    P2 = (np.abs(Xf) ** 2) / (fs * win_pow)
    P = P2.copy()
    if nfft > 1:
        if nfft % 2 == 0:
            if P.shape[0] > 2:
                P[1:-1, :] *= 2
        else:
            if P.shape[0] > 1:
                P[1:, :] *= 2
    f = np.arange(P.shape[0]) * (fs / nfft)
    return t, f, P


def savefig(fig, path: Path, dpi: int = 200) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")


def main() -> int:
    import matplotlib.pyplot as plt

    out_dir = Path(__file__).resolve().parents[1] / "figs"

    fs = 8000.0
    Fp = 800.0
    Fs = 1200.0
    Fc = 0.5 * (Fp + Fs)
    N = 101
    n = np.arange(N)
    alpha = (N - 1) / 2.0

    # 1) ideal lpf sinc slice
    h_id = ideal_lpf(fs, Fc, N)
    fig, ax = plt.subplots(1, 1, figsize=(8.0, 3.2))
    ax.plot(n - alpha, h_id, lw=1.2)
    ax.grid(True)
    ax.set_xlabel("n-α")
    ax.set_ylabel("h_d[n]")
    ax.set_title(f"Ideal LPF impulse response slice (sinc): Fc={int(Fc)} Hz, N={N}")
    savefig(fig, out_dir / "ch10_ideal_lpf_sinc.png")
    plt.close(fig)

    # 2) windowed coefficients
    h_rect = fir_window_design(fs, Fc, N, "rect")
    h_hann = fir_window_design(fs, Fc, N, "hann")
    h_hamm = fir_window_design(fs, Fc, N, "hamming")
    fig, ax = plt.subplots(1, 1, figsize=(8.0, 3.2))
    ax.plot(n - alpha, h_rect, lw=1.0, label="rect")
    ax.plot(n - alpha, h_hann, lw=1.0, label="hann")
    ax.plot(n - alpha, h_hamm, lw=1.0, label="hamming")
    ax.grid(True)
    ax.set_xlabel("n-α")
    ax.set_ylabel("h[n]")
    ax.set_title(f"Windowed FIR coefficients (LPF): N={N}, Fc={int(Fc)} Hz (spec: Fp={int(Fp)}, Fs={int(Fs)})")
    ax.legend(loc="best")
    savefig(fig, out_dir / "ch10_coeffs_windows.png")
    plt.close(fig)

    # 3) step response (delay-compensated) for different windows
    Nstep = 600
    x_step = np.ones(Nstep, dtype=np.float64)

    def filt_causal(h: np.ndarray, x: np.ndarray) -> np.ndarray:
        y = np.convolve(x, h, mode="full")[: x.size]
        return y

    y_rect = filt_causal(h_rect, x_step)
    y_hann = filt_causal(h_hann, x_step)
    y_hamm = filt_causal(h_hamm, x_step)

    d_samp = int((N - 1) // 2)
    y_rect_a = y_rect[d_samp:]
    y_hann_a = y_hann[d_samp:]
    y_hamm_a = y_hamm[d_samp:]
    n_a = np.arange(y_rect_a.size)

    fig, ax = plt.subplots(1, 1, figsize=(8.0, 3.6))
    ax.plot(n_a, y_rect_a, lw=1.1, label="rect")
    ax.plot(n_a, y_hann_a, lw=1.1, label="hann")
    ax.plot(n_a, y_hamm_a, lw=1.1, label="hamming")
    ax.axhline(1.0, color="k", ls=":", lw=1.0, label="ideal (1)")
    ax.grid(True)
    ax.set_xlabel("n (delay-compensated)")
    ax.set_ylabel("y[n]")
    ax.set_title(f"Step response (ringing): N={N}, Fc={int(Fc)} Hz (spec: Fp={int(Fp)}, Fs={int(Fs)})")
    ax.set_xlim(0, 220)
    ax.legend(loc="best")
    savefig(fig, out_dir / "ch10_step_response_windows.png")
    plt.close(fig)

    # Frequency axis
    f = np.linspace(0.0, fs / 2.0, 4096)
    Hrect = freq_response(fs, h_rect, f)
    Hhann = freq_response(fs, h_hann, f)
    Hhamm = freq_response(fs, h_hamm, f)

    # 3) window compare magnitude
    fig, ax = plt.subplots(1, 1, figsize=(8.0, 3.6))
    ax.plot(f, 20 * np.log10(np.abs(Hrect) + 1e-20), lw=1.1, label="rect")
    ax.plot(f, 20 * np.log10(np.abs(Hhann) + 1e-20), lw=1.1, label="hann")
    ax.plot(f, 20 * np.log10(np.abs(Hhamm) + 1e-20), lw=1.1, label="hamming")
    ax.grid(True)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Magnitude [dB]")
    ax.set_title(f"Window comparison (LPF): N={N}, Fc={int(Fc)} Hz (spec: Fp={int(Fp)}, Fs={int(Fs)})")
    ax.set_xlim(0, fs / 2)
    ax.set_ylim(-120, 5)
    ax.axvline(Fp, color="k", ls="--")
    ax.axvline(Fs, color="k", ls="--")
    ax.legend(loc="lower left")
    savefig(fig, out_dir / "ch10_window_compare_mag.png")
    plt.close(fig)

    # 4) tap count effect (hann)
    Ns = [31, 63, 127]
    fig, ax = plt.subplots(1, 1, figsize=(8.0, 3.6))
    for Ni in Ns:
        hi = fir_window_design(fs, Fc, Ni, "hann")
        Hi = freq_response(fs, hi, f)
        ax.plot(f, 20 * np.log10(np.abs(Hi) + 1e-20), lw=1.1, label=f"N={Ni}")
    ax.grid(True)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Magnitude [dB]")
    ax.set_title(f"Tap count effect (Hann window): Fc={int(Fc)} Hz (spec: Fp={int(Fp)}, Fs={int(Fs)})")
    ax.set_xlim(0, fs / 2)
    ax.set_ylim(-120, 5)
    ax.axvline(Fp, color="k", ls="--")
    ax.axvline(Fs, color="k", ls="--")
    ax.legend(loc="lower left")
    savefig(fig, out_dir / "ch10_tap_count_effect.png")
    plt.close(fig)

    # 5) linear phase + group delay (hann)
    # Phase / group delay are meaningful mainly in the passband; in deep stopbands
    # |H|≈0 so ∂φ/∂ω becomes numerically unstable and the plot gets unreadable.
    mask_pb = f <= Fp
    f_pb = f[mask_pb]
    H_pb = Hhann[mask_pb]
    phi_pb = np.unwrap(np.angle(H_pb))
    gd_pb = group_delay_from_H(fs, f_pb, H_pb)

    fig, ax = plt.subplots(2, 1, figsize=(8.0, 5.4), sharex=True)
    ax[0].plot(f_pb, phi_pb, lw=1.1)
    ax[0].grid(True)
    ax[0].set_ylabel("Phase [rad]")
    ax[0].set_title("Unwrapped phase of a linear-phase FIR (passband zoom)")

    ax[1].plot(f_pb, gd_pb, lw=1.1)
    ax[1].grid(True)
    ax[1].set_xlabel("Frequency [Hz]")
    ax[1].set_ylabel("Group delay [samples]")
    ax[1].set_title(f"Group delay in passband: (N-1)/2 = {(N-1)/2:.1f} samples")
    ax[1].set_xlim(0, Fp)
    gd0 = (N - 1) / 2.0
    ax[1].set_ylim(gd0 - 6.0, gd0 + 6.0)
    fig.tight_layout()
    savefig(fig, out_dir / "ch10_linear_phase_gd.png")
    plt.close(fig)

    # 6) FIR: direct vs FFT(h)
    Nfft = 4096
    Hfft = np.fft.rfft(h_hann, n=Nfft)
    f_fft = np.arange(Hfft.size) * (fs / Nfft)
    Hdir = freq_response(fs, h_hann, f_fft)

    fig, ax = plt.subplots(1, 1, figsize=(8.0, 3.6))
    ax.plot(f_fft, 20 * np.log10(np.abs(Hdir) + 1e-20), lw=1.2, label="direct eval")
    ax.plot(f_fft, 20 * np.log10(np.abs(Hfft) + 1e-20), lw=1.0, ls="--", label="FFT(h)")
    ax.grid(True)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Magnitude [dB]")
    ax.set_title("FIR: direct evaluation equals FFT(h) (up to numerical error)")
    ax.set_xlim(0, fs / 2)
    ax.set_ylim(-120, 5)
    ax.legend(loc="lower left")
    savefig(fig, out_dir / "ch10_fir_direct_vs_fft.png")
    plt.close(fig)

    # 7) spec attenuation bar (window compare)
    def stopband_atten_db(H: np.ndarray) -> float:
        mask = f >= Fs
        max_db = np.max(20 * np.log10(np.abs(H[mask]) + 1e-20))
        return float(-max_db)

    att = [stopband_atten_db(Hrect), stopband_atten_db(Hhann), stopband_atten_db(Hhamm)]
    fig, ax = plt.subplots(1, 1, figsize=(6.8, 3.2))
    ax.bar(["rect", "hann", "hamming"], att)
    ax.grid(True, axis="y")
    ax.set_ylabel("Stopband attenuation [dB]")
    ax.set_title(f"Measured stopband attenuation (N={N}, Fp={int(Fp)}, Fs={int(Fs)})")
    savefig(fig, out_dir / "ch10_specs_atten_compare.png")
    plt.close(fig)

    # 8) least squares design vs window method
    h_ls = fir_ls_lpf_design(fs, Fp, Fs, N, grid_n=4096, Wp=1.0, Ws=80.0)
    Hls = freq_response(fs, h_ls, f)

    def stopband_atten_db_range(H: np.ndarray, f_axis: np.ndarray, f_start: float) -> float:
        mask = f_axis >= f_start
        max_db = np.max(20 * np.log10(np.abs(H[mask]) + 1e-20))
        return float(-max_db)

    As_win = stopband_atten_db_range(Hhann, f, Fs)
    As_ls = stopband_atten_db_range(Hls, f, Fs)

    fig, ax = plt.subplots(1, 1, figsize=(8.0, 3.6))
    ax.plot(f, 20 * np.log10(np.abs(Hhann) + 1e-20), lw=1.1, label="window (hann)")
    ax.plot(f, 20 * np.log10(np.abs(Hls) + 1e-20), lw=1.1, label="least squares (Ws=80)")
    ax.grid(True)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Magnitude [dB]")
    ax.set_title(f"Window(Hann) vs LS: N={N} (stopband: {As_win:.1f} dB vs {As_ls:.1f} dB)")
    ax.set_xlim(0, fs / 2)
    ax.set_ylim(-120, 5)
    ax.axvline(Fp, color="k", ls="--")
    ax.axvline(Fs, color="k", ls="--")
    ax.legend(loc="lower left")
    savefig(fig, out_dir / "ch10_window_vs_ls_mag.png")
    plt.close(fig)

    # 9) LS weight tradeoff
    Ws_list = [10.0, 30.0, 100.0]
    H_list = []
    for Ws_i in Ws_list:
        h_i = fir_ls_lpf_design(fs, Fp, Fs, N, grid_n=4096, Wp=1.0, Ws=Ws_i)
        H_list.append(freq_response(fs, h_i, f))

    fig, ax = plt.subplots(2, 1, figsize=(8.0, 6.2))
    for H_i, Ws_i in zip(H_list, Ws_list, strict=True):
        ax[0].plot(f, 20 * np.log10(np.abs(H_i) + 1e-20), lw=1.1, label=f"Ws={int(Ws_i)}")
    ax[0].grid(True)
    ax[0].set_ylabel("Magnitude [dB]")
    ax[0].set_title("Passband zoom")
    ax[0].set_xlim(0, Fp)
    ax[0].set_ylim(-1, 1)
    ax[0].legend(loc="best")

    for H_i, Ws_i in zip(H_list, Ws_list, strict=True):
        ax[1].plot(f, 20 * np.log10(np.abs(H_i) + 1e-20), lw=1.1, label=f"Ws={int(Ws_i)}")
    ax[1].grid(True)
    ax[1].set_xlabel("Frequency [Hz]")
    ax[1].set_ylabel("Magnitude [dB]")
    ax[1].set_title("Stopband zoom")
    ax[1].set_xlim(Fs, fs / 2)
    ax[1].set_ylim(-120, 5)

    # Avoid title/xlabel collisions between stacked axes.
    fig.tight_layout()
    savefig(fig, out_dir / "ch10_ls_weight_tradeoff.png")
    plt.close(fig)

    # 8) LPF/HPF/BPF/BSF examples (same N/window)
    Nex = 101
    win = "hann"

    hL = fir_window_design(fs, 1000.0, Nex, win)
    # HPF via spectral inversion: delta - LPF
    n_ex = np.arange(Nex)
    alpha_ex = (Nex - 1) // 2
    d = np.zeros(Nex)
    d[alpha_ex] = 1.0
    hLP_1000 = fir_window_design(fs, 1000.0, Nex, win)
    hH = d - hLP_1000
    # normalize at Nyquist
    g_pi = np.sum(hH * ((-1.0) ** n_ex))
    hH = hH / g_pi

    # BPF/BSF via differences of LPFs
    hLP_1400 = fir_window_design(fs, 1400.0, Nex, win)
    hLP_600 = fir_window_design(fs, 600.0, Nex, win)
    hB = hLP_1400 - hLP_600
    # scale at center frequency
    f0 = 1000.0
    n = np.arange(Nex)
    omega0 = 2.0 * np.pi * f0 / fs
    g0 = np.sum(hB * np.exp(-1j * omega0 * n))
    hB = hB / np.abs(g0)

    hS = d - hB
    hS = hS / np.sum(hS)

    f_ex = np.linspace(0.0, fs / 2.0, 4096)
    HL = freq_response(fs, hL, f_ex)
    HH = freq_response(fs, hH, f_ex)
    HB = freq_response(fs, hB, f_ex)
    HS = freq_response(fs, hS, f_ex)

    fig, ax = plt.subplots(2, 2, figsize=(9.0, 6.2), sharex=True, sharey=True)
    ax[0, 0].plot(f_ex, 20 * np.log10(np.abs(HL) + 1e-20), lw=1.1)
    ax[0, 0].grid(True)
    ax[0, 0].set_title("LPF")

    ax[0, 1].plot(f_ex, 20 * np.log10(np.abs(HH) + 1e-20), lw=1.1)
    ax[0, 1].grid(True)
    ax[0, 1].set_title("HPF")

    ax[1, 0].plot(f_ex, 20 * np.log10(np.abs(HB) + 1e-20), lw=1.1)
    ax[1, 0].grid(True)
    ax[1, 0].set_title("BPF")

    ax[1, 1].plot(f_ex, 20 * np.log10(np.abs(HS) + 1e-20), lw=1.1)
    ax[1, 1].grid(True)
    ax[1, 1].set_title("BSF")

    for a in ax.ravel():
        a.set_xlim(0, fs / 2)
        a.set_ylim(-120, 5)
        a.set_xlabel("Frequency [Hz]")
        a.set_ylabel("Magnitude [dB]")
    fig.suptitle(f"Band type examples (window={win}, N={Nex})")
    savefig(fig, out_dir / "ch10_band_filters_mag.png")
    plt.close(fig)

    # 8) mini project: band extraction by LPF
    T = 2.0
    Nsamp = int(round(T * fs))
    t = np.arange(Nsamp) / fs
    x_wanted = np.sin(2 * np.pi * 300.0 * t) * (t < 1.0) + 0.7 * np.sin(2 * np.pi * 300.0 * t) * (t >= 1.0)
    rng = np.random.default_rng(0)
    x_interf = 0.6 * np.sin(2 * np.pi * 2200.0 * t) + 0.25 * rng.standard_normal(Nsamp)
    x = x_wanted + x_interf

    # FIR filtering (use convolution for portability)
    y = np.convolve(x, h_hann, mode="same")

    fig, ax = plt.subplots(2, 1, figsize=(8.0, 5.0), sharex=True)
    ax[0].plot(t, x, lw=1.0)
    ax[0].grid(True)
    ax[0].set_ylabel("Amplitude")
    ax[0].set_title("Input signal (wanted + interference)")
    ax[1].plot(t, y, lw=1.0)
    ax[1].grid(True)
    ax[1].set_xlabel("Time [s]")
    ax[1].set_ylabel("Amplitude")
    ax[1].set_title("Output after FIR low-pass (window method)")
    savefig(fig, out_dir / "ch10_project_time_before_after.png")
    plt.close(fig)

    # PSD before/after (Welch-like)
    fP, Pxx = psd_welch_like(x, fs, L=512, hop=256, nfft=2048)
    _, Pyy = psd_welch_like(y, fs, L=512, hop=256, nfft=2048)

    fig, ax = plt.subplots(1, 1, figsize=(8.0, 3.6))
    ax.plot(fP, 10 * np.log10(Pxx + 1e-20), lw=1.1, label="input")
    ax.plot(fP, 10 * np.log10(Pyy + 1e-20), lw=1.1, label="output")
    ax.grid(True)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("PSD [dB/Hz]")
    ax.set_title("PSD before/after filtering (Welch-like average)")
    ax.set_xlim(0, fs / 2)
    ax.legend(loc="best")
    savefig(fig, out_dir / "ch10_project_psd_before_after.png")
    plt.close(fig)

    # Spectrogram (PSD)
    tS, fS, Pst = stft_psd(x, fs, L=256, hop=128, nfft=1024)
    _, _, PstY = stft_psd(y, fs, L=256, hop=128, nfft=1024)

    def spec_image(ax, t_axis, f_axis, P, title: str, dr_db: float = 80.0) -> None:
        Z = 10.0 * np.log10(P + 1e-20)
        z_max = np.max(Z)
        Z = np.maximum(Z, z_max - dr_db)
        im = ax.imshow(
            Z,
            origin="lower",
            aspect="auto",
            extent=[t_axis[0], t_axis[-1], f_axis[0], f_axis[-1]],
            cmap="viridis",
        )
        ax.set_ylabel("Frequency [Hz]")
        ax.set_title(title)
        return im

    fig, ax = plt.subplots(2, 1, figsize=(8.0, 6.2), sharex=True, sharey=True)
    im0 = spec_image(ax[0], tS, fS, Pst, "Input spectrogram (PSD)")
    im1 = spec_image(ax[1], tS, fS, PstY, "Output spectrogram (PSD)")
    ax[1].set_xlabel("Time [s]")
    cb = fig.colorbar(im1, ax=ax.ravel().tolist())
    cb.set_label("PSD [dB/Hz]")
    savefig(fig, out_dir / "ch10_project_spectrogram_before_after.png")
    plt.close(fig)

    return 0


if __name__ == "__main__":
    os.environ.setdefault("MPLCONFIGDIR", str(Path.cwd() / "build" / "mplconfig"))
    raise SystemExit(main())


