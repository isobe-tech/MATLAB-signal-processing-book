#!/usr/bin/env python3
"""
Generate Chapter 9 figures (ch09_*.png) using NumPy + Matplotlib.

This is optional; the main book workflow uses MATLAB chapter scripts.
It is useful in environments where MATLAB is unavailable, and to ensure
LaTeX can embed the figures.

Outputs:
  matlab-signal-processing-book/figs/ch09_*.png
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np


def iir1_sp(x: np.ndarray, a: float, *, dtype=np.float64) -> np.ndarray:
    """Match sp_iir1: y[n]=(1-a)x[n]+a y[n-1], y[-1]=0."""
    x = np.asarray(x, dtype=dtype).reshape(-1)
    y = np.zeros_like(x)
    yp = dtype(0)
    one = dtype(1)
    a = dtype(a)
    for n in range(x.size):
        yc = (one - a) * x[n] + a * yp
        y[n] = yc
        yp = yc
    return y


def iir2_df1(x: np.ndarray, b: np.ndarray, a: np.ndarray, *, dtype=np.float64) -> np.ndarray:
    x = np.asarray(x, dtype=dtype).reshape(-1)
    b = np.asarray(b, dtype=dtype).reshape(-1)
    a = np.asarray(a, dtype=dtype).reshape(-1)
    assert b.size == 3 and a.size == 3 and a[0] != 0

    b0, b1, b2 = b
    a0, a1, a2 = a
    y = np.zeros_like(x)
    x1 = dtype(0)
    x2 = dtype(0)
    y1 = dtype(0)
    y2 = dtype(0)
    for n in range(x.size):
        yc = (b0 * x[n] + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0
        y[n] = yc
        x2, x1 = x1, x[n]
        y2, y1 = y1, yc
    return y


def iir2_df2(x: np.ndarray, b: np.ndarray, a: np.ndarray, *, dtype=np.float64) -> np.ndarray:
    x = np.asarray(x, dtype=dtype).reshape(-1)
    b = np.asarray(b, dtype=dtype).reshape(-1)
    a = np.asarray(a, dtype=dtype).reshape(-1)
    assert b.size == 3 and a.size == 3 and a[0] != 0

    b0, b1, b2 = b
    a0, a1, a2 = a
    y = np.zeros_like(x)
    w1 = dtype(0)
    w2 = dtype(0)
    for n in range(x.size):
        w0 = (x[n] - a1 * w1 - a2 * w2) / a0
        yc = b0 * w0 + b1 * w1 + b2 * w2
        y[n] = yc
        w2, w1 = w1, w0
    return y


def freq_response(fs: float, b: np.ndarray, a: np.ndarray, *, n: int = 2000) -> tuple[np.ndarray, np.ndarray]:
    fs = float(fs)
    b = np.asarray(b, dtype=np.complex128).reshape(-1)
    a = np.asarray(a, dtype=np.complex128).reshape(-1)
    f = np.linspace(0.0, fs / 2.0, n)
    omega = 2.0 * np.pi * f / fs
    mb = np.arange(b.size)
    ma = np.arange(a.size)
    B = np.exp(-1j * omega[:, None] * mb[None, :]) @ b
    A = np.exp(-1j * omega[:, None] * ma[None, :]) @ a
    H = B / A
    return f, H


def group_delay_from_H(fs: float, f: np.ndarray, H: np.ndarray) -> np.ndarray:
    fs = float(fs)
    f = np.asarray(f, dtype=np.float64).reshape(-1)
    omega = 2.0 * np.pi * f / fs
    phi = np.unwrap(np.angle(H))
    # d/domega
    gd = -np.gradient(phi, omega)
    return gd


def quantize_uniform(x: np.ndarray, bits: int, full_scale: float) -> np.ndarray:
    """Mimic sp_quantize_uniform (round + saturate)."""
    x = np.asarray(x, dtype=np.float64)
    bits = int(bits)
    full_scale = float(full_scale)
    step = (2.0 * full_scale) / (2.0**bits)
    code = np.round(x / step)
    code_min = -(2 ** (bits - 1))
    code_max = +(2 ** (bits - 1))
    code = np.clip(code, code_min, code_max)
    return code * step


def savefig(fig, path: Path, dpi: int = 200) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")


def plot_zplane(ax, b: np.ndarray, a: np.ndarray, *, title: str = "", axis_lim: float = 1.4) -> None:
    b = np.asarray(b, dtype=np.float64).reshape(-1)
    a = np.asarray(a, dtype=np.float64).reshape(-1)
    z = np.roots(b) if b.size > 1 else np.array([], dtype=np.complex128)
    p = np.roots(a) if a.size > 1 else np.array([], dtype=np.complex128)

    th = np.linspace(0, 2 * np.pi, 361)
    ax.plot(np.cos(th), np.sin(th), "k-", lw=1.0, label="unit circle")
    if z.size:
        ax.plot(z.real, z.imag, "bo", mfc="none", mew=1.2, label="zeros")
    if p.size:
        ax.plot(p.real, p.imag, "rx", mew=1.2, label="poles")
    ax.axhline(0, color="k", ls=":", lw=1.0)
    ax.axvline(0, color="k", ls=":", lw=1.0)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-axis_lim, axis_lim)
    ax.set_ylim(-axis_lim, axis_lim)
    ax.grid(True)
    ax.set_xlabel("Re{z}")
    ax.set_ylabel("Im{z}")
    ax.set_title(title)


def main() -> int:
    import matplotlib.pyplot as plt

    out_dir = Path(__file__).resolve().parents[1] / "figs"
    fs = 8000.0

    # 1) sp_iir1 impulse: stable vs unstable
    N = 200
    x = np.zeros(N)
    x[0] = 1.0
    a_st = 0.9
    a_un = 1.01
    y_st = iir1_sp(x, a_st)
    y_un = iir1_sp(x, a_un)

    fig, ax = plt.subplots(2, 1, figsize=(8.0, 5.0), sharex=True)
    n = np.arange(N)
    ax[0].stem(n, y_st, basefmt=" ")
    ax[0].grid(True)
    ax[0].set_ylabel("h[n]")
    ax[0].set_title(f"1st-order IIR impulse response (stable): a={a_st:.3f}")
    ax[1].stem(n, y_un, basefmt=" ")
    ax[1].grid(True)
    ax[1].set_xlabel("n")
    ax[1].set_ylabel("h[n]")
    ax[1].set_title(f"1st-order IIR impulse response (unstable): a={a_un:.3f}")
    savefig(fig, out_dir / "ch09_iir1_impulse_stability.png")
    plt.close(fig)

    # 2) ring-down vs pole radius
    theta = 0.20 * np.pi
    rs = [0.8, 0.95, 0.99]
    N = 400
    x = np.zeros(N)
    x[0] = 1.0

    fig, ax = plt.subplots(len(rs), 1, figsize=(8.0, 6.2), sharex=True)
    n = np.arange(N)
    for i, r in enumerate(rs):
        a = np.array([1.0, -2.0 * r * np.cos(theta), r**2])
        b = np.array([1.0, 0.0, 0.0])
        y = iir2_df1(x, b, a)
        ax[i].plot(n, y, lw=1.2)
        ax[i].grid(True)
        ax[i].set_ylabel("h[n]")
        ax[i].set_title(f"2nd-order ring-down: r={r:.2f}, theta={theta/np.pi:.2f} pi")
    ax[-1].set_xlabel("n")
    savefig(fig, out_dir / "ch09_pole_radius_ringdown.png")
    plt.close(fig)

    # 3) resonator: z-plane + magnitude response
    r = 0.98
    a = np.array([1.0, -2.0 * r * np.cos(theta), r**2])
    b = np.array([1.0, 0.0, 0.0])
    f, H = freq_response(fs, b, a, n=2000)

    fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.2))
    plot_zplane(ax[0], b, a, title="Resonator (poles/zeros)", axis_lim=1.4)
    ax[0].legend(loc="best")
    ax[1].plot(f, 20 * np.log10(np.abs(H) + 1e-20), lw=1.2)
    ax[1].grid(True)
    ax[1].set_xlabel("Frequency [Hz]")
    ax[1].set_ylabel("Magnitude [dB]")
    ax[1].set_title("Magnitude response |H(e^{jΩ})|")
    ax[1].set_xlim(0, fs / 2)
    savefig(fig, out_dir / "ch09_resonator_pz_mag.png")
    plt.close(fig)

    # 4) notch: zeros on unit circle
    theta0 = 0.25 * np.pi
    r0 = 0.98
    b = np.array([1.0, -2.0 * np.cos(theta0), 1.0])
    a = np.array([1.0, -2.0 * r0 * np.cos(theta0), r0**2])
    f, H = freq_response(fs, b, a, n=2000)

    fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.2))
    plot_zplane(ax[0], b, a, title="Notch (poles/zeros)", axis_lim=1.4)
    ax[0].legend(loc="best")
    ax[1].plot(f, 20 * np.log10(np.abs(H) + 1e-20), lw=1.2)
    ax[1].grid(True)
    ax[1].set_xlabel("Frequency [Hz]")
    ax[1].set_ylabel("Magnitude [dB]")
    ax[1].set_title("Notch magnitude response")
    ax[1].set_xlim(0, fs / 2)
    savefig(fig, out_dir / "ch09_notch_pz_mag.png")
    plt.close(fig)

    # 5) all-pass: flat magnitude + group delay
    theta1 = 0.30 * np.pi
    r1 = 0.90
    a = np.array([1.0, -2.0 * r1 * np.cos(theta1), r1**2])
    b = np.array([a[2], a[1], a[0]])  # flip(a)
    f, H = freq_response(fs, b, a, n=2000)
    gd = group_delay_from_H(fs, f, H)

    fig, ax = plt.subplots(2, 1, figsize=(8.0, 6.0), sharex=True)
    ax[0].plot(f, 20 * np.log10(np.abs(H) + 1e-20), lw=1.2)
    ax[0].grid(True)
    ax[0].set_ylabel("Magnitude [dB]")
    ax[0].set_title("All-pass magnitude (should be ~0 dB)")
    # Force a readable scale: the true curve is ~0 dB and numerical noise can
    # otherwise trigger scientific-notation offsets (e.g. 1e-14).
    ax[0].set_ylim(-0.2, 0.2)
    ax[0].ticklabel_format(axis="y", style="plain", useOffset=False)
    ax[1].plot(f, gd, lw=1.2)
    ax[1].grid(True)
    ax[1].set_xlabel("Frequency [Hz]")
    ax[1].set_ylabel("Group delay [samples]")
    ax[1].set_title("All-pass group delay")
    ax[1].set_xlim(0, fs / 2)
    savefig(fig, out_dir / "ch09_allpass_gd.png")
    plt.close(fig)

    # 6) coefficient quantization: pole/zero shift + mag response change
    thetaQ = 0.25 * np.pi
    rQ = 0.995
    b = np.array([1.0, -2.0 * np.cos(thetaQ), 1.0])
    a = np.array([1.0, -2.0 * rQ * np.cos(thetaQ), rQ**2])
    bits = 10
    full_scale = 2.0
    bq = quantize_uniform(b, bits, full_scale)
    aq = a.copy()
    aq[1:] = quantize_uniform(a[1:], bits, full_scale)

    z = np.roots(b)
    p = np.roots(a)
    zq = np.roots(bq)
    pq = np.roots(aq)

    f, H = freq_response(fs, b, a, n=2000)
    _, Hq = freq_response(fs, bq, aq, n=2000)

    fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.2))
    # z-plane overlay
    th = np.linspace(0, 2 * np.pi, 361)
    ax[0].plot(np.cos(th), np.sin(th), "k-", lw=1.0, label="unit circle")
    ax[0].plot(z.real, z.imag, "bo", mfc="none", mew=1.2, label="zeros (orig)")
    ax[0].plot(p.real, p.imag, "bx", mew=1.2, label="poles (orig)")
    ax[0].plot(zq.real, zq.imag, "ro", mfc="none", mew=1.2, label="zeros (quant)")
    ax[0].plot(pq.real, pq.imag, "rx", mew=1.2, label="poles (quant)")
    ax[0].axhline(0, color="k", ls=":", lw=1.0)
    ax[0].axvline(0, color="k", ls=":", lw=1.0)
    ax[0].set_aspect("equal", adjustable="box")
    ax[0].set_xlim(-1.4, 1.4)
    ax[0].set_ylim(-1.4, 1.4)
    ax[0].grid(True)
    ax[0].set_xlabel("Re{z}")
    ax[0].set_ylabel("Im{z}")
    ax[0].set_title(f"Pole/zero shift by quantization ({bits} bits)")
    ax[0].legend(loc="best")

    ax[1].plot(f, 20 * np.log10(np.abs(H) + 1e-20), lw=1.2, label="original")
    ax[1].plot(f, 20 * np.log10(np.abs(Hq) + 1e-20), lw=1.2, label="quantized")
    ax[1].grid(True)
    ax[1].set_xlabel("Frequency [Hz]")
    ax[1].set_ylabel("Magnitude [dB]")
    ax[1].set_title("Magnitude response: original vs quantized")
    ax[1].set_xlim(0, fs / 2)
    ax[1].legend(loc="best")
    savefig(fig, out_dir / "ch09_coeff_quantization_effect.png")
    plt.close(fig)

    # 7) structure + single precision: df1 vs df2 (error vs double)
    rS = 0.9995
    thetaS = 0.30 * np.pi
    b = np.array([1.0, 0.0, 0.0])
    a = np.array([1.0, -2.0 * rS * np.cos(thetaS), rS**2])

    rng = np.random.default_rng(0)
    N = 6000
    x = rng.standard_normal(N)

    y_ref = iir2_df1(x, b, a, dtype=np.float64)
    y_df1 = iir2_df1(x, b, a, dtype=np.float32).astype(np.float64)
    y_df2 = iir2_df2(x, b, a, dtype=np.float32).astype(np.float64)
    e1 = y_df1 - y_ref
    e2 = y_df2 - y_ref

    fig, ax = plt.subplots(2, 1, figsize=(8.0, 6.0), sharex=True)
    idx = np.arange(800)
    ax[0].plot(idx, y_ref[:800], "k-", lw=1.0, label="double ref")
    ax[0].plot(idx, y_df1[:800], "b--", lw=1.0, label="df1 single")
    ax[0].plot(idx, y_df2[:800], "r--", lw=1.0, label="df2 single")
    ax[0].grid(True)
    ax[0].set_ylabel("y[n]")
    ax[0].set_title("Output waveform (reference vs single-precision)")
    ax[0].legend(loc="best")
    ax[1].plot(idx, e1[:800], "b-", lw=1.0, label="df1 single - ref")
    ax[1].plot(idx, e2[:800], "r-", lw=1.0, label="df2 single - ref")
    ax[1].grid(True)
    ax[1].set_xlabel("n")
    ax[1].set_ylabel("error")
    ax[1].set_title("Error relative to double")
    ax[1].legend(loc="best")
    savefig(fig, out_dir / "ch09_structure_numeric_error.png")
    plt.close(fig)

    # 8) direct unit-circle evaluation vs FFT of truncated impulse response
    rF = 0.99
    thetaF = 0.20 * np.pi
    b = np.array([1.0, 0.0, 0.0])
    a = np.array([1.0, -2.0 * rF * np.cos(thetaF), rF**2])

    M = 256
    Nfft = 4096
    x = np.zeros(M)
    x[0] = 1.0
    h = iir2_df1(x, b, a)

    Hfft = np.fft.rfft(h, n=Nfft)
    f = np.arange(Hfft.size) * (fs / Nfft)
    f2, Hdir = freq_response(fs, b, a, n=Hfft.size)
    # Use the same frequency grid as rfft
    if np.max(np.abs(f2 - f)) > 1e-12:
        # Recompute direct H on the rfft grid
        omega = 2.0 * np.pi * f / fs
        mb = np.arange(b.size)
        ma = np.arange(a.size)
        B = np.exp(-1j * omega[:, None] * mb[None, :]) @ b.astype(np.complex128)
        A = np.exp(-1j * omega[:, None] * ma[None, :]) @ a.astype(np.complex128)
        Hdir = B / A

    fig, ax = plt.subplots(1, 1, figsize=(8.0, 3.8))
    ax.plot(f, 20 * np.log10(np.abs(Hdir) + 1e-20), lw=1.2, label="direct unit-circle eval")
    ax.plot(f, 20 * np.log10(np.abs(Hfft) + 1e-20), lw=1.2, label="FFT(truncated h)")
    ax.grid(True)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Magnitude [dB]")
    ax.set_xlim(0, fs / 2)
    ax.set_title(f"Direct H(e^{{jΩ}}) vs FFT of truncated h[n] (M={M}, Nfft={Nfft})")
    ax.legend(loc="best")
    savefig(fig, out_dir / "ch09_freqresp_direct_vs_fft.png")
    plt.close(fig)

    return 0


if __name__ == "__main__":
    # Ensure Matplotlib uses a writable config dir in sandboxes
    os.environ.setdefault("MPLCONFIGDIR", str(Path.cwd() / "build" / "mplconfig"))
    raise SystemExit(main())


