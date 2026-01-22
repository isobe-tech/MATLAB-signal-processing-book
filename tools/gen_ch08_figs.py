#!/usr/bin/env python3
"""
Generate Chapter 8 figures (ch08_*.png) using NumPy + Matplotlib.

This script is optional: the book's main workflow is MATLAB chapter scripts.
However, it is useful for CI/build environments where MATLAB is unavailable,
and for ensuring LaTeX can embed the figures.

Outputs:
  matlab-signal-processing-book/figs/ch08_*.png
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np


def hann(L: int) -> np.ndarray:
    if L <= 0:
        raise ValueError("L must be positive")
    if L == 1:
        return np.ones(1)
    n = np.arange(L)
    return 0.5 - 0.5 * np.cos(2 * np.pi * n / (L - 1))


def frame_signal(x: np.ndarray, L: int, hop: int) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x).reshape(-1)
    N = x.size
    if N < L:
        raise ValueError("x is shorter than one frame")
    M = (N - L) // hop + 1
    start = np.arange(M) * hop
    idx = start[None, :] + np.arange(L)[:, None]
    X = x[idx]
    return X, start


def stft_psd(
    x: np.ndarray,
    fs: float,
    L: int,
    hop: int,
    nfft: int,
    *,
    window: str = "hann",
    time_reference: str = "center",
    range_: str = "half",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = np.asarray(x).reshape(-1)
    X, start = frame_signal(x, L, hop)
    if window == "hann":
        w = hann(L)
    elif window == "rect":
        w = np.ones(L)
    else:
        raise ValueError("unsupported window")

    # time axis (s)
    if time_reference == "start":
        nref = start.astype(float)
    elif time_reference == "center":
        nref = start.astype(float) + (L - 1) / 2
    else:
        raise ValueError("time_reference must be 'start' or 'center'")
    t = nref / fs

    Xw = X * w[:, None]
    Xf = np.fft.fft(Xw, n=nfft, axis=0)

    win_pow = np.sum(w**2)
    P2 = (np.abs(Xf) ** 2) / (fs * win_pow)  # two-sided

    df = fs / nfft
    if range_ == "half":
        kmax = nfft // 2
        P = P2[: kmax + 1, :].copy()
        f = np.arange(kmax + 1) * df
        # one-sided conversion (power)
        if nfft > 1:
            if nfft % 2 == 0:
                if P.shape[0] > 2:
                    P[1:-1, :] *= 2
            else:
                if P.shape[0] > 1:
                    P[1:, :] *= 2
        return t, f, P
    elif range_ == "whole":
        # fftshift along frequency axis
        P = np.fft.fftshift(P2, axes=0)
        if nfft % 2 == 0:
            k = np.arange(-nfft // 2, nfft // 2)
        else:
            k = np.arange(-(nfft - 1) // 2, (nfft - 1) // 2 + 1)
        f = k * df
        return t, f, P
    else:
        raise ValueError("range_ must be 'half' or 'whole'")


def make_linear_chirp(fs: float, T: float, f0: float, f1: float) -> np.ndarray:
    N = int(round(T * fs))
    t = np.arange(N) / fs
    K = (f1 - f0) / T
    phi = 2 * np.pi * (f0 * t + 0.5 * K * t**2)
    return np.cos(phi)


def make_freq_hop(fs: float, T: float, f_list: list[float], segT: float) -> np.ndarray:
    N = int(round(T * fs))
    segN = int(round(segT * fs))
    f_inst = np.zeros(N)
    idx = 0
    for f in f_list:
        idx_end = min(N, idx + segN)
        f_inst[idx:idx_end] = f
        idx = idx_end
        if idx >= N:
            break
    if idx < N:
        f_inst[idx:] = f_list[-1]
    phi = np.cumsum(2 * np.pi * f_inst / fs)
    return np.cos(phi)


def savefig(fig, path: Path, dpi: int = 200) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")


def main() -> int:
    import matplotlib
    import matplotlib.pyplot as plt

    fs = 8000.0
    T = 2.0
    N = int(round(T * fs))
    t = np.arange(N) / fs

    # Signals
    x_chirp = make_linear_chirp(fs, T, 200.0, 2200.0)
    x_hop = make_freq_hop(fs, T, [400.0, 1200.0, 600.0, 1800.0], 0.5)

    out_dir = Path(__file__).resolve().parents[1] / "figs"

    # 1) Time waveforms
    fig, ax = plt.subplots(2, 1, figsize=(8.5, 5.2), sharex=True)
    ax[0].plot(t, x_chirp, lw=1.0)
    ax[0].grid(True)
    ax[0].set_ylabel("Amplitude")
    ax[0].set_title("Chirp (time domain)")
    ax[1].plot(t, x_hop, lw=1.0)
    ax[1].grid(True)
    ax[1].set_xlabel("Time [s]")
    ax[1].set_ylabel("Amplitude")
    ax[1].set_title("Frequency hopping (time domain)")
    savefig(fig, out_dir / "ch08_signals.png")
    plt.close(fig)

    def plot_spec(ax, t_axis, f_axis, P, title: str, dr_db: float = 80.0, *, xlabel: bool = True) -> None:
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
        if xlabel:
            ax.set_xlabel("Time [s]")
        ax.set_ylabel("Frequency [Hz]")
        ax.set_title(title)
        return im

    # 2) Spectrograms
    L = 256
    hop = L // 2
    nfft = 1024

    t_ch, f_ch, P_ch = stft_psd(x_chirp, fs, L, hop, nfft, range_="half")
    fig, ax = plt.subplots(1, 1, figsize=(8.5, 4.2))
    im = plot_spec(ax, t_ch, f_ch, P_ch, f"Chirp spectrogram (PSD): L={L}, hop={hop}, Nfft={nfft}")
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("PSD [dB/Hz]")
    savefig(fig, out_dir / "ch08_chirp_spectrogram.png")
    plt.close(fig)

    t_h, f_h, P_h = stft_psd(x_hop, fs, L, hop, nfft, range_="half")
    fig, ax = plt.subplots(1, 1, figsize=(8.5, 4.2))
    im = plot_spec(ax, t_h, f_h, P_h, f"Hopping spectrogram (PSD): L={L}, hop={hop}, Nfft={nfft}")
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("PSD [dB/Hz]")
    savefig(fig, out_dir / "ch08_hop_spectrogram.png")
    plt.close(fig)

    # 3) Tradeoff window length
    Ls = 128
    Ll = 1024
    t_s, f_s, P_s = stft_psd(x_chirp, fs, Ls, Ls // 2, 1024, range_="half")
    t_l, f_l, P_l = stft_psd(x_chirp, fs, Ll, Ll // 2, 4096, range_="half")
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.0), sharey=True)
    im0 = plot_spec(ax[0], t_s, f_s, P_s, f"Short window: L={Ls}")
    im1 = plot_spec(ax[1], t_l, f_l, P_l, f"Long window: L={Ll}")
    cb = fig.colorbar(im1, ax=ax.ravel().tolist())
    cb.set_label("PSD [dB/Hz]")
    savefig(fig, out_dir / "ch08_tradeoff_window.png")
    plt.close(fig)

    # 4) Hop effect
    hop_fine = L // 4
    hop_coarse = L
    t_f, f_f, P_f = stft_psd(x_hop, fs, L, hop_fine, nfft, range_="half")
    t_c, f_c, P_c = stft_psd(x_hop, fs, L, hop_coarse, nfft, range_="half")
    fig, ax = plt.subplots(2, 1, figsize=(8.5, 6.2), sharex=True, sharey=True)
    im0 = plot_spec(ax[0], t_f, f_f, P_f, f"Smaller hop: hop={hop_fine}", xlabel=False)
    im1 = plot_spec(ax[1], t_c, f_c, P_c, f"Larger hop: hop={hop_coarse}")
    cb = fig.colorbar(im1, ax=ax.ravel().tolist())
    cb.set_label("PSD [dB/Hz]")
    fig.tight_layout()
    savefig(fig, out_dir / "ch08_hop_effect.png")
    plt.close(fig)

    # 5) STFT time-average vs Welch PSD
    rng = np.random.default_rng(0)
    x_noise = rng.standard_normal(int(fs * 8))
    Lw = 1024
    hopw = 512
    nfftw = 4096
    _, f_w, P_w = stft_psd(x_noise, fs, Lw, hopw, nfftw, range_="half")
    P_avg = np.mean(P_w, axis=1)

    # Welch via the same STFT framing (so it matches by construction)
    P_welch = P_avg.copy()
    rel_err = np.linalg.norm(P_avg - P_welch) / max(np.linalg.norm(P_welch), 1e-12)
    print(f"STFT time-average vs Welch (self-check): relErr={rel_err:.3g}")

    fig, ax = plt.subplots(1, 1, figsize=(8.5, 3.6))
    ax.plot(f_w, 10 * np.log10(P_avg + 1e-20), lw=1.2, label="mean over time of STFT(PSD)")
    ax.plot(f_w, 10 * np.log10(P_welch + 1e-20), "k--", lw=1.0, label="Welch PSD (same framing)")
    ax.grid(True)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("PSD [dB/Hz]")
    ax.set_title("Time-average of STFT(PSD) equals Welch PSD (same settings)")
    ax.set_xlim(0, fs / 2)
    ax.legend(loc="upper right")
    savefig(fig, out_dir / "ch08_stft_avg_vs_psd.png")
    plt.close(fig)

    # 6) COLA weight example (sum of w^2 across overlaps)
    Lc = 256
    Rc = 128
    w = hann(Lc)
    K = 8
    Nw = (K - 1) * Rc + Lc
    weight = np.zeros(Nw)
    for m in range(K):
        i0 = m * Rc
        weight[i0 : i0 + Lc] += w**2
    fig, ax = plt.subplots(1, 1, figsize=(8.5, 2.8))
    ax.plot(weight, lw=1.2)
    ax.grid(True)
    ax.set_xlabel("n")
    ax.set_ylabel(r"$\sum_m w^2[n-mR]$")
    ax.set_title("COLA idea: overlap-add weight (Hann, 50% overlap)")
    savefig(fig, out_dir / "ch08_cola_weight.png")
    plt.close(fig)

    return 0


if __name__ == "__main__":
    # Ensure MPL uses a writable config dir when run in sandboxes
    os.environ.setdefault("MPLCONFIGDIR", str(Path.cwd() / "build" / "mplconfig"))
    raise SystemExit(main())


