%% 第6章 DFTの直接実装とFFTによる検証
% このスクリプトは「Run All」で図まで再現できる形にする。

clear; close all; clc;

% ルート（bookRoot）を推定して src/ を path に追加
thisFile = mfilename("fullpath");
bookRoot = fileparts(fileparts(thisFile));
addpath(fullfile(bookRoot, "src"));
sp_book_style(); % White background + readable fonts for LaTeX

% 出力先（図）
figDir = fullfile(bookRoot, "figs");
if ~exist(figDir, "dir")
    mkdir(figDir);
end

%% Parameters (keep N modest: direct DFT is O(N^2))
fs = 8000;   % [Hz]
N  = 256;    % samples (DFT length)
n  = (0:N-1).';

%% 0) Quick checks: direct DFT / matrix DFT must match fft
rng(0);
xChk = randn(N, 1) + 1j*randn(N, 1);

X0 = sp_dft(xChk);
X1 = fft(xChk);
relErr = norm(X0 - X1) / max(norm(X1), 1);
fprintf("sp_dft vs fft: relErr=%.3g\n", relErr);
assert(relErr < 1e-11, "sp_dft mismatch.");

[W, Winv] = sp_dft_matrix(N);
X2 = W * xChk;
relErr = norm(X2 - X1) / max(norm(X1), 1);
fprintf("sp_dft_matrix vs fft: relErr=%.3g\n", relErr);
assert(relErr < 1e-11, "sp_dft_matrix mismatch.");

xRec = Winv * X1;
relErr = norm(xRec - xChk) / max(norm(xChk), 1);
fprintf("inverse (Winv) recovery: relErr=%.3g\n", relErr);
assert(relErr < 1e-11, "inverse DFT mismatch.");

%% 1) Ideal case: complex exponential -> one bin (delta in DFT)
k0 = 14;                         % choose an integer bin
f0 = k0 * fs / N;                % [Hz] exactly on a bin
xExp = exp(1j * 2*pi*(k0/N) * n); % e^{j2pi k0 n/N}

X = fft(xExp);
mag = abs(fftshift(X)) / N;
[fWhole, ~] = sp_fft_axis(fs, N, "Range", "whole");

figure;
stem(fWhole, mag, "filled", "MarkerSize", 2);
grid on;
xlabel("Frequency [Hz]");
ylabel("Magnitude (|X|/N)");
title(sprintf("Complex exponential: x[n]=exp(j2\\pi k0 n/N), k0=%d (f0=%.1f Hz)", k0, f0));
xlim([-fs/2, fs/2 - fs/N]);
exportgraphics(gcf, fullfile(figDir, "ch06_complex_exponential_delta.png"), "Resolution", 200, "BackgroundColor", "white");

%% 2) Real cosine: two-sided spectrum (+f and -f)
xCos = cos(2*pi*f0*n/fs);
X = fft(xCos);

% two-sided (fftshift)
mag2 = abs(fftshift(X)) / N;

figure;
stem(fWhole, mag2, "filled", "MarkerSize", 2);
grid on;
xlabel("Frequency [Hz]");
ylabel("Magnitude (two-sided, |X|/N)");
title(sprintf("Real cosine: two symmetric peaks at \\pm f0 (f0=%.1f Hz)", f0));
xlim([-fs/2, fs/2 - fs/N]);
exportgraphics(gcf, fullfile(figDir, "ch06_real_cosine_two_sided.png"), "Resolution", 200, "BackgroundColor", "white");

%% 3) Frequency axis + scaling: half (one-sided) vs whole (two-sided)
[fHalf, ~] = sp_fft_axis(fs, N, "Range", "half");
kMax = numel(fHalf) - 1;

magHalf = abs(X) / N;
magHalf = magHalf(1:kMax+1);
if mod(N, 2) == 0
    if numel(magHalf) > 2
        magHalf(2:end-1) = 2 * magHalf(2:end-1);
    end
else
    if numel(magHalf) > 1
        magHalf(2:end) = 2 * magHalf(2:end);
    end
end

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
stem(fHalf, magHalf, "filled", "MarkerSize", 2);
grid on;
xlabel("Frequency [Hz]");
ylabel("Magnitude (one-sided)");
title("One-sided spectrum (0..fs/2): non-DC bins are doubled");

nexttile;
stem(fWhole, mag2, "filled", "MarkerSize", 2);
grid on;
xlabel("Frequency [Hz]");
ylabel("Magnitude (two-sided)");
title("Two-sided spectrum (-fs/2..fs/2): peaks split into +f and -f");
xlim([-fs/2, fs/2 - fs/N]);

exportgraphics(gcf, fullfile(figDir, "ch06_freq_axis_half_vs_whole.png"), "Resolution", 200, "BackgroundColor", "white");

%% 4) Leakage: bin-aligned vs bin-mismatched (finite-length observation)
fMis = 440; % [Hz] intentionally off-bin for N=256, fs=8000 (Δf=31.25 Hz)
xAligned = cos(2*pi*f0*n/fs);
xMis = cos(2*pi*fMis*n/fs);

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
sp_plot_mag_spec(fs, xAligned, "NewFigure", false, "Db", true, "OneSided", true, ...
    "Title", sprintf("Bin-aligned tone: f0=%.1f Hz (k0=%d)", f0, k0));
ylim([-90, 5]);

nexttile;
sp_plot_mag_spec(fs, xMis, "NewFigure", false, "Db", true, "OneSided", true, ...
    "Title", sprintf("Bin-mismatched tone: f0=%.1f Hz (leakage)", fMis));
ylim([-90, 5]);

exportgraphics(gcf, fullfile(figDir, "ch06_leakage_bin_mismatch.png"), "Resolution", 200, "BackgroundColor", "white");

%% 5) Zero padding: smoother-looking spectrum, but the "width" is not improved
Nfft1 = N;
Nfft2 = 8*N;

[f1, mag1] = one_sided_amp_spectrum(fs, xMis, Nfft1);
[f2, mag2p] = one_sided_amp_spectrum(fs, xMis, Nfft2);

figure;
plot(f1, 20*log10(mag1 + 1e-15), "o-", "LineWidth", 1.0, "MarkerSize", 3);
hold on;
plot(f2, 20*log10(mag2p + 1e-15), "-", "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]");
ylabel("Amplitude [dB]");
title(sprintf("Zero padding: N=%d vs Nfft=%d (same data, denser frequency grid)", N, Nfft2));
legend(["N (coarse bins)", "zero padded (dense bins)"], "Location", "northeast");
ylim([-90, 5]);
exportgraphics(gcf, fullfile(figDir, "ch06_zero_padding_interp.png"), "Resolution", 200, "BackgroundColor", "white");

%% 6) DFT matrix visualization: real part looks like a bank of cosines
Nshow = 64;
[Wshow, ~] = sp_dft_matrix(Nshow);

figure;
imagesc(real(Wshow));
axis image;
colormap(gray);
cb = colorbar;
try
    cb.Color = [0 0 0];
catch
end
xlabel("n (column)");
ylabel("k (row)");
title("Real part of the DFT matrix W (cosine patterns)");
exportgraphics(gcf, fullfile(figDir, "ch06_dft_matrix_real.png"), "Resolution", 200, "BackgroundColor", "white");

%% 7) Orthogonality: W^H W becomes (almost) diagonal
[Wu, ~] = sp_dft_matrix(Nshow, "Unitary", true);
G = abs(Wu' * Wu);

figure;
imagesc(G);
axis image;
colormap(gray);
cb = colorbar;
try
    cb.Color = [0 0 0];
catch
end
caxis([0, 1]);
xlabel("k");
ylabel("k");
title("|W^H W| for unitary DFT matrix (should be I)");
exportgraphics(gcf, fullfile(figDir, "ch06_dft_orthogonality.png"), "Resolution", 200, "BackgroundColor", "white");

%% 8) Timing: O(N^2) vs O(N log N)
Nlist = [16, 32, 64, 128, 256];
tDft = zeros(size(Nlist));
tFft = zeros(size(Nlist));

for i = 1:numel(Nlist)
    Ni = Nlist(i);
    rng(0);
    x = randn(Ni, 1) + 1j*randn(Ni, 1);

    if exist("timeit", "file")
        tDft(i) = timeit(@() sp_dft(x));
        tFft(i) = timeit(@() fft(x));
    else
        rep = 5;
        tic;
        for r = 1:rep %#ok<NASGU>
            sp_dft(x);
        end
        tDft(i) = toc / rep;

        tic;
        for r = 1:rep %#ok<NASGU>
            fft(x);
        end
        tFft(i) = toc / rep;
    end

    fprintf("timing: N=%4d, sp_dft=%.3g s, fft=%.3g s\n", Ni, tDft(i), tFft(i));
end

figure;
loglog(Nlist, tDft, "o-", "LineWidth", 1.2, "MarkerSize", 4);
hold on;
loglog(Nlist, tFft, "o-", "LineWidth", 1.2, "MarkerSize", 4);
grid on;
xlabel("N");
ylabel("time [s]");
title("Timing comparison: direct DFT vs FFT");
legend(["sp\_dft (O(N^2))", "fft (fast)"], "Location", "northwest");
exportgraphics(gcf, fullfile(figDir, "ch06_dft_vs_fft_timing.png"), "Resolution", 200, "BackgroundColor", "white");

%% Local helper: one-sided amplitude spectrum (for real signals)
function [f, amp] = one_sided_amp_spectrum(fs, x, Nfft)
x = x(:);
N = numel(x);
X = fft(x, Nfft);
amp = abs(X) / N;

kMax = floor(Nfft/2);
amp = amp(1:kMax+1);
f = (0:kMax).' * (fs / Nfft);

if mod(Nfft, 2) == 0
    if numel(amp) > 2
        amp(2:end-1) = 2 * amp(2:end-1);
    end
else
    if numel(amp) > 1
        amp(2:end) = 2 * amp(2:end);
    end
end
end
