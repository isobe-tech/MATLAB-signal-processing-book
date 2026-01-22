%% 第7章 窓関数とスペクトル推定：漏れ・分解能・分散
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

%% Parameters (leakage + window demo)
fs = 8000;   % [Hz]
N  = 512;    % samples for short demo
n  = (0:N-1).';
t  = n / fs;

fTone = 440; % [Hz] intentionally off-bin (Δf = fs/N = 15.625 Hz)
xTone = cos(2*pi*fTone*t);

wRect = sp_window_rect(N);
wHann = sp_window_hann(N);
wHamm = sp_window_hamming(N);

cgRect = mean(wRect);
cgHann = mean(wHann);
cgHamm = mean(wHamm);

NfftLeak = 8192;

%% 0) Window shapes in time domain (what does "tapering" mean?)
Nw = 64;
nW = (0:Nw-1).';
wR = sp_window_rect(Nw);
wH = sp_window_hann(Nw);
wM = sp_window_hamming(Nw);

figure;
plot(nW, wR, "LineWidth", 1.2); hold on;
plot(nW, wH, "LineWidth", 1.2);
plot(nW, wM, "LineWidth", 1.2);
grid on;
xlabel("n");
ylabel("w[n]");
title("Window shapes in time domain: rect vs Hann vs Hamming");
legend(["rect", "hann", "hamming"], "Location", "best");

exportgraphics(gcf, fullfile(figDir, "ch07_windows_time.png"), "Resolution", 200, "BackgroundColor", "white");

%% 1) Leakage shape: rectangular vs Hann vs Hamming (raw vs coherent-gain normalized)
[fLeak, aRectRaw] = one_sided_amp_spectrum(fs, xTone .* wRect, NfftLeak);
[~,     aHannRaw] = one_sided_amp_spectrum(fs, xTone .* wHann, NfftLeak);
[~,     aHammRaw] = one_sided_amp_spectrum(fs, xTone .* wHamm, NfftLeak);

aRect = aRectRaw / cgRect;
aHann = aHannRaw / cgHann;
aHamm = aHammRaw / cgHamm;

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
plot(fLeak, 20*log10(aRectRaw + 1e-15), "LineWidth", 1.2); hold on;
plot(fLeak, 20*log10(aHannRaw + 1e-15), "LineWidth", 1.2);
plot(fLeak, 20*log10(aHammRaw + 1e-15), "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]");
ylabel("Amplitude [dB]");
title("Leakage (raw): window changes the peak height as well as sidelobes");
legend(["rect (raw)", "hann (raw)", "hamming (raw)"], "Location", "northeast");
xlim([0, 2000]);
ylim([-120, 5]);

nexttile;
plot(fLeak, 20*log10(aRect + 1e-15), "LineWidth", 1.2); hold on;
plot(fLeak, 20*log10(aHann + 1e-15), "LineWidth", 1.2);
plot(fLeak, 20*log10(aHamm + 1e-15), "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]");
ylabel("Amplitude [dB]");
title("Leakage (coherent-gain normalized): compare sidelobes fairly");
legend(["rect", "hann", "hamming"], "Location", "northeast");
xlim([0, 2000]);
ylim([-120, 5]);

exportgraphics(gcf, fullfile(figDir, "ch07_leakage_window_norm.png"), "Resolution", 200, "BackgroundColor", "white");

%% 2) Window spectra: main lobe vs sidelobes (normalize DC to 0 dB)
[fWin, Wr] = one_sided_window_spectrum(fs, wRect, NfftLeak);
[~,    Wh] = one_sided_window_spectrum(fs, wHann, NfftLeak);
[~,    Wm] = one_sided_window_spectrum(fs, wHamm, NfftLeak);

figure;
plot(fWin, 20*log10(Wr + 1e-15), "LineWidth", 1.2); hold on;
plot(fWin, 20*log10(Wh + 1e-15), "LineWidth", 1.2);
plot(fWin, 20*log10(Wm + 1e-15), "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]");
ylabel("Normalized magnitude [dB]");
title("Window spectra (DC normalized): main lobe width vs sidelobe level");
legend(["rect", "hann", "hamming"], "Location", "northeast");
xlim([0, 1000]);
ylim([-120, 5]);

exportgraphics(gcf, fullfile(figDir, "ch07_window_spectra.png"), "Resolution", 200, "BackgroundColor", "white");

%% 2b) Scalloping loss: peak drop vs bin offset (complex tone)
Nsc = 256;
nSc = (0:Nsc-1).';
k0 = 25; % reference bin
delta = linspace(0, 0.5, 201); % [bin] offset (symmetry: +/-delta is the same)

wRectSc = sp_window_rect(Nsc);
wHannSc = sp_window_hann(Nsc);
wHammSc = sp_window_hamming(Nsc);

cgRectSc = mean(wRectSc);
cgHannSc = mean(wHannSc);
cgHammSc = mean(wHammSc);

slRect = zeros(size(delta));
slHann = zeros(size(delta));
slHamm = zeros(size(delta));

for i = 1:numel(delta)
    d = delta(i);
    x = exp(1j * 2*pi*((k0 + d)/Nsc) * nSc); % unit-amplitude complex tone

    Xr = fft(x .* wRectSc);
    Xh = fft(x .* wHannSc);
    Xm = fft(x .* wHammSc);

    % peak amplitude estimate (DFT bin peak), coherent-gain normalized
    ar = max(abs(Xr)) / Nsc / cgRectSc;
    ah = max(abs(Xh)) / Nsc / cgHannSc;
    am = max(abs(Xm)) / Nsc / cgHammSc;

    slRect(i) = 20*log10(ar + 1e-15);
    slHann(i) = 20*log10(ah + 1e-15);
    slHamm(i) = 20*log10(am + 1e-15);
end

figure;
plot(delta, slRect, "LineWidth", 1.2); hold on;
plot(delta, slHann, "LineWidth", 1.2);
plot(delta, slHamm, "LineWidth", 1.2);
grid on;
xlabel("Bin offset |\delta| [bin]");
ylabel("Peak amplitude [dB]");
title("Scalloping loss: peak amplitude vs bin offset (coherent-gain normalized)");
legend(["rect", "hann", "hamming"], "Location", "southwest");
ylim([-6, 1]);

exportgraphics(gcf, fullfile(figDir, "ch07_scalloping_loss.png"), "Resolution", 200, "BackgroundColor", "white");

%% 3) PSD: single periodogram vs averaged PSD (white noise)
rng(0);
Nlong = fs * 8; % 8 s
xNoise = randn(Nlong, 1);

frameLength = 1024;
hop = frameLength/2;
nfftPsd = 4096;

% single-frame PSD (periodogram-like)
[fP, P1] = sp_psd_avg(xNoise(1:frameLength), fs, ...
    "FrameLength", frameLength, "Hop", frameLength, "Window", "hann", "Nfft", nfftPsd, "Range", "half");

% averaged PSD (Welch-like)
[~, Pavg, infoAvg] = sp_psd_avg(xNoise, fs, ...
    "FrameLength", frameLength, "Hop", hop, "Window", "hann", "Nfft", nfftPsd, "Range", "half");

df = infoAvg.df;
pTime = mean(xNoise.^2);
pFreq = sum(Pavg) * df;
relErr = abs(pFreq - pTime) / max(pTime, 1e-12);
fprintf("PSD power check (white noise): time=%.4g, freq=%.4g, relErr=%.3g\n", pTime, pFreq, relErr);
assert(relErr < 0.05, "PSD scaling/power check failed (unexpected mismatch).");

Ptheory = (2*pTime/fs) * ones(size(fP)); % one-sided (approx; DC/Nyquist detail ignored)

figure;
plot(fP, 10*log10(P1 + 1e-20), "LineWidth", 1.0); hold on;
plot(fP, 10*log10(Pavg + 1e-20), "LineWidth", 1.2);
plot(fP, 10*log10(Ptheory + 1e-20), "k--", "LineWidth", 1.0);
grid on;
xlabel("Frequency [Hz]");
ylabel("PSD [dB/Hz]");
title(sprintf("PSD variance reduction: single vs averaged (nFrames=%d)", infoAvg.nFrames));
legend(["single (periodogram)", "averaged (Welch-like)", "theory (flat)"], "Location", "northeast");
xlim([0, fs/2]);

exportgraphics(gcf, fullfile(figDir, "ch07_psd_single_vs_avg.png"), "Resolution", 200, "BackgroundColor", "white");

%% 3b) PSD averaging levels: how many averages are "enough"?
% Build per-frame periodograms once, then average the first M frames.
wWelch = sp_window_hann(frameLength);
winPow = sum(wWelch.^2);
[Xn, ~] = sp_frame_signal(xNoise, frameLength, hop, "PadEnd", false);
Xf = fft(Xn .* wWelch, nfftPsd, 1); % nfft x nFrames
P2_frames = abs(Xf).^2 / (fs * winPow); % two-sided PSD per frame

nFrames = size(P2_frames, 2);
Mlist = [1, 4, 16, 64];
Mlist = Mlist(Mlist <= nFrames);

figure;
colors = lines(numel(Mlist));
for i = 1:numel(Mlist)
    M = Mlist(i);
    P2m = mean(P2_frames(:, 1:M), 2);
    [fM, Pm] = one_sided_psd_from_two_sided(fs, P2m);
    plot(fM, 10*log10(Pm + 1e-20), "LineWidth", 1.2, "Color", colors(i, :)); hold on;
end
grid on;
xlabel("Frequency [Hz]");
ylabel("PSD [dB/Hz]");
title("PSD averaging levels: more averages -> smoother (lower variance)");
legend("M=1", "M=4", "M=16", "M=64", "Location", "northeast");
xlim([0, fs/2]);

exportgraphics(gcf, fullfile(figDir, "ch07_psd_avg_levels.png"), "Resolution", 200, "BackgroundColor", "white");

%% 4) Tradeoff: resolution vs variance (two close tones + noise)
rng(1);
Nsig = fs * 2; % 2 s
n2 = (0:Nsig-1).';
t2 = n2 / fs;
f1 = 440;   % [Hz]
f2 = 470;   % [Hz] close-by tone

xTwo = cos(2*pi*f1*t2) + 0.8*cos(2*pi*f2*t2) + 0.3*randn(Nsig, 1);

% Short frames: more averaging -> smoother, but poorer resolution
[fS, PS] = sp_psd_avg(xTwo, fs, ...
    "FrameLength", 256, "Hop", 128, "Window", "hann", "Nfft", 2048, "Range", "half");

% Long frames: better resolution, less averaging -> noisier
[fL, PL, infoL] = sp_psd_avg(xTwo, fs, ...
    "FrameLength", 2048, "Hop", 1024, "Window", "hann", "Nfft", 2048, "Range", "half");

figure;
plot(fS, 10*log10(PS + 1e-20), "LineWidth", 1.2); hold on;
plot(fL, 10*log10(PL + 1e-20), "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]");
ylabel("PSD [dB/Hz]");
title(sprintf("Tradeoff: short frames vs long frames (long nFrames=%d)", infoL.nFrames));
legend(["short L=256 (more avg)", "long L=2048 (better resolution)"], "Location", "northeast");
xlim([300, 700]);

exportgraphics(gcf, fullfile(figDir, "ch07_psd_tradeoff_two_tones.png"), "Resolution", 200, "BackgroundColor", "white");

%% 5) Framing as a matrix: each column is a frame (bridge to STFT)
frameVisLen = 128;
hopVis = 32;
xVis = xTwo(1:2048);
Xframes = sp_frame_signal(xVis, frameVisLen, hopVis, "PadEnd", false);

figure;
imagesc(Xframes);
axis tight;
colormap(gray);
cb = colorbar;
try
    cb.Color = [0 0 0];
catch
end
xlabel("frame index (column)");
ylabel("sample index within a frame (row)");
title(sprintf("Frame matrix: size %d x %d (L x nFrames)", size(Xframes, 1), size(Xframes, 2)));

exportgraphics(gcf, fullfile(figDir, "ch07_frame_matrix.png"), "Resolution", 200, "BackgroundColor", "white");

%% Local helpers
function [f, amp] = one_sided_amp_spectrum(fs, x, Nfft)
x = x(:);
N = numel(x);
X = fft(x, Nfft);

amp = abs(X) / N;
kMax = floor(Nfft/2);
amp = amp(1:kMax+1);
f = (0:kMax).' * (fs / Nfft);

if Nfft > 1
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
end

function [f, mag] = one_sided_window_spectrum(fs, w, Nfft)
w = w(:);
W = fft(w, Nfft);
mag = abs(W) / max(sum(w), 1e-12); % normalize DC to 1
kMax = floor(Nfft/2);
mag = mag(1:kMax+1);
f = (0:kMax).' * (fs / Nfft);
end

function [f, P1] = one_sided_psd_from_two_sided(fs, P2)
P2 = P2(:);
nfft = numel(P2);
kMax = floor(nfft/2);
P1 = P2(1:kMax+1);
f = (0:kMax).' * (fs / nfft);

if nfft > 1
    if mod(nfft, 2) == 0
        if numel(P1) > 2
            P1(2:end-1) = 2 * P1(2:end-1);
        end
    else
        if numel(P1) > 1
            P1(2:end) = 2 * P1(2:end);
        end
    end
end
end
