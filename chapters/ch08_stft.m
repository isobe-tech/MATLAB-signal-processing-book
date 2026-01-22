%% 第8章 STFT（スペクトログラム）の実装：非定常信号の解析
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

%% Parameters
fs = 8000;          % [Hz]
T  = 2.0;           % [s]
N  = round(T * fs); % samples
[t, n] = sp_time_axis(fs, N); %#ok<ASGLU>

%% 1) Example non-stationary signals: chirp and frequency hopping
f0 = 200;   % [Hz]
f1 = 2200;  % [Hz]
xChirp = make_linear_chirp(fs, T, f0, f1);

fSeq = [400, 1200, 600, 1800];  % [Hz]
segT = 0.5;                     % [s] each segment
xHop = make_freq_hop(fs, T, fSeq, segT);

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
% Plot a short-time zoom; the full 2 s waveform is too dense to read.
tWinChirp = [0, 0.05]; % [s]
idxCh = (t >= tWinChirp(1)) & (t <= tWinChirp(2));
plot(t(idxCh), xChirp(idxCh), "LineWidth", 1.0);
grid on;
xlabel("Time [s]");
ylabel("Amplitude");
ylim([-1.2, 1.2]);
xlim(tWinChirp);
title(sprintf("Chirp (waveform zoom): t=%.3f..%.3f s, f(t)=%d..%d Hz (over %.1f s)", tWinChirp(1), tWinChirp(2), f0, f1, T));

nexttile;
% Zoom around the first hop boundary (shows frequency change clearly).
tHop0 = segT;               % first boundary [s]
tWinHop = [tHop0 - 0.01, tHop0 + 0.01]; % 20 ms window
idxHop = (t >= tWinHop(1)) & (t <= tWinHop(2));
plot(t(idxHop), xHop(idxHop), "LineWidth", 1.0);
grid on;
xlabel("Time [s]");
ylabel("Amplitude");
ylim([-1.2, 1.2]);
xlim(tWinHop);
% Mark the hop timing and annotate frequencies.
xline(tHop0, "k--", "hop", "LabelOrientation", "horizontal", "LabelVerticalAlignment", "middle");
text(tHop0 - 0.009, 1.05, sprintf("%d Hz", fSeq(1)), "HorizontalAlignment", "left");
text(tHop0 + 0.001, 1.05, sprintf("%d Hz", fSeq(2)), "HorizontalAlignment", "left");
title("Frequency hopping (waveform zoom around a hop)");

exportgraphics(gcf, fullfile(figDir, "ch08_signals.png"), "Resolution", 200, "BackgroundColor", "white");

%% 2) Spectrograms (PSD) for chirp and hopping signals
L = 256;
hop = L/2;
Nfft = 1024;

[tCh, fCh, Pch] = sp_stft(xChirp, fs, ...
    "FrameLength", L, "Hop", hop, "Window", "hann", "Nfft", Nfft, ...
    "Range", "half", "PadEnd", false, "TimeReference", "center", "Output", "psd");

figure;
sp_spectrogram_plot(tCh, fCh, Pch, ...
    "Title", sprintf("Chirp spectrogram (PSD): L=%d, hop=%d, Nfft=%d", L, hop, Nfft), ...
    "DynamicRangeDb", 80, "NewFigure", false);
exportgraphics(gcf, fullfile(figDir, "ch08_chirp_spectrogram.png"), "Resolution", 200, "BackgroundColor", "white");

[tH, fH, Ph] = sp_stft(xHop, fs, ...
    "FrameLength", L, "Hop", hop, "Window", "hann", "Nfft", Nfft, ...
    "Range", "half", "PadEnd", false, "TimeReference", "center", "Output", "psd");

figure;
sp_spectrogram_plot(tH, fH, Ph, ...
    "Title", sprintf("Hopping spectrogram (PSD): L=%d, hop=%d, Nfft=%d", L, hop, Nfft), ...
    "DynamicRangeDb", 80, "NewFigure", false);
exportgraphics(gcf, fullfile(figDir, "ch08_hop_spectrogram.png"), "Resolution", 200, "BackgroundColor", "white");

%% 3) Tradeoff: window length changes time/frequency resolution
Lshort = 128;
Llong  = 1024;

[tS, fS, PS] = sp_stft(xChirp, fs, ...
    "FrameLength", Lshort, "Hop", Lshort/2, "Window", "hann", "Nfft", 1024, ...
    "Range", "half", "PadEnd", false, "TimeReference", "center", "Output", "psd");

[tL, fL, PL] = sp_stft(xChirp, fs, ...
    "FrameLength", Llong, "Hop", Llong/2, "Window", "hann", "Nfft", 4096, ...
    "Range", "half", "PadEnd", false, "TimeReference", "center", "Output", "psd");

figure;
tiledlayout(1, 2, "Padding", "compact", "TileSpacing", "compact");

nexttile;
sp_spectrogram_plot(tS, fS, PS, ...
    "Title", sprintf("Short window: L=%d (better time resolution)", Lshort), ...
    "DynamicRangeDb", 80, "NewFigure", false);

nexttile;
sp_spectrogram_plot(tL, fL, PL, ...
    "Title", sprintf("Long window: L=%d (better frequency resolution)", Llong), ...
    "DynamicRangeDb", 80, "NewFigure", false);

exportgraphics(gcf, fullfile(figDir, "ch08_tradeoff_window.png"), "Resolution", 200, "BackgroundColor", "white");

%% 4) Hop effect: same window, different time sampling of the spectrogram
hopFine = L/4;
hopCoarse = L; % no overlap

[tF, fF, PF] = sp_stft(xHop, fs, ...
    "FrameLength", L, "Hop", hopFine, "Window", "hann", "Nfft", Nfft, ...
    "Range", "half", "PadEnd", false, "TimeReference", "center", "Output", "psd");

[tC, fC, PC] = sp_stft(xHop, fs, ...
    "FrameLength", L, "Hop", hopCoarse, "Window", "hann", "Nfft", Nfft, ...
    "Range", "half", "PadEnd", false, "TimeReference", "center", "Output", "psd");

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

axTop = nexttile;
sp_spectrogram_plot(tF, fF, PF, ...
    "Title", sprintf("Smaller hop: hop=%d (denser time grid)", hopFine), ...
    "DynamicRangeDb", 80, "NewFigure", false);
% Avoid label/title collisions between stacked axes (top x-label is redundant).
xlabel(axTop, "");
xticklabels(axTop, []);

nexttile;
sp_spectrogram_plot(tC, fC, PC, ...
    "Title", sprintf("Larger hop: hop=%d (sparser time grid)", hopCoarse), ...
    "DynamicRangeDb", 80, "NewFigure", false);

exportgraphics(gcf, fullfile(figDir, "ch08_hop_effect.png"), "Resolution", 200, "BackgroundColor", "white");

%% 5) STFT time-average vs Welch PSD (connection to Chapter 7)
rng(0);
Nnoise = fs * 8;
xNoise = randn(Nnoise, 1);

Lw = 1024;
hopw = 512;
Nfftw = 4096;

[tW, fW, Pw] = sp_stft(xNoise, fs, ...
    "FrameLength", Lw, "Hop", hopw, "Window", "hann", "Nfft", Nfftw, ...
    "Range", "half", "PadEnd", false, "TimeReference", "center", "Output", "psd");

PavgStft = mean(Pw, 2);
[fWelch, Pwelch] = sp_psd_avg(xNoise, fs, ...
    "FrameLength", Lw, "Hop", hopw, "Window", "hann", "Nfft", Nfftw, ...
    "Range", "half", "PadEnd", false);

relErr = norm(PavgStft - Pwelch) / max(norm(Pwelch), 1e-12);
fprintf("STFT time-average vs Welch PSD: relErr=%.3g\n", relErr);
assert(relErr < 1e-12, "STFT average PSD mismatch (should match sp_psd_avg).");

figure;
plot(fW, 10*log10(PavgStft + 1e-20), "LineWidth", 1.2); hold on;
plot(fWelch, 10*log10(Pwelch + 1e-20), "k--", "LineWidth", 1.0);
grid on;
xlabel("Frequency [Hz]");
ylabel("PSD [dB/Hz]");
title("Time-average of STFT(PSD) equals Welch PSD (same framing/windowing)");
legend(["mean over time of STFT(PSD)", "sp\_psd\_avg"], "Location", "northeast");
xlim([0, fs/2]);

exportgraphics(gcf, fullfile(figDir, "ch08_stft_avg_vs_psd.png"), "Resolution", 200, "BackgroundColor", "white");

%% 6) COLA weight example (sum of w^2 across overlaps) [preview for inverse STFT]
Lc = 256;
Rc = 128;
w = sp_window_hann(Lc);

K = 8;
Nw = (K-1)*Rc + Lc;
weight = zeros(Nw, 1);
for m = 0:K-1
    idx = (1:Lc) + m*Rc;
    weight(idx) = weight(idx) + w.^2;
end

figure;
plot(weight, "LineWidth", 1.2);
grid on;
xlabel("n");
ylabel("sum(w^2) across overlaps");
title("COLA idea: overlap-add weight (Hann, 50% overlap)");
exportgraphics(gcf, fullfile(figDir, "ch08_cola_weight.png"), "Resolution", 200, "BackgroundColor", "white");

%% Local helpers (no toolboxes)
function x = make_linear_chirp(fs, T, f0, f1)
N = round(T * fs);
t = (0:N-1).' / fs;
K = (f1 - f0) / T; % [Hz/s]
phi = 2*pi * (f0 * t + 0.5 * K * t.^2);
x = cos(phi);
end

function x = make_freq_hop(fs, T, fList, segT)
N = round(T * fs);
fInst = zeros(N, 1);
segN = round(segT * fs);

idx = 1;
for i = 1:numel(fList)
    idxEnd = min(N, idx + segN - 1);
    fInst(idx:idxEnd) = fList(i);
    idx = idxEnd + 1;
    if idx > N
        break;
    end
end
if idx <= N
    % Repeat last frequency if there is remaining tail.
    fInst(idx:N) = fList(min(end, numel(fList)));
end

phi = cumsum(2*pi * fInst / fs);
x = cos(phi);
end

