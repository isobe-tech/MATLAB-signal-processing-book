%% 第10章 FIRフィルタ設計（窓法）と総合プロジェクト
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

%% Parameters (example spec)
fs = 8000;                 % [Hz]
Fp = 800;                  % passband edge [Hz]
Fs = 1200;                 % stopband edge [Hz]
N  = 101;                  % taps (odd)
Fc = (Fp + Fs) / 2;        % design cutoff placed at the middle of transition band [Hz]

%% 1) Ideal LPF impulse response (sinc slice)
hIdeal = sp_fir_ideal_lpf(fs, Fc, N);
n = (0:N-1).';
alpha = (N-1)/2;

figure;
plot(n - alpha, hIdeal, "LineWidth", 1.2);
grid on;
xlabel("n-\alpha");
ylabel("h_d[n]");
title(sprintf("Ideal LPF impulse response slice (sinc): Fc=%d Hz, N=%d", Fc, N));
exportgraphics(gcf, fullfile(figDir, "ch10_ideal_lpf_sinc.png"), "Resolution", 200, "BackgroundColor", "white");

%% 2) Windowed FIR coefficients (rect / Hann / Hamming)
[hRect, infoR] = sp_fir_window_design(fs, "lpf", Fc, N, "Window", "rect");
[hHann, infoHn] = sp_fir_window_design(fs, "lpf", Fc, N, "Window", "hann");
[hHamm, infoHm] = sp_fir_window_design(fs, "lpf", Fc, N, "Window", "hamming");

figure;
plot(n - alpha, hRect, "LineWidth", 1.0); hold on;
plot(n - alpha, hHann, "LineWidth", 1.0);
plot(n - alpha, hHamm, "LineWidth", 1.0);
grid on;
xlabel("n-\alpha");
ylabel("h[n]");
title(sprintf("Windowed FIR coefficients (LPF): N=%d, Fc=%d Hz (spec: Fp=%d, Fs=%d)", N, Fc, Fp, Fs));
legend([infoR.windowName, infoHn.windowName, infoHm.windowName], "Location", "best");
hold off;
exportgraphics(gcf, fullfile(figDir, "ch10_coeffs_windows.png"), "Resolution", 200, "BackgroundColor", "white");

%% 3) Step response (ringing / Gibbs) for different windows
Nstep = 600;
xStep = ones(Nstep, 1);

yRectStep = filter(hRect, 1, xStep);
yHannStep = filter(hHann, 1, xStep);
yHammStep = filter(hHamm, 1, xStep);

d = (N-1)/2; % group delay [samples]
yRectA = yRectStep(d+1:end);
yHannA = yHannStep(d+1:end);
yHammA = yHammStep(d+1:end);
nA = (0:numel(yRectA)-1).';

figure;
plot(nA, yRectA, "LineWidth", 1.1); hold on;
plot(nA, yHannA, "LineWidth", 1.1);
plot(nA, yHammA, "LineWidth", 1.1);
yline(1, "k:", "ideal (1)");
grid on;
xlabel("n (delay-compensated)");
ylabel("y[n]");
title(sprintf("Step response (ringing): N=%d, Fc=%d Hz (spec: Fp=%d, Fs=%d)", N, Fc, Fp, Fs));
legend(["rect", "hann", "hamming"], "Location", "best");
xlim([0, 80]);
hold off;
exportgraphics(gcf, fullfile(figDir, "ch10_step_response_windows.png"), "Resolution", 200, "BackgroundColor", "white");

%% 4) Window comparison in magnitude response
[f, Hrect] = sp_freq_response(fs, hRect, "N", 4096, "Range", "half");
[~, Hhann] = sp_freq_response(fs, hHann, "F", f, "Range", "half");
[~, Hhamm] = sp_freq_response(fs, hHamm, "F", f, "Range", "half");

figure;
plot(f, 20*log10(abs(Hrect) + eps), "LineWidth", 1.1); hold on;
plot(f, 20*log10(abs(Hhann) + eps), "LineWidth", 1.1);
plot(f, 20*log10(abs(Hhamm) + eps), "LineWidth", 1.1);
grid on;
xlabel("Frequency [Hz]");
ylabel("Magnitude [dB]");
title(sprintf("Window comparison (LPF): N=%d, Fc=%d Hz (spec: Fp=%d, Fs=%d)", N, Fc, Fp, Fs));
xlim([0, fs/2]);
ylim([-120, 5]);
xline(Fp, "k--", "Fp");
xline(Fs, "k--", "Fs");
legend(["rect", "hann", "hamming"], "Location", "southwest");
hold off;
exportgraphics(gcf, fullfile(figDir, "ch10_window_compare_mag.png"), "Resolution", 200, "BackgroundColor", "white");

%% 5) Tap count effect (transition width vs N)
Ns = [31, 63, 127];
figure;
for i = 1:numel(Ns)
    Ni = Ns(i);
    [hi, ~] = sp_fir_window_design(fs, "lpf", Fc, Ni, "Window", "hann");
    [fi, Hi] = sp_freq_response(fs, hi, "N", 4096, "Range", "half");
    plot(fi, 20*log10(abs(Hi) + eps), "LineWidth", 1.1); hold on;
end
grid on;
xlabel("Frequency [Hz]");
ylabel("Magnitude [dB]");
title(sprintf("Tap count effect (Hann window): Fc=%d Hz (spec: Fp=%d, Fs=%d)", Fc, Fp, Fs));
xlim([0, fs/2]);
ylim([-120, 5]);
xline(Fp, "k--", "Fp");
xline(Fs, "k--", "Fs");
legend("N=31","N=63","N=127", "Location", "southwest");
hold off;
exportgraphics(gcf, fullfile(figDir, "ch10_tap_count_effect.png"), "Resolution", 200, "BackgroundColor", "white");

%% 6) Linear phase and group delay (passband zoom)
% Phase/group delay are most meaningful in the passband; deep stopbands can
% make phase unwrapping and numerical derivatives look noisy.
[fAll, Hall] = sp_freq_response(fs, hHann, "N", 4096, "Range", "half");
maskPb = fAll <= Fp;
fPb = fAll(maskPb);
HPb = Hall(maskPb);
[~, gd] = sp_group_delay(fs, hHann, "F", fPb, "Range", "half");

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
plot(fPb, unwrap(angle(HPb)), "LineWidth", 1.1);
grid on;
ylabel("Phase [rad]");
title("Unwrapped phase of a linear-phase FIR (passband zoom)");
xlim([0, Fp]);

nexttile;
plot(fPb, gd, "LineWidth", 1.1);
grid on;
xlabel("Frequency [Hz]");
ylabel("Group delay [samples]");
title(sprintf("Group delay (should be ~constant): (N-1)/2 = %.1f samples", (N-1)/2));
xlim([0, Fp]);
gd0 = (N-1)/2;
ylim([gd0-6, gd0+6]);
exportgraphics(gcf, fullfile(figDir, "ch10_linear_phase_gd.png"), "Resolution", 200, "BackgroundColor", "white");

%% 7) FIR frequency response: direct evaluation vs FFT (for FIR they match)
Nfft = 4096;
Hfft = fft(hHann, Nfft);
kMax = floor(Nfft/2);
fFft = (0:kMax).' * (fs/Nfft);
Hfft = Hfft(1:kMax+1);
[~, Hdir] = sp_freq_response(fs, hHann, "F", fFft, "Range", "half");

figure;
plot(fFft, 20*log10(abs(Hdir) + eps), "LineWidth", 1.2); hold on;
plot(fFft, 20*log10(abs(Hfft) + eps), "--", "LineWidth", 1.0);
grid on;
xlabel("Frequency [Hz]");
ylabel("Magnitude [dB]");
title("FIR: direct unit-circle evaluation equals FFT of h[n] (up to numerical error)");
legend(["direct (sp\_freq\_response)", "FFT(h)"], "Location", "southwest");
xlim([0, fs/2]);
ylim([-120, 5]);
hold off;
exportgraphics(gcf, fullfile(figDir, "ch10_fir_direct_vs_fft.png"), "Resolution", 200, "BackgroundColor", "white");

%% 8) Spec evaluation example (ripple/atten) for window comparison
pb = [0, Fp];
sb = [Fs, fs/2];
sRect = sp_filter_specs_eval(f, Hrect, "Passband", pb, "Stopband", sb);
sHann = sp_filter_specs_eval(f, Hhann, "Passband", pb, "Stopband", sb);
sHamm = sp_filter_specs_eval(f, Hhamm, "Passband", pb, "Stopband", sb);

figure;
bar([sRect.stopbandAttenDb, sHann.stopbandAttenDb, sHamm.stopbandAttenDb]);
grid on;
xticklabels(["rect","hann","hamming"]);
ylabel("Stopband attenuation [dB]");
title(sprintf("Measured stopband attenuation (N=%d, Fp=%d, Fs=%d)", N, Fp, Fs));
exportgraphics(gcf, fullfile(figDir, "ch10_specs_atten_compare.png"), "Resolution", 200, "BackgroundColor", "white");

%% 9) Least squares FIR design (A\\d): comparison with window method
[hLs, ~] = sp_fir_ls_lpf(fs, Fp, Fs, N, "GridN", 4096, "Wp", 1, "Ws", 80, "Scale", true);
[~, Hls] = sp_freq_response(fs, hLs, "F", f, "Range", "half");
sWin = sp_filter_specs_eval(f, Hhann, "Passband", pb, "Stopband", sb);
sLs  = sp_filter_specs_eval(f, Hls,  "Passband", pb, "Stopband", sb);

figure;
plot(f, 20*log10(abs(Hhann) + eps), "LineWidth", 1.1); hold on;
plot(f, 20*log10(abs(Hls) + eps), "LineWidth", 1.1);
grid on;
xlabel("Frequency [Hz]");
ylabel("Magnitude [dB]");
title(sprintf("Window(Hann) vs LS: N=%d (stopband: %.1f dB vs %.1f dB)", N, sWin.stopbandAttenDb, sLs.stopbandAttenDb));
xlim([0, fs/2]);
ylim([-120, 5]);
xline(Fp, "k--", "Fp");
xline(Fs, "k--", "Fs");
legend(["window (hann)", "least squares (Ws=80)"], "Location", "southwest");
hold off;
exportgraphics(gcf, fullfile(figDir, "ch10_window_vs_ls_mag.png"), "Resolution", 200, "BackgroundColor", "white");

%% 10) LS weight tradeoff (passband vs stopband)
WsList = [10, 30, 100];
HlsAll = zeros(numel(f), numel(WsList));
for i = 1:numel(WsList)
    [hTmp, ~] = sp_fir_ls_lpf(fs, Fp, Fs, N, "GridN", 4096, "Wp", 1, "Ws", WsList(i), "Scale", true);
    [~, Htmp] = sp_freq_response(fs, hTmp, "F", f, "Range", "half");
    HlsAll(:, i) = Htmp;
end

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
for i = 1:numel(WsList)
    plot(f, 20*log10(abs(HlsAll(:, i)) + eps), "LineWidth", 1.1); hold on;
end
grid on;
ylabel("Magnitude [dB]");
title("Passband zoom");
xlim([0, Fp]);
ylim([-1, 1]);
legend("Ws=10","Ws=30","Ws=100", "Location", "best");
hold off;

nexttile;
for i = 1:numel(WsList)
    plot(f, 20*log10(abs(HlsAll(:, i)) + eps), "LineWidth", 1.1); hold on;
end
grid on;
xlabel("Frequency [Hz]");
ylabel("Magnitude [dB]");
title("Stopband zoom");
xlim([Fs, fs/2]);
ylim([-120, 5]);
hold off;

exportgraphics(gcf, fullfile(figDir, "ch10_ls_weight_tradeoff.png"), "Resolution", 200, "BackgroundColor", "white");

%% 11) LPF/HPF/BPF/BSF examples (same N/window)
Nex = 101;
win = "hann";

[hL, ~] = sp_fir_window_design(fs, "lpf", 1000, Nex, "Window", win);
[hH, ~] = sp_fir_window_design(fs, "hpf", 1000, Nex, "Window", win);
[hB, ~] = sp_fir_window_design(fs, "bpf", [600, 1400], Nex, "Window", win);
[hS, ~] = sp_fir_window_design(fs, "bsf", [600, 1400], Nex, "Window", win);

[fEx, HL] = sp_freq_response(fs, hL, "N", 4096, "Range", "half");
[~,  HH] = sp_freq_response(fs, hH, "F", fEx, "Range", "half");
[~,  HB] = sp_freq_response(fs, hB, "F", fEx, "Range", "half");
[~,  HS] = sp_freq_response(fs, hS, "F", fEx, "Range", "half");

figure;
tiledlayout(2, 2, "Padding", "compact", "TileSpacing", "compact");

nexttile;
plot(fEx, 20*log10(abs(HL) + eps), "LineWidth", 1.1); grid on;
xlabel("Frequency [Hz]"); ylabel("Magnitude [dB]"); title("LPF");
xlim([0, fs/2]); ylim([-120, 5]);

nexttile;
plot(fEx, 20*log10(abs(HH) + eps), "LineWidth", 1.1); grid on;
xlabel("Frequency [Hz]"); ylabel("Magnitude [dB]"); title("HPF");
xlim([0, fs/2]); ylim([-120, 5]);

nexttile;
plot(fEx, 20*log10(abs(HB) + eps), "LineWidth", 1.1); grid on;
xlabel("Frequency [Hz]"); ylabel("Magnitude [dB]"); title("BPF");
xlim([0, fs/2]); ylim([-120, 5]);

nexttile;
plot(fEx, 20*log10(abs(HS) + eps), "LineWidth", 1.1); grid on;
xlabel("Frequency [Hz]"); ylabel("Magnitude [dB]"); title("BSF");
xlim([0, fs/2]); ylim([-120, 5]);

exportgraphics(gcf, fullfile(figDir, "ch10_band_filters_mag.png"), "Resolution", 200, "BackgroundColor", "white");

%% 12) Mini project: band extraction (low-pass) on a synthetic non-stationary signal
T = 2.0;            % [s]
Nsamp = round(T*fs);
[t, ~] = sp_time_axis(fs, Nsamp);

% Wanted component: 300 Hz tone that changes amplitude halfway
xWanted = sin(2*pi*300*t) .* (t < 1.0) + 0.7*sin(2*pi*300*t) .* (t >= 1.0);
% Interference: high-frequency tone + white noise
rng(0);
xInterf = 0.6*sin(2*pi*2200*t) + 0.25*randn(Nsamp, 1);
x = xWanted + xInterf;

y = filter(hHann, 1, x); % FIR filtering (allowed)

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");
nexttile;
plot(t, x, "LineWidth", 1.0);
grid on;
xlabel("Time [s]");
ylabel("Amplitude");
title("Input signal (wanted + interference)");
nexttile;
plot(t, y, "LineWidth", 1.0);
grid on;
xlabel("Time [s]");
ylabel("Amplitude");
title("Output after FIR low-pass (window method)");
exportgraphics(gcf, fullfile(figDir, "ch10_project_time_before_after.png"), "Resolution", 200, "BackgroundColor", "white");

% PSD before/after (reuse Chapter 7 scaling)
[fP, Pxx] = sp_psd_avg(x, fs, "FrameLength", 512, "Hop", 256, "Window", "hann", "Nfft", 2048, "Range", "half");
[~,  Pyy] = sp_psd_avg(y, fs, "FrameLength", 512, "Hop", 256, "Window", "hann", "Nfft", 2048, "Range", "half");

figure;
plot(fP, 10*log10(Pxx + 1e-20), "LineWidth", 1.1); hold on;
plot(fP, 10*log10(Pyy + 1e-20), "LineWidth", 1.1);
grid on;
xlabel("Frequency [Hz]");
ylabel("PSD [dB/Hz]");
title("PSD before/after filtering (Welch-like average)");
legend(["input", "output"], "Location", "best");
xlim([0, fs/2]);
hold off;
exportgraphics(gcf, fullfile(figDir, "ch10_project_psd_before_after.png"), "Resolution", 200, "BackgroundColor", "white");

% Spectrogram (optional: reuse Chapter 8)
[tS, fS, Pst] = sp_stft(x, fs, "FrameLength", 256, "Hop", 128, "Window", "hann", "Nfft", 1024, ...
    "Range", "half", "PadEnd", false, "TimeReference", "center", "Output", "psd");
[~,  ~, PstY] = sp_stft(y, fs, "FrameLength", 256, "Hop", 128, "Window", "hann", "Nfft", 1024, ...
    "Range", "half", "PadEnd", false, "TimeReference", "center", "Output", "psd");

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");
nexttile;
sp_spectrogram_plot(tS, fS, Pst, "Title", "Input spectrogram (PSD)", "DynamicRangeDb", 80, "NewFigure", false);
nexttile;
sp_spectrogram_plot(tS, fS, PstY, "Title", "Output spectrogram (PSD)", "DynamicRangeDb", 80, "NewFigure", false);
exportgraphics(gcf, fullfile(figDir, "ch10_project_spectrogram_before_after.png"), "Resolution", 200, "BackgroundColor", "white");
