%% 第5章 周波数応答とDTFT：時間領域と周波数領域の対応
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
fs = 8000;   % [Hz]
N  = 2400;   % samples (0.3 s)
[t, n] = sp_time_axis(fs, N); %#ok<ASGLU>

%% 0) Quick checks: sp_freq_response must match the direct DTFT evaluation
rng(0);
b = randn(9, 1);
fCheck = linspace(0, fs/2, 1234).';
[fOut, H] = sp_freq_response(fs, b, "F", fCheck);

Omega = 2*pi*fOut/fs;
m = 0:numel(b)-1;
Href = exp(-1j * (Omega * m)) * b;

relErr = norm(H - Href) / max(norm(Href), 1);
fprintf("sp_freq_response check: relErr=%.3g\n", relErr);
assert(relErr < 1e-12, "sp_freq_response mismatch.");

%% 1) Complex exponential is (almost) an eigenfunction of an LTI system
M = 31;
hMA = ones(M, 1) / M;
f0 = 440; % [Hz]
Omega0 = 2*pi*f0/fs;
x = exp(1j*Omega0*n); % complex exponential

yFull = sp_convolve(x, hMA);
y = yFull(1:N);

[~, H0] = sp_freq_response(fs, hMA, "F", f0);
yPred = H0 * x;

idxSteady = (n >= (M-1)) & (n <= (M-1) + 500);
err = y(idxSteady) - yPred(idxSteady);
relErr = norm(err) / max(norm(yPred(idxSteady)), 1);
fprintf("complex exponential eigen check: relErr=%.3g\n", relErr);

idxPlot = (t >= 0.05) & (t <= 0.08);

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
plot(t(idxPlot), real(x(idxPlot)), "k--", "LineWidth", 1.0);
hold on;
plot(t(idxPlot), real(y(idxPlot)), "LineWidth", 1.2);
grid on;
xlabel("Time [s]"); ylabel("Real part");
title(sprintf("LTI response to a complex exponential (f0 = %.0f Hz, moving average M=%d)", f0, M));
legend(["Re{x}", "Re{y}"], "Location", "southeast");

nexttile;
plot(n(idxSteady), real(err), "LineWidth", 1.0);
hold on;
plot(n(idxSteady), imag(err), "LineWidth", 1.0);
grid on;
xlabel("n"); ylabel("error");
title(sprintf("y - H(f0) x  (steady region), relErr=%.3g", relErr));
legend(["Re", "Im"], "Location", "southeast");

exportgraphics(gcf, fullfile(figDir, "ch05_eigenfunction.png"), "Resolution", 200, "BackgroundColor", "white");

%% 2) "Time-short <-> Frequency-wide" intuition using moving averages
Mshort = 5;
Mlong = 51;
hShort = ones(Mshort, 1) / Mshort;
hLong = ones(Mlong, 1) / Mlong;

f = linspace(0, fs/2, 2500).';
[~, Hs] = sp_freq_response(fs, hShort, "F", f);
[~, Hl] = sp_freq_response(fs, hLong, "F", f);

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
hold on;
stem(0:Mshort-1, hShort, "filled", "MarkerSize", 2);
stem(0:Mlong-1, hLong, "filled", "MarkerSize", 2);
grid on;
xlabel("m"); ylabel("h[m]");
title("Impulse responses (short vs long)");
legend(["M=5", "M=51"], "Location", "northeast");

nexttile;
plot(f, 20*log10(abs(Hs) + 1e-12), "LineWidth", 1.2);
hold on;
plot(f, 20*log10(abs(Hl) + 1e-12), "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]"); ylabel("|H| [dB]");
title("Magnitude responses: longer averaging -> stronger low-pass");
legend(["M=5", "M=51"], "Location", "northeast");
ylim([-80, 5]);

exportgraphics(gcf, fullfile(figDir, "ch05_time_short_freq_wide.png"), "Resolution", 200, "BackgroundColor", "white");

%% 2b) DTFT periodicity: H repeats every fs (Hz) / 2*pi (rad/sample)
% NOTE: In discrete-time, "frequency in Hz" is only unique up to +/-fs/2.
% This plot intentionally shows repetition to build that intuition.
h = hLong;
fPer = linspace(-fs, fs, 4000).';
[~, Hper] = sp_freq_response(fs, h, "F", fPer);

figure;
plot(fPer, 20*log10(abs(Hper) + 1e-12), "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]"); ylabel("|H| [dB]");
title("DTFT periodicity in Hz: repeats every fs (aliasing of the frequency axis)");
xline(-fs/2, "k--", "LineWidth", 1.0);
xline(+fs/2, "k--", "LineWidth", 1.0);
xlim([-fs, fs]);
ylim([-80, 5]);

exportgraphics(gcf, fullfile(figDir, "ch05_dtft_periodicity.png"), "Resolution", 200, "BackgroundColor", "white");

%% 3) Frequency sampling: direct evaluation vs FFT samples (preview for Ch.6)
Nfft = 4096;
h = hLong; % use moving average as a clean example
Hfft = fft(h, Nfft);
fFft = (0:Nfft-1).' * (fs / Nfft);
[fDir, Hdir] = sp_freq_response(fs, h, "F", fFft);

magDiffDb = 20*log10(abs(Hdir) + 1e-15) - 20*log10(abs(Hfft) + 1e-15);
maxAbsDiffDb = max(abs(magDiffDb));
fprintf("DTFT vs FFT samples (moving average): max |magDiff| = %.3g dB\n", maxAbsDiffDb);

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
plot(fFft(1:Nfft/2+1), 20*log10(abs(Hdir(1:Nfft/2+1)) + 1e-12), "LineWidth", 1.2);
hold on;
plot(fFft(1:Nfft/2+1), 20*log10(abs(Hfft(1:Nfft/2+1)) + 1e-12), "--", "LineWidth", 1.0);
grid on;
xlabel("Frequency [Hz]"); ylabel("|H| [dB]");
title(sprintf("Direct evaluation vs FFT sampling (Nfft=%d)", Nfft));
legend(["direct (DTFT)", "fft(h,Nfft)"], "Location", "northeast");
ylim([-80, 5]);

nexttile;
plot(fFft(1:Nfft/2+1), magDiffDb(1:Nfft/2+1), "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]"); ylabel("difference [dB]");
title(sprintf("Magnitude difference (should be ~0), max=%.3g dB", maxAbsDiffDb));

exportgraphics(gcf, fullfile(figDir, "ch05_freq_sampling_vs_fft.png"), "Resolution", 200, "BackgroundColor", "white");

%% 4) Phase and pure delay: magnitude is flat, phase slope encodes delay
d = 12; % samples
hDelay = [zeros(d, 1); 1];

f = linspace(0, fs/2, 2500).';
[~, Hd] = sp_freq_response(fs, hDelay, "F", f);
[~, gd] = sp_group_delay(fs, hDelay, "F", f);

figure;
tiledlayout(3, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
plot(f, 20*log10(abs(Hd) + 1e-12), "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]"); ylabel("|H| [dB]");
title(sprintf("Pure delay: magnitude is flat (d=%d samples)", d));
ylim([-1, 1]);

nexttile;
plot(f, unwrap(angle(Hd)), "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]"); ylabel("phase [rad]");
title("Unwrapped phase is (almost) a straight line");

nexttile;
plot(f, gd, "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]"); ylabel("group delay [samples]");
title("Group delay is constant");
ylim([d-2, d+2]);

exportgraphics(gcf, fullfile(figDir, "ch05_delay_phase_gd.png"), "Resolution", 200, "BackgroundColor", "white");

%% 4b) Phase wrapping vs unwrapping (why unwrap matters)
f = linspace(0, fs/2, 2500).';
[~, Hd] = sp_freq_response(fs, hDelay, "F", f);
phiWrapped = angle(Hd);
phiUnwrapped = unwrap(phiWrapped);

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
plot(f, phiWrapped, "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]"); ylabel("angle(H) [rad]");
title("Wrapped phase (principal value in (-pi, pi])");

nexttile;
plot(f, phiUnwrapped, "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]"); ylabel("unwrap(angle(H)) [rad]");
title("Unwrapped phase (continuous slope is visible)");

exportgraphics(gcf, fullfile(figDir, "ch05_phase_wrap_unwrap.png"), "Resolution", 200, "BackgroundColor", "white");

%% 5) Linear-phase vs non-linear-phase: symmetry matters
hSym = [0.1; 0.2; 0.4; 0.2; 0.1];
hAsym = [0.4; 0.2; 0.1; 0.1; 0.2];

f = linspace(0, fs/2, 2500).';
[~, Hsym] = sp_freq_response(fs, hSym, "F", f);
[~, Hasym] = sp_freq_response(fs, hAsym, "F", f);
[~, gdSym] = sp_group_delay(fs, hSym, "F", f);
[~, gdAsym] = sp_group_delay(fs, hAsym, "F", f);

figure;
tiledlayout(3, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
plot(f, 20*log10(abs(Hsym) + 1e-12), "LineWidth", 1.2);
hold on;
plot(f, 20*log10(abs(Hasym) + 1e-12), "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]"); ylabel("|H| [dB]");
title("Magnitude responses");
legend(["symmetric (linear-phase)", "asymmetric"], "Location", "northeast");
ylim([-60, 10]);

nexttile;
plot(f, unwrap(angle(Hsym)), "LineWidth", 1.2);
hold on;
plot(f, unwrap(angle(Hasym)), "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]"); ylabel("phase [rad]");
title("Phase responses (unwrapped)");

nexttile;
plot(f, gdSym, "LineWidth", 1.2);
hold on;
plot(f, gdAsym, "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]"); ylabel("group delay [samples]");
title("Group delay: symmetric FIR -> almost constant");
ylim([0, 6]);

exportgraphics(gcf, fullfile(figDir, "ch05_linear_phase_vs_nonlinear.png"), "Resolution", 200, "BackgroundColor", "white");

%% 5b) Difference filter: high-pass / "derivative" intuition
hDiff = [1; -1];
hAvg2 = [0.5; 0.5];

f = linspace(0, fs/2, 2500).';
[~, Hdiff] = sp_freq_response(fs, hDiff, "F", f);
[~, Havg2] = sp_freq_response(fs, hAvg2, "F", f);
[~, gdDiff] = sp_group_delay(fs, hDiff, "F", f, "MinMag", 1e-3);

figure;
tiledlayout(3, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
hold on;
stem(0:numel(hDiff)-1, hDiff, "filled", "MarkerSize", 3);
stem(0:numel(hAvg2)-1, hAvg2, "filled", "MarkerSize", 3);
grid on;
xlabel("m"); ylabel("h[m]");
title("2-tap FIR examples: difference vs moving average");
legend(["difference [1, -1]", "average [1/2, 1/2]"], "Location", "northeast");

nexttile;
plot(f, 20*log10(abs(Hdiff) + 1e-12), "LineWidth", 1.2);
hold on;
plot(f, 20*log10(abs(Havg2) + 1e-12), "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]"); ylabel("|H| [dB]");
title("Magnitude responses (difference is high-pass)");
legend(["difference", "average"], "Location", "northeast");
ylim([-80, 10]);

nexttile;
plot(f, gdDiff, "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]"); ylabel("group delay [samples]");
title("Group delay of the difference filter (undefined near |H|~0)");
ylim([0, 2]);

exportgraphics(gcf, fullfile(figDir, "ch05_difference_filter.png"), "Resolution", 200, "BackgroundColor", "white");

%% 6) Same magnitude, different phase: reversing h keeps |H| but changes the waveform
rng(1);
h0 = randn(9, 1);
h0 = h0 / sum(abs(h0));
h1 = flipud(h0);

f = linspace(0, fs/2, 2500).';
[~, H0] = sp_freq_response(fs, h0, "F", f);
[~, H1] = sp_freq_response(fs, h1, "F", f);
magGap = max(abs(abs(H0) - abs(H1)));
fprintf("reverse-h magnitude gap: %.3g\n", magGap);

% input: short burst to make time-domain differences visible
tb = 0.12;
burst = (t >= 0.06) & (t <= 0.06 + tb);
xBurst = zeros(N, 1);
L = nnz(burst);
k = (0:L-1).';
wHann = 0.5 - 0.5*cos(2*pi*k/(L-1));
xBurst(burst) = sin(2*pi*440*t(burst)) .* wHann;

y0 = conv(xBurst, h0, "same");
y1 = conv(xBurst, h1, "same");

idx = (t >= 0.04) & (t <= 0.20);

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
plot(f, 20*log10(abs(H0) + 1e-12), "LineWidth", 1.2);
hold on;
plot(f, 20*log10(abs(H1) + 1e-12), "--", "LineWidth", 1.0);
grid on;
xlabel("Frequency [Hz]"); ylabel("|H| [dB]");
title(sprintf("Same magnitude (h and flipud(h)), max |Δ| = %.3g", magGap));
legend(["h", "flipud(h)"], "Location", "northeast");

nexttile;
plot(t(idx), xBurst(idx), "k--", "LineWidth", 1.0);
hold on;
plot(t(idx), y0(idx), "LineWidth", 1.2);
plot(t(idx), y1(idx), "LineWidth", 1.2);
grid on;
xlabel("Time [s]"); ylabel("Amplitude");
title("Time-domain outputs can differ even when |H| is the same");
legend(["x (burst)", "y (h)", "y (flipud(h))"], "Location", "southeast");

exportgraphics(gcf, fullfile(figDir, "ch05_same_mag_diff_phase.png"), "Resolution", 200, "BackgroundColor", "white");

%% 7) Finite-length observation: shorter time -> wider spectrum (preview for Ch.7)
f0 = 440;
Nshort = 256;
Nlong = 2048;
Nfft = 16384;

nS = (0:Nshort-1).';
nL = (0:Nlong-1).';
xS = sin(2*pi*f0*nS/fs);
xL = sin(2*pi*f0*nL/fs);

XS = fft(xS, Nfft);
XL = fft(xL, Nfft);

kMax = floor(Nfft/2);
fAxis = (0:kMax).' * (fs / Nfft);
maxAbsXS = max(abs(XS));
maxAbsXL = max(abs(XL));
maxAbsXS = max(maxAbsXS, 1e-12);
maxAbsXL = max(maxAbsXL, 1e-12);

magSdB = 20*log10(abs(XS(1:kMax+1)) ./ maxAbsXS + 1e-12);
magLdB = 20*log10(abs(XL(1:kMax+1)) ./ maxAbsXL + 1e-12);

idx = (fAxis >= f0 - 250) & (fAxis <= f0 + 250);

figure;
plot(fAxis(idx), magSdB(idx), "LineWidth", 1.2);
hold on;
plot(fAxis(idx), magLdB(idx), "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]"); ylabel("Normalized magnitude [dB]");
title("Finite-length tone: shorter observation -> wider spectrum (rectangular window)");
legend([sprintf("N=%d", Nshort), sprintf("N=%d", Nlong)], "Location", "southeast");
xlim([200, 700]);
ylim([-80, 5]);

exportgraphics(gcf, fullfile(figDir, "ch05_finite_length_tone.png"), "Resolution", 200, "BackgroundColor", "white");

%% 8) IIR frequency response (preview): first-order smoothing is not linear-phase
alist = [0.2, 0.8, 0.95];
f = linspace(0, fs/2, 2500).';

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
hold on;
for k = 1:numel(alist)
    a = alist(k);
    b = (1-a);
    denom = [1; -a];
    [~, H] = sp_freq_response(fs, b, "a", denom, "F", f);
    plot(f, 20*log10(abs(H) + 1e-12), "LineWidth", 1.2);
end
grid on;
xlabel("Frequency [Hz]"); ylabel("|H| [dB]");
title("First-order IIR magnitude response  H(z)=(1-a)/(1-a z^{-1})");
legend(["a=0.2", "a=0.8", "a=0.95"], "Location", "northeast");
ylim([-60, 5]);

nexttile;
hold on;
for k = 1:numel(alist)
    a = alist(k);
    b = (1-a);
    denom = [1; -a];
    [~, gd] = sp_group_delay(fs, b, "a", denom, "F", f, "MinMag", 1e-4);
    plot(f, gd, "LineWidth", 1.2);
end
grid on;
xlabel("Frequency [Hz]"); ylabel("group delay [samples]");
title("Group delay is frequency-dependent (waveform can distort)");
legend(["a=0.2", "a=0.8", "a=0.95"], "Location", "northeast");

exportgraphics(gcf, fullfile(figDir, "ch05_iir1_freqresp.png"), "Resolution", 200, "BackgroundColor", "white");
