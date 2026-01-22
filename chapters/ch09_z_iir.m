%% 第9章 Z変換・極零・安定性：IIRフィルタの基礎
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
fs = 8000;  % [Hz]

%% 1) Stability in time domain (1st-order IIR impulse response)
N = 200;
x = zeros(N, 1);
x(1) = 1; % impulse

aStable = 0.9;
aUnstable = 1.01; % intentionally unstable (|a|>1)
yStable = sp_iir1(x, aStable);
yUnstable = sp_iir1(x, aUnstable);

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
stem(0:N-1, yStable, "filled");
grid on;
xlabel("n");
ylabel("h[n]");
title(sprintf("1st-order IIR impulse response (stable): a=%.3f", aStable));

nexttile;
stem(0:N-1, yUnstable, "filled");
grid on;
xlabel("n");
ylabel("h[n]");
title(sprintf("1st-order IIR impulse response (unstable): a=%.3f", aUnstable));

exportgraphics(gcf, fullfile(figDir, "ch09_iir1_impulse_stability.png"), "Resolution", 200, "BackgroundColor", "white");

%% 2) Complex poles: ring-down gets longer as |p| -> 1
theta = 0.20 * pi; % rad/sample (peak around theta*fs/(2pi) Hz)
rs = [0.8, 0.95, 0.99];

N = 400;
x = zeros(N, 1);
x(1) = 1;

figure;
tiledlayout(numel(rs), 1, "Padding", "compact", "TileSpacing", "compact");
for i = 1:numel(rs)
    r = rs(i);
    a = [1, -2*r*cos(theta), r^2];
    b = [1, 0, 0]; % all-pole resonator
    y = sp_iir2(x, b, a, "Structure", "df1");

    nexttile;
    plot(0:N-1, y, "LineWidth", 1.2);
    grid on;
    xlabel("n");
    ylabel("h[n]");
    title(sprintf("2nd-order ring-down: r=%.2f, theta=%.2f pi", r, theta/pi));
end
exportgraphics(gcf, fullfile(figDir, "ch09_pole_radius_ringdown.png"), "Resolution", 200, "BackgroundColor", "white");

%% 3) Pole/zero plot and magnitude response (resonator)
r = 0.98;
a = [1, -2*r*cos(theta), r^2];
b = [1, 0, 0];

[f, H] = sp_freq_response(fs, b, "a", a, "N", 2000, "Range", "half");

% Make the figure wide enough so the z-plane plot (axis equal) does not
% shrink vertically within a 1x2 tiledlayout.
figure("Units", "pixels", "Position", [100 100 980 420]);
tiledlayout(1, 2, "Padding", "compact", "TileSpacing", "compact");

nexttile;
sp_zplane_plot(b, a, "NewFigure", false, "Title", "Resonator (poles/zeros)", "AxisLimit", 1.4);

nexttile;
plot(f, 20*log10(abs(H) + eps), "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]");
ylabel("Magnitude [dB]");
title("Magnitude response |H(e^{j\Omega})|");
xlim([0, fs/2]);

exportgraphics(gcf, fullfile(figDir, "ch09_resonator_pz_mag.png"), "Resolution", 200, "BackgroundColor", "white");

%% 4) Notch filter: zeros on the unit circle, poles slightly inside
theta0 = 0.25 * pi;
r0 = 0.98;
b = [1, -2*cos(theta0), 1];
a = [1, -2*r0*cos(theta0), r0^2];

[f, H] = sp_freq_response(fs, b, "a", a, "N", 2000, "Range", "half");

% Make the figure wide enough so the z-plane plot (axis equal) does not
% shrink vertically within a 1x2 tiledlayout.
figure("Units", "pixels", "Position", [100 100 980 420]);
tiledlayout(1, 2, "Padding", "compact", "TileSpacing", "compact");

nexttile;
sp_zplane_plot(b, a, "NewFigure", false, "Title", "Notch (poles/zeros)", "AxisLimit", 1.4);

nexttile;
plot(f, 20*log10(abs(H) + eps), "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]");
ylabel("Magnitude [dB]");
title("Notch magnitude response");
xlim([0, fs/2]);

exportgraphics(gcf, fullfile(figDir, "ch09_notch_pz_mag.png"), "Resolution", 200, "BackgroundColor", "white");

%% 5) All-pass example: flat magnitude, non-trivial group delay
theta1 = 0.30 * pi;
r1 = 0.90;
a = [1, -2*r1*cos(theta1), r1^2];
b = [a(3), a(2), a(1)]; % flip(a) for real all-pass (2nd-order)

[f, H] = sp_freq_response(fs, b, "a", a, "N", 2000, "Range", "half");
[~, gd] = sp_group_delay(fs, b, "a", a, "F", f, "Range", "half");

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
plot(f, 20*log10(abs(H) + eps), "LineWidth", 1.2);
grid on;
ylabel("Magnitude [dB]");
title("All-pass magnitude (should be ~0 dB)");
xlim([0, fs/2]);
% Force a readable scale around 0 dB to avoid scientific-notation offsets.
ylim([-0.2, 0.2]);
ax = gca;
ax.YAxis.Exponent = 0;

nexttile;
plot(f, gd, "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]");
ylabel("Group delay [samples]");
title("All-pass group delay");
xlim([0, fs/2]);

exportgraphics(gcf, fullfile(figDir, "ch09_allpass_gd.png"), "Resolution", 200, "BackgroundColor", "white");

%% 6) Coefficient quantization moves poles/zeros
thetaQ = 0.25 * pi;
rQ = 0.995;  % near-unit poles: sensitive to coefficient perturbations
b = [1, -2*cos(thetaQ), 1];
a = [1, -2*rQ*cos(thetaQ), rQ^2];

bits = 10;
fullScale = 2;
bq = sp_quantize_uniform(b(:), bits, "FullScale", fullScale);
aq12 = sp_quantize_uniform(a(2:3).', bits, "FullScale", fullScale);
aq = a(:);
aq(2:3) = aq12(:);

z = roots(b);
p = roots(a);
zq = roots(bq);
pq = roots(aq);

[f, H] = sp_freq_response(fs, b, "a", a, "N", 2000, "Range", "half");
[~, Hq] = sp_freq_response(fs, bq, "a", aq, "F", f, "Range", "half");

% Make the figure wide enough so the z-plane plot (axis equal) does not
% shrink vertically within a 1x2 tiledlayout.
figure("Units", "pixels", "Position", [100 100 980 420]);
tiledlayout(1, 2, "Padding", "compact", "TileSpacing", "compact");

nexttile;
hold on;
th = linspace(0, 2*pi, 361);
plot(cos(th), sin(th), "k-", "LineWidth", 1.0);
plot(real(z), imag(z), "bo", "LineWidth", 1.2);
plot(real(p), imag(p), "bx", "LineWidth", 1.2);
plot(real(zq), imag(zq), "ro", "LineWidth", 1.2);
plot(real(pq), imag(pq), "rx", "LineWidth", 1.2);
axis equal;
grid on;
xlim([-1.4 1.4]); ylim([-1.4 1.4]);
xlabel("Re\{z\}"); ylabel("Im\{z\}");
title(sprintf("Pole/zero shift by quantization (%d bits)", bits));
% Keep legend inside to avoid shrinking the plot area in tiledlayout.
legend(["unit circle", "zeros (orig)", "poles (orig)", "zeros (quant)", "poles (quant)"], "Location", "northwest");
hold off;

nexttile;
plot(f, 20*log10(abs(H) + eps), "LineWidth", 1.2); hold on;
plot(f, 20*log10(abs(Hq) + eps), "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]");
ylabel("Magnitude [dB]");
title("Magnitude response: original vs quantized");
legend(["original", "quantized"], "Location", "best");
xlim([0, fs/2]);
hold off;

exportgraphics(gcf, fullfile(figDir, "ch09_coeff_quantization_effect.png"), "Resolution", 200, "BackgroundColor", "white");

%% 7) Structure and finite precision: df1 vs df2 in single precision
rS = 0.9995;
thetaS = 0.30 * pi;
b = [1, 0, 0];
a = [1, -2*rS*cos(thetaS), rS^2];

rng(0);
N = 6000;
x = randn(N, 1);

yRef = sp_iir2(x, b, a, "Structure", "df1", "Arithmetic", "double");
yDf1 = sp_iir2(x, b, a, "Structure", "df1", "Arithmetic", "single");
yDf2 = sp_iir2(x, b, a, "Structure", "df2", "Arithmetic", "single");

e1 = double(yDf1) - double(yRef);
e2 = double(yDf2) - double(yRef);

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
plot(yRef(1:800), "k-", "LineWidth", 1.0); hold on;
plot(double(yDf1(1:800)), "b--", "LineWidth", 1.0);
plot(double(yDf2(1:800)), "r--", "LineWidth", 1.0);
grid on;
xlabel("n");
ylabel("y[n]");
title("Output waveform (reference vs single-precision)");
legend(["double ref", "df1 single", "df2 single"], "Location", "best");
hold off;

nexttile;
plot(e1(1:800), "b-", "LineWidth", 1.0); hold on;
plot(e2(1:800), "r-", "LineWidth", 1.0);
grid on;
xlabel("n");
ylabel("error");
title("Error relative to double");
legend(["df1 single - ref", "df2 single - ref"], "Location", "best");
hold off;

exportgraphics(gcf, fullfile(figDir, "ch09_structure_numeric_error.png"), "Resolution", 200, "BackgroundColor", "white");

%% 8) Frequency response: direct unit-circle evaluation vs FFT of truncated impulse response
rF = 0.99;
thetaF = 0.20 * pi;
b = [1, 0, 0];
a = [1, -2*rF*cos(thetaF), rF^2];

M = 256;           % truncation length of the impulse response
Nfft = 4096;       % FFT length (zero padding for a smooth curve)
x = zeros(M, 1);
x(1) = 1;
h = sp_iir2(x, b, a, "Structure", "df1");

Hfft = fft(h, Nfft);
kMax = floor(Nfft/2);
fFft = (0:kMax).' * (fs / Nfft);
Hfft = Hfft(1:kMax+1);

[fDir, Hdir] = sp_freq_response(fs, b, "a", a, "F", fFft, "Range", "half");
assert(max(abs(fDir - fFft)) < 1e-12);

figure;
plot(fDir, 20*log10(abs(Hdir) + eps), "LineWidth", 1.2); hold on;
plot(fDir, 20*log10(abs(Hfft) + eps), "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]");
ylabel("Magnitude [dB]");
title(sprintf("Direct H(e^{j\\Omega}) vs FFT of truncated h[n] (M=%d, Nfft=%d)", M, Nfft));
legend(["direct unit-circle eval", "FFT(truncated h)"], "Location", "best");
xlim([0, fs/2]);
hold off;

exportgraphics(gcf, fullfile(figDir, "ch09_freqresp_direct_vs_fft.png"), "Resolution", 200, "BackgroundColor", "white");
