%% 第2章 離散時間信号：サンプリング・量子化・時間軸生成
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

%% Sampling: fs を変えると「何が変わるか」
f0 = 440;             % [Hz]
tObs = 0.02;          % [s] 20 ms
fsList = [8000, 2000, 800];

tFine = linspace(0, tObs, 2500).';
xFine = sin(2*pi*f0*tFine);

figure;
tiledlayout(numel(fsList), 1, "Padding", "compact", "TileSpacing", "compact");
for i = 1:numel(fsList)
    fs = fsList(i);
    N = max(1, round(tObs * fs));
    [t, ~] = sp_time_axis(fs, N);
    x = sin(2*pi*f0*t);

    nexttile;
    plot(tFine, xFine, "-", "Color", [0.65 0.65 0.65], "LineWidth", 1.0);
    hold on;
    plot(t, x, "o-", "LineWidth", 1.2, "MarkerSize", 3);
    grid on;
    xlabel("Time [s]");
    ylabel("Amplitude");
    title(sprintf("Sampling (f0=%g Hz): fs=%g Hz, N=%d, Ts=%.3g ms", f0, fs, N, 1e3/fs));
end
exportgraphics(gcf, fullfile(figDir, "ch02_sampling_fs_compare.png"), "Resolution", 200, "BackgroundColor", "white");

%% Aliasing experiment: 折り返し周波数を「図＋数値」で確認する
f0 = 900;        % [Hz]
N = 4096;        % frequency resolution: Δf = fs/N
fsList = [2000, 1500];

figure;
tiledlayout(2, 2, "Padding", "compact", "TileSpacing", "compact");
for i = 1:numel(fsList)
    fs = fsList(i);
    [t, ~] = sp_time_axis(fs, N);
    x = sin(2*pi*f0*t);

    fAlias = abs(mod(f0 + fs/2, fs) - fs/2);

    % Time-domain (zoom in)
    nexttile;
    idx = (t <= 0.02);
    plot(t(idx), x(idx), "o-", "LineWidth", 1.2, "MarkerSize", 2.5);
    grid on;
    xlabel("Time [s]");
    ylabel("Amplitude");
    title(sprintf("Time waveform: fs=%g Hz (Nyquist=%g Hz)", fs, fs/2));

    % Spectrum (quick look)
    nexttile;
    sp_plot_mag_spec(fs, x, "NewFigure", false, "Db", false, ...
        "Title", sprintf("Spectrum: f0=%g Hz -> alias≈%g Hz", f0, fAlias));
    xlim([0 1000]);

    % Peak frequency (numerical check)
    X = fft(x);
    mag = abs(X) / N;
    kMax = floor(N/2);
    mag1 = mag(1:kMax+1);
    f = (0:kMax).' * (fs / N);
    [~, kHat] = max(mag1);
    fHat = f(kHat);

    fprintf("alias check: fs=%g Hz, predicted=%g Hz, measured peak=%g Hz\n", fs, fAlias, fHat);
    assert(abs(fHat - fAlias) < 2*(fs/N), "Aliasing peak mismatch (resolution too coarse?).");
end
exportgraphics(gcf, fullfile(figDir, "ch02_aliasing.png"), "Resolution", 200, "BackgroundColor", "white");

%% Uniform quantization: bits を変えたときの見え方と SNR
fs = 8000;
N = 4096;
f0 = 440;
[t, ~] = sp_time_axis(fs, N); %#ok<NASGU>
x = sin(2*pi*f0*(0:N-1).'/fs);

fullScale = 1;
bitsList = [4, 6, 8, 10, 12];
snrMeas = zeros(size(bitsList));
snrTheo = 6.02*bitsList + 1.76;

for i = 1:numel(bitsList)
    bits = bitsList(i);
    [xq, qerr] = sp_quantize_uniform(x, bits, "FullScale", fullScale);
    snrMeas(i) = 10*log10(mean(x.^2) / mean(qerr.^2));
end

figure;
plot(bitsList, snrMeas, "o-", "LineWidth", 1.2);
hold on;
plot(bitsList, snrTheo, "--", "LineWidth", 1.2);
grid on;
xlabel("Bits");
ylabel("SNR [dB]");
legend("measured", "rule-of-thumb (6.02B+1.76)", "Location", "northwest");
title("Quantization SNR vs. bits");
exportgraphics(gcf, fullfile(figDir, "ch02_quant_snr_vs_bits.png"), "Resolution", 200, "BackgroundColor", "white");

% Waveform example (low bits)
bits = 4;
[xq4, qerr4, step4] = sp_quantize_uniform(x, bits, "FullScale", fullScale); %#ok<ASGLU>
figure;
idx = (1:120).';
tt = (idx-1) / fs;
plot(tt, x(idx), "-", "Color", [0.65 0.65 0.65], "LineWidth", 1.2);
hold on;
stairs(tt, xq4(idx), "-", "LineWidth", 1.4);
plot(tt, xq4(idx), "o", "LineWidth", 1.0, "MarkerSize", 3.0);
grid on;
xlabel("Time [s]");
ylabel("Amplitude");
title(sprintf("Quantization example: bits=%d, step=%.4g", bits, step4));
legend("original", "quantized", "Location", "southwest");
exportgraphics(gcf, fullfile(figDir, "ch02_quant_waveform_bits4.png"), "Resolution", 200, "BackgroundColor", "white");

% Summary figure (for LaTeX): waveform + SNR vs bits
figure("Units", "pixels", "Position", [100 100 980 420]);
tiledlayout(1, 2, "Padding", "compact", "TileSpacing", "compact");

nexttile;
plot(tt, x(idx), "-", "Color", [0.65 0.65 0.65], "LineWidth", 1.2);
hold on;
stairs(tt, xq4(idx), "-", "LineWidth", 1.4);
plot(tt, xq4(idx), "o", "LineWidth", 1.0, "MarkerSize", 3.0);
grid on;
xlabel("Time [s]");
ylabel("Amplitude");
title(sprintf("Quantization example: bits=%d, step=%.4g", bits, step4));
legend("original", "quantized", "Location", "southwest");
hold off;

nexttile;
plot(bitsList, snrMeas, "o-", "LineWidth", 1.2); hold on;
plot(bitsList, snrTheo, "--", "LineWidth", 1.2);
grid on;
xlabel("Bits");
ylabel("SNR [dB]");
legend("measured", "rule-of-thumb (6.02B+1.76)", "Location", "northwest");
title("Quantization SNR vs. bits");
hold off;

exportgraphics(gcf, fullfile(figDir, "ch02_quant_summary.png"), "Resolution", 200, "BackgroundColor", "white");

%% Clipping: 飽和が作る高調波歪み
fs = 8000;
N = 4096;
f0 = 440;
[t, ~] = sp_time_axis(fs, N); %#ok<ASGLU>

amp = 1.3;
xIn = amp * sin(2*pi*f0*t);
y = sp_clip(xIn, -1, 1);

clipRatio = mean(abs(xIn) > 1);
fprintf("clipping ratio (|x|>1): %.3g\n", clipRatio);

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
idx = (t <= 0.02);
plot(t(idx), xIn(idx), "-", "LineWidth", 1.2);
hold on;
plot(t(idx), y(idx), "-", "LineWidth", 1.2);
grid on;
xlabel("Time [s]");
ylabel("Amplitude");
title(sprintf("Clipping (amp=%.2f): time waveform (zoom)", amp));
legend("input", "clipped", "Location", "southwest");

nexttile;
sp_plot_mag_spec(fs, xIn, "NewFigure", false, "Db", true, "Title", "input spectrum");
hold on;
sp_plot_mag_spec(fs, y, "NewFigure", false, "Db", true, "Title", "clipped spectrum");
title("Spectrum (dB): input vs clipped");
legend("input", "clipped", "Location", "southwest");

exportgraphics(gcf, fullfile(figDir, "ch02_clipping.png"), "Resolution", 200, "BackgroundColor", "white");
