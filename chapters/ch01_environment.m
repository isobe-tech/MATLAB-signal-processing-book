%% 第1章 MATLABによる信号処理実験環境の構築
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
fs = 8000;    % [Hz]
N = 1024;     % number of samples
f0 = 440;     % [Hz]
snrDb = 10;   % [dB]

%% Experiment: sine + noise (reproducible)
rng(0); % 図を本文と一致させたいので固定する
[t, n] = sp_time_axis(fs, N); %#ok<NASGU>
[x, ~, xClean, w] = sp_sine_with_noise(fs, N, f0, snrDb); %#ok<ASGLU>

sp_plot_time(t, x, "Title", sprintf("Time waveform (fs=%g Hz, N=%d, f0=%g Hz, SNR=%g dB)", fs, N, f0, snrDb));
exportgraphics(gcf, fullfile(figDir, "ch01_time.png"), "Resolution", 200, "BackgroundColor", "white");

sp_plot_mag_spec(fs, x, "Title", "Magnitude spectrum (quick look)");
exportgraphics(gcf, fullfile(figDir, "ch01_mag.png"), "Resolution", 200, "BackgroundColor", "white");

%% Summary figure: time + spectrum (for LaTeX)
figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
sp_plot_time(t, x, ...
    "NewFigure", false, ...
    "Title", sprintf("Time waveform (fs=%g Hz, N=%d, f0=%g Hz, SNR=%g dB)", fs, N, f0, snrDb));
xlim([0, min(t(end), 0.03)]); % show the first ~30 ms for readability

nexttile;
sp_plot_mag_spec(fs, x, "NewFigure", false, "Title", "Magnitude spectrum (quick look)");
xline(f0, "k--", "f0");

exportgraphics(gcf, fullfile(figDir, "ch01_time_mag.png"), "Resolution", 200, "BackgroundColor", "white");

%% Basic statistics
m = mean(x);
sd = std(x);
rmsVal = sqrt(mean(abs(x).^2));
fprintf("mean=%.3g, std=%.3g, rms=%.3g\n", m, sd, rmsVal);

%% Checks: energy identity (two equivalent calculations)
E1 = x' * x;
E2 = sum(abs(x).^2);
relErr = abs(E1 - E2) / max(E2, 1);
fprintf("energy rel error = %.3g\n", relErr);
assert(relErr < 1e-12, "Energy mismatch is too large.");

%% Mini-check: transpose vs conjugate transpose
z = [1+2j; 3-4j];
zT = z.'; %#ok<NASGU>
zH = z';  %#ok<NASGU>
