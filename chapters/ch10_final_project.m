%% 第10章 総合プロジェクト：帯域抽出（設計→評価→処理→説明）
% 配布コードは .m のみで完結させる。
% このファイルは「自分で設計を回す」ためのテンプレート。
%
% 目的：
% - 仕様（通過域/阻止域/リップル/減衰）を自分で言葉と数で定義する
% - 窓法でFIR係数を作る（ツールボックス無し）
% - 周波数応答を評価し，仕様を満たすか確認する
% - 信号へ適用し，スペクトル/PSD/スペクトログラムで改善を示す
% - 「なぜそのパラメータにしたか」を文章で説明できるようにする

clear; close all; clc;

% ルート（bookRoot）を推定して src/ を path に追加
thisFile = mfilename("fullpath");
bookRoot = fileparts(fileparts(thisFile));
addpath(fullfile(bookRoot, "src"));

%% 1) Signal to process (choose one)
% TODO: 実データがある場合は data/ から読み込む（今は合成信号で練習）
fs = 8000;   % [Hz]
T = 2.0;     % [s]
N = round(T * fs);
[t, ~] = sp_time_axis(fs, N);

% Example: wanted band + interference band + noise
rng(0);
xWanted = sin(2*pi*300*t) .* (t < 1.0) + 0.7*sin(2*pi*300*t) .* (t >= 1.0);
xInterf = 0.6*sin(2*pi*2200*t) + 0.25*randn(N, 1);
x = xWanted + xInterf;

%% 2) Spec definition (YOU decide)
% TODO: 「何を通す/落とすか」をHzで決める
Fp = 800;     % passband edge [Hz] (example)
Fs = 1200;    % stopband edge [Hz] (example)
Ntap = 101;   % taps (odd)
win = "hann"; % "rect"|"hann"|"hamming"

%% 3) Design by window method (no toolboxes)
[h, info] = sp_fir_window_design(fs, "lpf", Fp, Ntap, "Window", win, "Scale", true);
disp(info);

%% 4) Evaluate frequency response and specs
[f, H] = sp_freq_response(fs, h, "N", 4096, "Range", "half");
pb = [0, Fp];
sb = [Fs, fs/2];
spec = sp_filter_specs_eval(f, H, "Passband", pb, "Stopband", sb);
disp(spec);

figure;
plot(f, 20*log10(abs(H) + eps), "LineWidth", 1.2);
grid on;
xlabel("Frequency [Hz]");
ylabel("Magnitude [dB]");
title(sprintf("Designed FIR (window=%s, N=%d)", string(win), Ntap));
xlim([0, fs/2]);
ylim([-120, 5]);
xline(Fp, "k--", "Fp");
xline(Fs, "k--", "Fs");

%% 5) Apply filter
% TODO: convでも良いが、長い信号には filter が便利
y = filter(h, 1, x);

%% 5.1) Delay compensation and edge handling (for time-domain comparisons)
% Type-I linear-phase FIR (odd length) has group delay d=(Ntap-1)/2 samples.
% If you compare waveforms or compute time-domain errors/SNR, align first.
d = (Ntap - 1) / 2;
yA = y(d+1:end);
xA = x(1:end-d);

% Guard region to avoid transients (choose based on Ntap)
guard = Ntap;
if numel(xA) > 2*guard
    xE = xA(guard+1:end-guard);
    yE = yA(guard+1:end-guard);
else
    xE = xA;
    yE = yA;
end

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");
nexttile;
plot(t, x, "LineWidth", 1.0); grid on;
xlabel("Time [s]"); ylabel("Amplitude"); title("Input");
nexttile;
plot(t, y, "LineWidth", 1.0); grid on;
xlabel("Time [s]"); ylabel("Amplitude"); title("Output");

%% 6) Verification (PSD / spectrogram)
[fP, Pxx] = sp_psd_avg(x, fs, "FrameLength", 512, "Hop", 256, "Window", "hann", "Nfft", 2048, "Range", "half");
[~,  Pyy] = sp_psd_avg(y, fs, "FrameLength", 512, "Hop", 256, "Window", "hann", "Nfft", 2048, "Range", "half");

figure;
plot(fP, 10*log10(Pxx + 1e-20), "LineWidth", 1.1); hold on;
plot(fP, 10*log10(Pyy + 1e-20), "LineWidth", 1.1);
grid on;
xlabel("Frequency [Hz]");
ylabel("PSD [dB/Hz]");
title("PSD before/after filtering");
legend(["input", "output"], "Location", "best");
xlim([0, fs/2]);
hold off;

% Optional: STFT spectrogram (Chapter 8)
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

%% 6.1) Optional: SNR (only if you know the wanted component)
% For real data you usually don't know xWanted. For synthetic exercises you can.
if exist("xWanted", "var") && numel(xWanted) == numel(x)
    xWantedA = xWanted(1:end-d);
    if numel(xWantedA) > 2*guard
        xWantedE = xWantedA(guard+1:end-guard);
    else
        xWantedE = xWantedA;
    end
    nIn = xE - xWantedE;
    nOut = yE - xWantedE;
    snrInDb = 10*log10(sum(xWantedE.^2) / max(sum(nIn.^2), 1e-20));
    snrOutDb = 10*log10(sum(xWantedE.^2) / max(sum(nOut.^2), 1e-20));
    fprintf("SNR in/out: %.2f dB -> %.2f dB (improvement %.2f dB)\\n", snrInDb, snrOutDb, snrOutDb - snrInDb);
end

%% 7) Report (write in your own words)
% TODO:
% - なぜ Fp/Fs/Ntap/window を選んだか
% - 測定した ripple/atten は何dBで、仕様を満たすか
% - 図（周波数応答/PSD/スペクトログラム）が設計意図と整合しているか


