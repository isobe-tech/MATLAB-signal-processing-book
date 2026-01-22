%% tests/test_ch08_stft.m
% 第8章のユーティリティ関数の簡易テスト（assertベース）

clear; clc;

thisFile = mfilename("fullpath");
bookRoot = fileparts(fileparts(thisFile));
addpath(fullfile(bookRoot, "src"));

fs = 8000;

%% sp_stft: size and axis checks
rng(0);
N = 3000;
x = randn(N, 1);

L = 256;
hop = 128;
nfft = 512;

[t, f, P, info] = sp_stft(x, fs, ...
    "FrameLength", L, "Hop", hop, "Window", "hann", "Nfft", nfft, ...
    "Range", "half", "PadEnd", false, "TimeReference", "center", "Output", "psd");

nFramesExp = floor((N - L) / hop) + 1;
assert(numel(t) == nFramesExp);
assert(numel(f) == floor(nfft/2) + 1);
assert(isequal(size(P), [numel(f), numel(t)]));
assert(all(P(:) >= 0));

% Time axis reference
[tStart, ~, ~] = sp_stft(x, fs, ...
    "FrameLength", L, "Hop", hop, "Window", "hann", "Nfft", nfft, ...
    "Range", "half", "PadEnd", false, "TimeReference", "start", "Output", "psd");
assert(abs(tStart(1) - 0) < 1e-12);
assert(abs(t(1) - ((L-1)/2)/fs) < 1e-12);

% Output="complex" should return complex matrix
[~, ~, Xc] = sp_stft(x, fs, ...
    "FrameLength", L, "Hop", hop, "Window", "hann", "Nfft", nfft, ...
    "Range", "half", "PadEnd", false, "TimeReference", "center", "Output", "complex");
assert(~isreal(Xc));

%% STFT(PSD) time-average equals sp_psd_avg (same framing/window)
rng(1);
x = randn(fs * 4, 1);

L = 512;
hop = 256;
nfft = 1024;

[~, f, Pstft] = sp_stft(x, fs, ...
    "FrameLength", L, "Hop", hop, "Window", "hann", "Nfft", nfft, ...
    "Range", "half", "PadEnd", false, "TimeReference", "center", "Output", "psd");

Pavg = mean(Pstft, 2);
[f2, Pwelch] = sp_psd_avg(x, fs, ...
    "FrameLength", L, "Hop", hop, "Window", "hann", "Nfft", nfft, ...
    "Range", "half", "PadEnd", false);

assert(max(abs(f - f2)) < 1e-12);
relErr = norm(Pavg - Pwelch) / max(norm(Pwelch), 1e-12);
assert(relErr < 1e-12);

disp("OK: Chapter 8 tests passed.");





