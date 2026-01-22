%% tests/test_ch07_window_psd.m
% 第7章のユーティリティ関数の簡易テスト（assertベース）

clear; clc;

thisFile = mfilename("fullpath");
bookRoot = fileparts(fileparts(thisFile));
addpath(fullfile(bookRoot, "src"));

fs = 8000;

%% Window functions: shape and basic properties
N = 16;
wRect = sp_window_rect(N);
wHann = sp_window_hann(N);
wHamm = sp_window_hamming(N);

assert(isequal(size(wRect), [N, 1]));
assert(isequal(size(wHann), [N, 1]));
assert(isequal(size(wHamm), [N, 1]));

assert(all(wRect == 1));
assert(max(abs(wHann - flipud(wHann))) < 1e-12);  % symmetric
assert(max(abs(wHamm - flipud(wHamm))) < 1e-12);  % symmetric

% Hann endpoints are (almost) zero for N>1
assert(abs(wHann(1)) < 1e-12 && abs(wHann(end)) < 1e-12);

%% sp_frame_signal: framing without padding
x = (1:10).';
[X, startIdx] = sp_frame_signal(x, 4, 2, "PadEnd", false);
assert(isequal(startIdx, [1; 3; 5; 7]));
assert(isequal(size(X), [4, 4]));
assert(isequal(X(:, 1), [1; 2; 3; 4]));
assert(isequal(X(:, 4), [7; 8; 9; 10]));

%% sp_frame_signal: framing with padding
x = (1:9).';
[X, startIdx] = sp_frame_signal(x, 4, 3, "PadEnd", true, "PadValue", 0);
assert(isequal(startIdx, [1; 4; 7]));
assert(isequal(size(X), [4, 3]));
assert(isequal(X(:, 3), [7; 8; 9; 0]));

%% sp_psd_avg: power check (frequency integral ~= time-domain mean-square)
rng(0);
N = 4096;
x = randn(N, 1);

[f, Pxx, info] = sp_psd_avg(x, fs, ...
    "FrameLength", 256, "Hop", 128, "Window", "hann", "Nfft", 256, "Range", "half");

assert(isfield(info, "df"));
assert(numel(f) == floor(info.nfft/2) + 1);
assert(all(f >= 0) && abs(f(end) - fs/2) < 1e-12);

pTime = mean(x.^2);
pFreq = sum(Pxx) * info.df;
relErr = abs(pFreq - pTime) / max(pTime, 1e-12);
assert(relErr < 0.05);

%% sp_psd_avg: specifying window vector should match specifying by name (within tiny numerical tolerance)
w = sp_window_hann(256);
[~, Pname] = sp_psd_avg(x, fs, "FrameLength", 256, "Hop", 128, "Window", "hann", "Nfft", 256, "Range", "half");
[~, Pvec]  = sp_psd_avg(x, fs, "FrameLength", 256, "Hop", 128, "Window", w,      "Nfft", 256, "Range", "half");
assert(norm(Pname - Pvec) / max(norm(Pname), 1e-12) < 1e-12);

disp("OK: Chapter 7 tests passed.");





