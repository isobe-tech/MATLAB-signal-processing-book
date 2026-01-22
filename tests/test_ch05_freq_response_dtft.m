%% tests/test_ch05_freq_response_dtft.m
% 第5章のユーティリティ関数の簡易テスト（assertベース）

clear; clc;

thisFile = mfilename("fullpath");
bookRoot = fileparts(fileparts(thisFile));
addpath(fullfile(bookRoot, "src"));

fs = 8000;
tol = 1e-12;

%% sp_freq_response (FIR): must match direct DTFT evaluation
rng(0);
b = randn(9, 1);
f = linspace(0, fs/2, 2000).';
[fOut, H] = sp_freq_response(fs, b, "F", f);

Omega = 2*pi*fOut/fs;
m = 0:numel(b)-1;
Href = exp(-1j * (Omega * m)) * b;
assert(norm(H - Href) / max(norm(Href), 1) < tol);

%% sp_freq_response (IIR): must match B/A on the unit circle
rng(1);
b = randn(5, 1);
a = [1; -0.3; 0.12]; % simple, stable-ish denominator for evaluation
f = linspace(0, fs/2, 2000).';
[fOut, H] = sp_freq_response(fs, b, "a", a, "F", f);

Omega = 2*pi*fOut/fs;
mb = 0:numel(b)-1;
ma = 0:numel(a)-1;
B = exp(-1j * (Omega * mb)) * b;
A = exp(-1j * (Omega * ma)) * a;
Href = B ./ A;
assert(norm(H - Href) / max(norm(Href), 1) < tol);

%% sp_group_delay: pure delay should be constant
d = 13;
hDelay = [zeros(d, 1); 1];
f = linspace(0, fs/2, 2000).';
[~, gd] = sp_group_delay(fs, hDelay, "F", f);

gdMid = gd(100:end-100);
assert(max(abs(gdMid - d)) < 1e-9);

%% sp_group_delay: symmetric FIR has (almost) constant delay where |H| is not tiny
M = 21;
h = ones(M, 1) / M;
[f, gd] = sp_group_delay(fs, h, "N", 4000);
[~, H] = sp_freq_response(fs, h, "F", f);
mask = abs(H) > 1e-2;
gdOk = gd(mask);
assert(max(abs(gdOk - (M-1)/2)) < 1e-2);

disp("OK: Chapter 5 tests passed.");

