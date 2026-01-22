%% tests/test_ch01_environment.m
% 第1章のユーティリティ関数の簡易テスト（assertベース）

clear; clc;

thisFile = mfilename("fullpath");
bookRoot = fileparts(fileparts(thisFile));
addpath(fullfile(bookRoot, "src"));

% sp_time_axis
fs = 8000;
N = 5;
[t, n] = sp_time_axis(fs, N);
assert(isequal(size(t), [N, 1]));
assert(isequal(size(n), [N, 1]));
assert(t(1) == 0);
assert(abs(t(2) - 1/fs) < 100*eps);

% sp_sine_with_noise: determinism with Seed
[x1, t1, xc1, w1, ns1] = sp_sine_with_noise(fs, 1024, 440, 20, "Seed", 123, "Amplitude", 2, "Phase", 0.1);
[x2, t2, xc2, w2, ns2] = sp_sine_with_noise(fs, 1024, 440, 20, "Seed", 123, "Amplitude", 2, "Phase", 0.1);
assert(isequal(x1, x2));
assert(isequal(t1, t2));
assert(isequal(xc1, xc2));
assert(isequal(w1, w2));
assert(ns1 == ns2);
assert(max(abs((xc1 + w1) - x1)) == 0);

sigPow = mean(abs(xc1).^2);
expectedNoisePow = sigPow / 10^(20/10);
expectedNs = sqrt(expectedNoisePow);
assert(abs(ns1 - expectedNs) < 100*eps(expectedNs));

% complex mode
[xc, ~, xcClean, wC, nsC] = sp_sine_with_noise(fs, 2048, 1000, 30, "Seed", 42, "Complex", true);
sigPowC = mean(abs(xcClean).^2);
expectedNoisePowC = sigPowC / 10^(30/10);
expectedNsC = sqrt(expectedNoisePowC/2);
assert(abs(nsC - expectedNsC) < 100*eps(expectedNsC));
assert(max(abs((xcClean + wC) - xc)) == 0);

disp("OK: Chapter 1 tests passed.");
