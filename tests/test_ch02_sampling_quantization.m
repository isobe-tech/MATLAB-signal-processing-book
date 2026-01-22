%% tests/test_ch02_sampling_quantization.m
% 第2章のユーティリティ関数の簡易テスト（assertベース）

clear; clc;

thisFile = mfilename("fullpath");
bookRoot = fileparts(fileparts(thisFile));
addpath(fullfile(bookRoot, "src"));

% sp_clip
x = [-2; -1; -0.5; 0; 0.5; 1; 2];
y = sp_clip(x, -1, 1);
assert(isequal(y, [-1; -1; -0.5; 0; 0.5; 1; 1]));

% sp_quantize_uniform
x = [-1.2; -1; -0.3; 0; 0.3; 1; 1.2];
[xq, qerr, step, code] = sp_quantize_uniform(x, 4, "FullScale", 1);
assert(abs(step - 0.125) < 100*eps(step));
assert(isequal(xq, [-1; -1; -0.25; 0; 0.25; 1; 1]));
assert(max(abs((xq - x) - qerr)) == 0);
assert(all(code == round(code)));
assert(all(code >= -8) && all(code <= 8));

disp("OK: Chapter 2 tests passed.");

