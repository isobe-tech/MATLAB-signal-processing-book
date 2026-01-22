%% tests/test_ch03_complex_phase.m
% 第3章のユーティリティ関数の簡易テスト（assertベース）

clear; clc;

thisFile = mfilename("fullpath");
bookRoot = fileparts(fileparts(thisFile));
addpath(fullfile(bookRoot, "src"));

%% sp_complex_tone: basic properties
fs = 8000;
N = 16;
f0 = 0;
phi = 0.2*pi;
A = 1.5;
[x, t, n] = sp_complex_tone(fs, N, f0, "Phase", phi, "Amplitude", A); %#ok<ASGLU>

assert(iscolumn(x) && numel(x) == N);
assert(iscolumn(t) && numel(t) == N);
assert(iscolumn(n) && numel(n) == N);
assert(max(abs(x - A*exp(1j*phi))) == 0, "DC complex tone should be constant.");

%% sp_complex_tone: sample-to-sample ratio is constant
fs = 8000;
N = 128;
f0 = 440;
phi = -0.3*pi;
x = sp_complex_tone(fs, N, f0, "Phase", phi, "Amplitude", 1);
ratio = x(2:end) ./ x(1:end-1);
omega = 2*pi*f0/fs;
target = exp(1j*omega);
assert(max(abs(ratio - target)) < 2e-13, "Successive-sample ratio mismatch.");

%% Orthogonality on a DFT grid (finite-length)
N0 = 64;
fs0 = N0; % so f0=k gives exp(j*2*pi*k*n/N0)
tol = 1e-12;

x1 = sp_complex_tone(fs0, N0, 1);
x3 = sp_complex_tone(fs0, N0, 3);

ip13 = x1' * x3;
ip11 = x1' * x1;

assert(abs(ip13) < tol, "Inner product between different bins should be ~0.");
assert(abs(ip11 - N0) < tol, "Inner product with itself should be N.");

disp("OK: Chapter 3 tests passed.");

