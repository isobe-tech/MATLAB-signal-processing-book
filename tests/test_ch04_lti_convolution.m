%% tests/test_ch04_lti_convolution.m
% 第4章のユーティリティ関数の簡易テスト（assertベース）

clear; clc;

thisFile = mfilename("fullpath");
bookRoot = fileparts(fileparts(thisFile));
addpath(fullfile(bookRoot, "src"));

tol = 1e-12;

%% sp_convolve vs conv (full)
rng(0);
x = randn(64, 1);
h = randn(17, 1);
y = sp_convolve(x, h);
yRef = conv(x, h);
assert(norm(y - yRef) / max(norm(yRef), 1) < tol);

%% sp_convmtx_toeplitz: y = H*x matches conv(x,h)
rng(1);
Nx = 32;
h = randn(9, 1);
x = randn(Nx, 1);
H = sp_convmtx_toeplitz(h, Nx);
yMat = H * x;
yRef = conv(x, h);
assert(max(abs(yMat - yRef)) < tol);

%% sp_moving_average: causal slice matches conv(...)(1:N)
rng(2);
N = 128;
x = randn(N, 1);
M = 11;
[y, h] = sp_moving_average(x, M);
yRef = conv(x, h);
yRef = yRef(1:N);
assert(max(abs(y - yRef)) < tol);

%% sp_iir1: matches its own difference equation
rng(3);
N = 100;
x = randn(N, 1);
a = 0.8;
y = sp_iir1(x, a);
yRef = zeros(N, 1);
for n = 1:N
    if n == 1
        yRef(n) = (1-a) * x(n);
    else
        yRef(n) = (1-a) * x(n) + a * yRef(n-1);
    end
end
assert(max(abs(y - yRef)) < tol);

disp("OK: Chapter 4 tests passed.");

