%% tests/test_ch06_dft_fft.m
% 第6章のユーティリティ関数の簡易テスト（assertベース）

clear; clc;

thisFile = mfilename("fullpath");
bookRoot = fileparts(fileparts(thisFile));
addpath(fullfile(bookRoot, "src"));

fs = 8000;
tol = 1e-11;

%% sp_dft must match fft (real)
rng(0);
N = 64;
x = randn(N, 1);
X = sp_dft(x);
Xref = fft(x);
assert(norm(X - Xref) / max(norm(Xref), 1) < tol);

%% sp_dft must match fft (complex)
rng(1);
N = 63;
x = randn(N, 1) + 1j*randn(N, 1);
X = sp_dft(x);
Xref = fft(x);
assert(norm(X - Xref) / max(norm(Xref), 1) < tol);

%% sp_dft_matrix: W*x equals fft(x), and inverse recovers x
rng(2);
N = 32;
x = randn(N, 1) + 1j*randn(N, 1);
[W, Winv] = sp_dft_matrix(N);
X = W * x;
Xref = fft(x);
assert(norm(X - Xref) / max(norm(Xref), 1) < tol);

xRec = Winv * X;
assert(norm(xRec - x) / max(norm(x), 1) < tol);

%% sp_dft_matrix unitary form: W'*W = I
[Wu, ~] = sp_dft_matrix(N, "Unitary", true);
Ierr = norm(Wu' * Wu - eye(N), "fro") / N;
assert(Ierr < 1e-12);

%% sp_fft_axis: half range (one-sided)
Nfft = 8;
[fHalf, kHalf] = sp_fft_axis(fs, Nfft, "Range", "half");
assert(numel(fHalf) == floor(Nfft/2) + 1);
assert(kHalf(1) == 0 && kHalf(end) == floor(Nfft/2));
assert(abs(fHalf(1) - 0) < 1e-12);
assert(abs(fHalf(end) - fs/2) < 1e-12);

%% sp_fft_axis: whole range (fftshift-style)
Nfft = 8;
[fWhole, kWhole] = sp_fft_axis(fs, Nfft, "Range", "whole");
assert(numel(fWhole) == Nfft);
assert(kWhole(1) == -Nfft/2 && kWhole(end) == Nfft/2 - 1);
assert(abs(fWhole(1) + fs/2) < 1e-12);
assert(abs(fWhole(end) - (fs/2 - fs/Nfft)) < 1e-12);

disp("OK: Chapter 6 tests passed.");

