%% tests/test_ch09_z_iir.m
% 第9章のユーティリティ関数の簡易テスト（assertベース）

clear; clc;

thisFile = mfilename("fullpath");
bookRoot = fileparts(fileparts(thisFile));
addpath(fullfile(bookRoot, "src"));

fs = 8000;

%% sp_iir2: matches MATLAB filter (double, DF-I)
rng(0);
N = 5000;
x = randn(N, 1);

theta = 0.25 * pi;
r = 0.9;
a = [1, -2*r*cos(theta), r^2];
b = [1, -2*cos(theta), 1]; % notch-like numerator (stable IIR)

yRef = filter(b, a, x); % allowed for verification/comparison
y = sp_iir2(x, b, a, "Structure", "df1", "Arithmetic", "double");

relErr = norm(y - yRef) / max(norm(yRef), 1e-12);
assert(relErr < 1e-12);

%% sp_iir2: DF-I and DF-II agree in double
y1 = sp_iir2(x, b, a, "Structure", "df1", "Arithmetic", "double");
y2 = sp_iir2(x, b, a, "Structure", "df2", "Arithmetic", "double");
assert(norm(y1 - y2) / max(norm(y1), 1e-12) < 1e-12);

%% Poles of a=[1, -2 r cos(theta), r^2] lie on radius r (conjugate pair)
p = roots(a);
assert(numel(p) == 2);
assert(max(abs(abs(p) - r)) < 1e-12);

%% 2nd-order all-pass: |H(e^{jΩ})| ~= 1 on the unit circle
thetaAp = 0.35 * pi;
rAp = 0.8;
aAp = [1, -2*rAp*cos(thetaAp), rAp^2];
bAp = [aAp(3), aAp(2), aAp(1)];

[f, H] = sp_freq_response(fs, bAp, "a", aAp, "N", 4096, "Range", "half");
assert(all(f >= 0) && abs(f(end) - fs/2) < 1e-12);
magErr = max(abs(abs(H) - 1));
assert(magErr < 1e-9);

disp("OK: Chapter 9 tests passed.");





