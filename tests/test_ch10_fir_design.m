%% tests/test_ch10_fir_design.m
% 第10章のユーティリティ関数の簡易テスト（assertベース）

clear; clc;

thisFile = mfilename("fullpath");
bookRoot = fileparts(fileparts(thisFile));
addpath(fullfile(bookRoot, "src"));

fs = 8000;

%% sp_fir_ideal_lpf: symmetry around center (odd N)
N = 101;
fc = 1000;
h = sp_fir_ideal_lpf(fs, fc, N);
assert(isequal(size(h), [N, 1]));
assert(max(abs(h - flipud(h))) < 1e-12);

%% sp_fir_window_design: LPF properties
[h, info] = sp_fir_window_design(fs, "lpf", fc, N, "Window", "hann", "Scale", true);
assert(isfield(info, "groupDelaySamples"));
assert(abs(info.groupDelaySamples - (N-1)/2) < 1e-12);
assert(max(abs(h - flipud(h))) < 1e-12); % linear phase (Type-I)

% DC gain ~ 1 (after scaling)
g0 = sum(h);
assert(abs(g0 - 1) < 1e-10);

%% FIR frequency response: sp_freq_response equals FFT(h) on the same grid
Nfft = 4096;
Hfft = fft(h, Nfft);
kMax = floor(Nfft/2);
f = (0:kMax).' * (fs/Nfft);
Hfft = Hfft(1:kMax+1);

[~, Hdir] = sp_freq_response(fs, h, "F", f, "Range", "half");
relErr = norm(abs(Hdir) - abs(Hfft)) / max(norm(abs(Hdir)), 1e-12);
assert(relErr < 1e-10);

%% Spec evaluation: stopband attenuation should be "decent" for a Hann window
pb = [0, fc];
sb = [fc + 400, fs/2];
spec = sp_filter_specs_eval(f, Hdir, "Passband", pb, "Stopband", sb);
assert(isfield(spec, "stopbandAttenDb"));
assert(spec.stopbandAttenDb > 20); % loose lower bound (depends on N/spec)

%% sp_fir_ls_lpf: basic properties (symmetry, DC gain)
Fp = 800;
Fs = 1200;
[hLs, infoLs] = sp_fir_ls_lpf(fs, Fp, Fs, N, "GridN", 4096, "Wp", 1, "Ws", 80, "Scale", true);
assert(isfield(infoLs, "method") && infoLs.method == "ls");
assert(max(abs(hLs - flipud(hLs))) < 1e-10);
assert(abs(sum(hLs) - 1) < 1e-10);

% Basic spec check (should not be terrible)
[~, Hls] = sp_freq_response(fs, hLs, "F", f, "Range", "half");
specLs = sp_filter_specs_eval(f, Hls, "Passband", [0, Fp], "Stopband", [Fs, fs/2]);
assert(isfield(specLs, "stopbandAttenDb"));
assert(specLs.stopbandAttenDb > 15);

disp("OK: Chapter 10 tests passed.");


