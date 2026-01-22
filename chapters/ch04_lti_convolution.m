%% 第4章 LTIシステムと畳み込み：FIR/IIRフィルタの基本構造
% このスクリプトは「Run All」で図まで再現できる形にする。

clear; close all; clc;

% ルート（bookRoot）を推定して src/ を path に追加
thisFile = mfilename("fullpath");
bookRoot = fileparts(fileparts(thisFile));
addpath(fullfile(bookRoot, "src"));
sp_book_style(); % White background + readable fonts for LaTeX

% 出力先（図）
figDir = fullfile(bookRoot, "figs");
if ~exist(figDir, "dir")
    mkdir(figDir);
end

%% Parameters
fs = 8000;     % [Hz]
N  = 1600;     % samples (0.2 s)
[t, n] = sp_time_axis(fs, N); %#ok<ASGLU>

%% 0) Quick checks: sp_convolve must match conv (full)
rng(0);
xTest = randn(64, 1);
hTest = randn(9, 1);
y0 = sp_convolve(xTest, hTest);
yRef = conv(xTest, hTest);
relErr = norm(y0 - yRef) / max(norm(yRef), 1);
fprintf("sp_convolve vs conv: relErr=%.3g\n", relErr);
assert(relErr < 1e-12, "sp_convolve mismatch.");

%% 1) Changing impulse response changes the output (same input, different h)
[x, ~] = sp_sine_with_noise(fs, N, 440, 10, "Amplitude", 1, "Phase", 0.1*pi, "Seed", 0); %#ok<ASGLU>

Mma = 31;
hMA = ones(Mma, 1) / Mma;
hDiff = [1; -1]; % simple difference

yMAFull = sp_convolve(x, hMA);
yDiffFull = sp_convolve(x, hDiff);
yMA = yMAFull(1:N);
yDiff = yDiffFull(1:N);

idx = (t <= 0.03); % show first 30 ms

figure;
tiledlayout(2, 2, "Padding", "compact", "TileSpacing", "compact");

nexttile;
plot(t(idx), x(idx), "LineWidth", 1.1);
grid on;
xlabel("Time [s]"); ylabel("x[n]");
title("Input (sine + noise)");

nexttile;
hold on;
stem(0:numel(hMA)-1, hMA, "filled", "MarkerSize", 2);
stem(0:numel(hDiff)-1, hDiff, "filled", "MarkerSize", 2);
grid on;
xlabel("m"); ylabel("h[m]");
legend("Moving average", "Difference", "Location", "northeast");
title("Impulse responses");

nexttile;
plot(t(idx), yMA(idx), "LineWidth", 1.1);
grid on;
xlabel("Time [s]"); ylabel("y[n]");
title(sprintf("Output: moving average (M=%d)", Mma));

nexttile;
plot(t(idx), yDiff(idx), "LineWidth", 1.1);
grid on;
xlabel("Time [s]"); ylabel("y[n]");
title("Output: difference (high-pass)");

exportgraphics(gcf, fullfile(figDir, "ch04_impulse_response_effect.png"), "Resolution", 200, "BackgroundColor", "white");

%% 2) Boundary effects: zero-padding shows up at the ends
M = 51;
h = ones(M, 1) / M;

xConst = ones(N, 1);
yFull = sp_convolve(xConst, h);
tFull = (0:numel(yFull)-1).' / fs;

T = (N-1) / fs;

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
plot(t, xConst, "--", "LineWidth", 1.0);
hold on;
plot(t, yFull(1:N), "-", "LineWidth", 1.2);
grid on;
xlabel("Time [s]"); ylabel("Amplitude");
title("Start transient (zero-padding before n=0)");
legend("x (constant)", "y (causal slice)", "Location", "southeast");
xlim([0, 0.03]);

nexttile;
plot(tFull, yFull, "-", "LineWidth", 1.2);
hold on;
xline(T, "k--", "LineWidth", 1.0);
grid on;
xlabel("Time [s]"); ylabel("Amplitude");
title("End transient in full convolution (zero-padding after the data)");
xlim([max(T-0.03, 0), tFull(end)]);

exportgraphics(gcf, fullfile(figDir, "ch04_boundary_effect.png"), "Resolution", 200, "BackgroundColor", "white");

%% Common input: unit step (used in multiple comparisons)
t0 = 0.04;
xStep = double(t >= t0);

%% 2b) Convolution output modes: full vs causal slice vs same (and filter)
Mmode = 31;
hMode = ones(Mmode, 1) / Mmode;

yFull = conv(xStep, hMode);
yCausal = yFull(1:N);
ySame = conv(xStep, hMode, "same");

yFilter = filter(hMode, 1, xStep);
relErr = norm(yCausal - yFilter) / max(norm(yFilter), 1);
fprintf("causal slice vs filter: relErr=%.3g\n", relErr);
assert(relErr < 1e-12, "filter(h,1,x) mismatch.");

idx = (t >= 0.02) & (t <= 0.12);

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

	nexttile;
	plot(t(idx), xStep(idx), "k--", "LineWidth", 1.0);
	hold on;
	plot(t(idx), yCausal(idx), "LineWidth", 1.2);
	plot(t(idx), ySame(idx), "LineWidth", 1.2);
	grid on;
	xlabel("Time [s]"); ylabel("Amplitude");
	title(sprintf("conv output modes (step input, moving average M=%d)", Mmode));
	lgd = legend(["x (step)", "yFull(1:N) (causal slice)", "conv(x,h,'same')"], "Location", "southeast");
	lgd.AutoUpdate = "off";
	xline(t0, "k:", "LineWidth", 1.0);

nexttile;
plot(t(idx), ySame(idx) - yCausal(idx), "LineWidth", 1.2);
grid on;
xlabel("Time [s]"); ylabel("difference");
title("Difference: same - causal slice (definition, not a bug)");
xline(t0, "k:", "LineWidth", 1.0);

exportgraphics(gcf, fullfile(figDir, "ch04_conv_modes.png"), "Resolution", 200, "BackgroundColor", "white");

%% 3) FIR example: moving average step response (smooth + delay)
Mlist = [5, 21, 51];
figure;
hold on;
plot(t, xStep, "k--", "LineWidth", 1.0);
for k = 1:numel(Mlist)
    Mk = Mlist(k);
    yk = sp_moving_average(xStep, Mk);
    plot(t, yk, "LineWidth", 1.2);
end
grid on;
xlabel("Time [s]"); ylabel("Amplitude");
title("Moving average: step response (causal output, length N)");
legend(["x (step)", "M=5", "M=21", "M=51"], "Location", "southeast");
xlim([0, 0.12]);

exportgraphics(gcf, fullfile(figDir, "ch04_moving_average_step.png"), "Resolution", 200, "BackgroundColor", "white");

%% 3b) Moving average on noisy sine: noise reduction vs delay
[xNoisy, ~, xClean] = sp_sine_with_noise(fs, N, 440, 0, "Amplitude", 1, "Phase", 0.1*pi, "Seed", 1);

Mlist = [5, 21, 51];
idx = (t >= 0.03) & (t <= 0.06); % zoom window

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
plot(t(idx), xNoisy(idx), "Color", [0.7, 0.7, 0.7], "LineWidth", 0.9);
hold on;
plot(t(idx), xClean(idx), "k--", "LineWidth", 1.1);
grid on;
xlabel("Time [s]"); ylabel("Amplitude");
title("Noisy input vs clean sine (zoom)");
legend("x (noisy)", "xClean", "Location", "southeast");

nexttile;
plot(t(idx), xClean(idx), "k--", "LineWidth", 1.1);
hold on;
for k = 1:numel(Mlist)
    Mk = Mlist(k);
    yk = sp_moving_average(xNoisy, Mk);
    plot(t(idx), yk(idx), "LineWidth", 1.1);
end
grid on;
xlabel("Time [s]"); ylabel("Amplitude");
title("Moving average outputs (noise is reduced, but delay increases)");
legend(["xClean", "M=5", "M=21", "M=51"], "Location", "southeast");

exportgraphics(gcf, fullfile(figDir, "ch04_moving_average_noise.png"), "Resolution", 200, "BackgroundColor", "white");

%% 3c) Moving average frequency response preview (magnitude)
Mlist = [5, 21, 51];
f = linspace(0, fs/2, 2000);
Omega = 2*pi*f/fs;

figure;
hold on;
for k = 1:numel(Mlist)
    Mk = Mlist(k);
    m = 0:(Mk-1);
    H = exp(-1j * (Omega(:) * m)) * (ones(Mk, 1) / Mk);
    magDb = 20*log10(abs(H) + 1e-12);
    plot(f, magDb, "LineWidth", 1.2);
end
grid on;
xlabel("Frequency [Hz]"); ylabel("|H| [dB]");
title("Moving average frequency response (preview)");
legend(["M=5", "M=21", "M=51"], "Location", "northeast");
xlim([0, fs/2]);
ylim([-80, 5]);

exportgraphics(gcf, fullfile(figDir, "ch04_moving_average_freqresp.png"), "Resolution", 200, "BackgroundColor", "white");

%% 4) IIR example: first-order step response (time constant intuition)
alist = [0.2, 0.8, 0.95];
figure;
hold on;
plot(t, xStep, "k--", "LineWidth", 1.0);
for k = 1:numel(alist)
    a = alist(k);
    y = sp_iir1(xStep, a);
    plot(t, y, "LineWidth", 1.2);
end
grid on;
xlabel("Time [s]"); ylabel("Amplitude");
title("First-order IIR: step response  y[n]=(1-a)x[n]+ay[n-1]");
legend(["x (step)", "a=0.2", "a=0.8", "a=0.95"], "Location", "southeast");
xlim([0, 0.12]);

exportgraphics(gcf, fullfile(figDir, "ch04_iir1_step.png"), "Resolution", 200, "BackgroundColor", "white");

%% 4b) IIR impulse response: "Infinite" tail (exponential decay)
K = 120; % show first K samples
xImp = zeros(N, 1);
xImp(1) = 1;

alist = [0.2, 0.8, 0.95];

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
hold on;
for k = 1:numel(alist)
    a = alist(k);
    h = sp_iir1(xImp, a);
    stem(0:K-1, h(1:K), "filled", "MarkerSize", 2);
end
grid on;
xlabel("n"); ylabel("h[n]");
title("First-order IIR impulse response  h[n]=(1-a)a^n  (zero initial condition)");
legend(["a=0.2", "a=0.8", "a=0.95"], "Location", "northeast");

nexttile;
hold on;
for k = 1:numel(alist)
    a = alist(k);
    h = sp_iir1(xImp, a);
    s = cumsum(h);
    plot(0:K-1, s(1:K), "LineWidth", 1.1);
end
grid on;
xlabel("n"); ylabel("sum_{k=0}^{n} h[k]");
title("Cumulative sum approaches 1 (step response connection)");
legend(["a=0.2", "a=0.8", "a=0.95"], "Location", "southeast");
ylim([-0.05, 1.05]);

exportgraphics(gcf, fullfile(figDir, "ch04_iir1_impulse.png"), "Resolution", 200, "BackgroundColor", "white");

%% 5) Toeplitz convolution matrix: visualization
Nx = 18;
M = 6;
hT = ones(M, 1) / M;
H = sp_convmtx_toeplitz(hT, Nx);

figure;
imagesc(H);
axis image;
colormap(gray);
cb = colorbar;
try
    cb.Color = [0 0 0];
catch
end
xlabel("column (x index)");
ylabel("row (y index)");
title("Toeplitz convolution matrix H (full, FIR)");

exportgraphics(gcf, fullfile(figDir, "ch04_toeplitz_conv_matrix.png"), "Resolution", 200, "BackgroundColor", "white");

%% 6) y = H*x equals conv(h,x)
xSmall = zeros(Nx, 1);
xSmall(3) = 1;
xSmall(8) = -0.6;
xSmall(14) = 0.4;

yMat = H * xSmall;
yConv = conv(xSmall, hT);
maxErr = max(abs(yMat - yConv));
fprintf("Toeplitz vs conv: maxErr=%.3g\n", maxErr);
assert(maxErr < 1e-12, "Toeplitz matrix mismatch.");

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
stem(0:numel(yConv)-1, yConv, "filled", "MarkerSize", 2);
hold on;
plot(0:numel(yMat)-1, yMat, "o", "MarkerSize", 3);
grid on;
xlabel("n"); ylabel("y[n]");
legend("conv(x,h)", "H*x", "Location", "northeast");
title("Full convolution: conv(x,h) vs H*x");

nexttile;
stem(0:numel(yConv)-1, yMat - yConv, "filled", "MarkerSize", 2);
grid on;
xlabel("n"); ylabel("error");
title("Difference (should be ~0)");

exportgraphics(gcf, fullfile(figDir, "ch04_toeplitz_vs_conv.png"), "Resolution", 200, "BackgroundColor", "white");
